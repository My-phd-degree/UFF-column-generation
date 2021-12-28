using VrpSolver, JuMP, ArgParse

include("data.jl")
include("gvrp_afs_tree.jl")
include("reduced_graph.jl")
include("lower_bounds.jl")
include("model.jl")
include("model_compact_with_arcs.jl")
include("model_y.jl")
include("model_compact_y.jl")
include("solution.jl")
include("preprocessings.jl")
include("SparseMaxFlowMinCut.jl")

function parse_commandline(args_array::Array{String,1}, appfolder::String)
    s = ArgParseSettings(usage="##### VRPSolver #####\n\n" *
	   "  On interactive mode, call main([\"arg1\", ..., \"argn\"])", exit_after_help=false)
    @add_arg_table s begin
        "instance"
            help = "Instance file path"
        "--model-type"
        help = "Model type (normal, compacted-with-arcs, y (only for non-consec instances), compacted-y)"
            default = "normal"
        "--preprocessings"
            help = "Instance edges preprocessings"
        "--non-consec" 
            help = "If the instance must be solved without using edges between AFSs"
            action = :store_true
        "--cfg", "-c"
            help = "Configuration file path"
            default = "$appfolder/../config/BWTSP.cfg"
        "--ub", "-u"
            help = "Upper bound (primal bound)"
            arg_type = Float64
            default = 10000000.0
        "--noround", "-r"
            help = "Does not round the distance matrix"
            action = :store_true
        "--sol", "-s"
            help = "Solution file path."
        "--out", "-o"
            help = "Path to write the solution found"
        "--tikz", "-t"
            help = "Path to write the TikZ figure of the solution found."
        "--nosolve", "-n"
            help = "Does not call the VRPSolver. Only to check or draw a given solution."
            action = :store_true
        "--batch", "-b"
            help = "batch file path"
        "--instance-type", "-i"
            help = "Select the instance type (EMH, Matheus, Andelmin-Bartolini)"
            default = "EMH"
        "--verbose", "-v"
            help = "true to see the logs, false otherwise"
            default = true
          end
   return parse_args(args_array, s)
end

function run_gvrp(app::Dict{String,Any})
  if app["verbose"]
    println("Application parameters:")
    for (arg, val) in app
      println("  $arg  =>  $(repr(val))")
    end
  end
  flush(stdout)

  instance_name = split(basename(app["instance"]), ".")[1] 
  if app["instance-type"] == "Matheus"
    data = readMatheusInstance(app)
  elseif app["instance-type"] == "EMH"
    data = readEMHInstance(app)
  elseif app["instance-type"] == "Andelmin-Bartolini"
    data = read_Andelmin_Bartolini_Instance(app)
  end

  if app["sol"] != nothing
    sol = read_Andelmin_Bartolini_Solution(app, data)
    checksolution(data, sol) # checks the solution feasibility
    app["ub"] = (sol.cost < app["ub"]) ? sol.cost : app["ub"] # update the upper bound if necessary
  end

  solution_found = false
  if !app["nosolve"]
    if app["model-type"] == "compacted-with-arcs"
      data.non_consec && error("The model compacted-with-arcs only can be executed without non_consec flag")
      (directedData, model, x, afss_pairs) = build_model_compact_with_arcs(data)
      optimizer = VrpOptimizer(model, app["cfg"], instance_name)
      set_cutoff!(optimizer, app["ub"])
      (status, solution_found) = optimize!(optimizer)
      if solution_found
        sol = getsolution_compact_with_arcs(data, directedData, optimizer, x, get_objective_value(optimizer), app, afss_pairs)
      end
    elseif app["model-type"] == "compacted-y"
      data.non_consec && error("The model compacted-y only can be executed without non_consec flag")
      (model, P, x, y) = build_model_compact_y(data)
      optimizer = VrpOptimizer(model, app["cfg"], instance_name)
      set_cutoff!(optimizer, app["ub"])
      (status, solution_found) = optimize!(optimizer)
      if solution_found
        sol = getsolution_compact_y(data, optimizer, P, x, y, get_objective_value(optimizer), app)
      end
    elseif app["model-type"] == "y"
      !data.non_consec && error("The model y only can be executed with non_consec flag")
      (model, x, y) = build_model_y(data)
      optimizer = VrpOptimizer(model, app["cfg"], instance_name)
      set_cutoff!(optimizer, app["ub"])
      (status, solution_found) = optimize!(optimizer)
      if solution_found
        sol = getsolution_y(data, optimizer, x, y, get_objective_value(optimizer), app)
      end
    else
      (model, x, y) = build_model(data)
      optimizer = VrpOptimizer(model, app["cfg"], instance_name)
      set_cutoff!(optimizer, app["ub"])
      (status, solution_found) = optimize!(optimizer)
      if solution_found
        sol = getsolution(data, optimizer, x, y, get_objective_value(optimizer), app)
      end
    end
  end
  println("########################################################")
  if solution_found || app["sol"] != nothing # Is there a solution?
    print_routes(data, sol)
    checksolution(data, sol)
    println("Cost $(sol.cost)")
    #checksolution(data, sol)
    if app["out"] != nothing
      writesolution(app["out"], data, sol)
    end
    if app["tikz"] != nothing
      if data.coord
        #drawsolution(app["tikz"], data, sol) # write tikz figure
        drawsolution(app["tikz"], data, sol.routes) # write tikz figure
      else
        println("TikZ figure ($(app["tikz"])) will not be generated, since the instance has no coordinates.")
      end
    end
  elseif !app["nosolve"]
    if status == :Optimal
      println("Problem infeasible")
    else
      println("Solution not found")
    end
  end
  println("########################################################")
end

function main(ARGS)
    appfolder = dirname(@__FILE__)
    app = parse_commandline(ARGS, appfolder)
    isnothing(app) && return
    if app["batch"] != nothing
        for line in readlines(app["batch"])
            if isempty(strip(line)) || strip(line)[1] == '#'
                continue
            end
            args_array = [String(s) for s in split(line)]
            app_line = parse_commandline(args_array, appfolder)
            run_gvrp(app_line)
        end
    else
        run_gvrp(app)
    end
end

# main()
if isempty(ARGS)
    main(["--help"])
else
    main(ARGS)
end
