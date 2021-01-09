using VrpSolver, JuMP, ArgParse
# using CPLEX
include("data.jl")
include("model.jl")
include("solution.jl")
include("SparseMaxFlowMinCut.jl")

function parse_commandline(args_array::Array{String,1}, appfolder::String)
    s = ArgParseSettings(usage="##### VRPSolver #####\n\n" *
	   "  On interactive mode, call main([\"arg1\", ..., \"argn\"])", exit_after_help=false)
    @add_arg_table s begin
        "instance"
            help = "Instance file path"
        "--preprocessings"
            help = "Instance edges preprocessings"
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
        "--instance_type", "-i"
            help = "Select the instance type (EMH, MATHEUS)"
            default = "EMH"
    end
   return parse_args(args_array, s)
end

function run_gvrp(app::Dict{String,Any})
    println("Application parameters:")
    for (arg, val) in app
        println("  $arg  =>  $(repr(val))")
    end
    flush(stdout)

    instance_name = split(basename(app["instance"]), ".")[1] 
    if app["instance_type"] == "Matheus"
      data = readMatheusInstance(app)
    else
      data = readEMHInstance(app)
    end

    if app["sol"] != nothing
        sol = readsolution(app)
        checksolution(data, sol, app) # checks the solution feasibility
        app["ub"] = (sol.cost < app["ub"]) ? sol.cost : app["ub"] # update the upper bound if necessary
    end

    solution_found = false
    if !app["nosolve"]
        (model, x) = build_model(data)

        #println(model.formulation)
        # enum_paths, complete_form = get_complete_formulation(model, app["cfg"])
        # complete_form.solver = CplexSolver() # set MIP solver
        # print_enum_paths(enum_paths)
        # println(complete_form)
        # solve(complete_form)
        # println("Objective value: $(getobjectivevalue(complete_form))\n")
        
        optimizer = VrpOptimizer(model, app["cfg"], instance_name)
        set_cutoff!(optimizer, app["ub"])
        
        (status, solution_found) = optimize!(optimizer)
        if solution_found
            sol = getsolution(data, optimizer, x, get_objective_value(optimizer), app)
        end
    end

    println("########################################################")
    if solution_found || app["sol"] != nothing # Is there a solution?
        print_routes(sol)
        #checksolution(data, sol)
        println("Cost $(sol.cost)")
        checksolution(data, sol)
        if app["out"] != nothing
            writesolution(app["out"], sol)
        end
        if app["tikz"] != nothing
            if data.coord
                drawsolution(app["tikz"], data, sol) # write tikz figure
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
