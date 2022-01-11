include("gvrp_afs_tree.jl")

using CPLEX

ed(i, j) = i < j ? (i, j) : (j, i)

function build_model(data::DataGVRP)
  n_2_path_cuts_gh = 0
  n_k_path_cuts_time = 0
  n_energy_cuts = 0
  E = edges(data) # set of edges of the input graph G′
  n = nb_vertices(data)
  V = [i for i in 1:n] # set of vertices of the input graph G′
  V′ = [i for i in 0:n] # V ⋃ {0}, where 0 is a dummy vertex
  depot_id = data.depot_id

  β = data.β
  C = data.C # Set of customers vertices
  C₀ = vcat([data.depot_id], data.C) # Set of customers vertices with depot
  F = data.F # Set of AFSs vertices
  T = data.T # Set of AFSs vertices
  F₀ = vcat([data.depot_id], data.F) # Set of AFSs vertices with depot

  #println("Customers ")
  #for i in data.C
  #  println(i, " ", data.G′.V′[i], ", δ($i) = $(δ(data, i))")
  #end
  #println("AFSs ")
  #for i in data.F
  #  println(i, " ", data.G′.V′[i], ", δ($i) = $(δ(data, i))")
  #end
  # println("F´: ", F´)
  # println("β: $β")
  # println("T: $T")
  # println("ρ: $(data.ρ)")
  # println("ε: $(data.ε)")
  # println("E: ", length(E))
  # for (i, j) in E
  #   println("d(($i, $j)) = ", d(data, (i, j)), " f(($i, $j)) = ", f(data, (i, j)), " t(($i, $j)) = ", t(data, (i, j)))
  # end

  # Formulation
  gvrp = VrpModel()
  @variable(gvrp.formulation, 2 * length(C) >= x[e in E] >= 0, Int)
  @variable(gvrp.formulation, 2 * length(C) >= y[i in F₀] >= 0, Int)
  @objective(gvrp.formulation, Min, sum((d(data, (i, j)) * x[(i, j)] for (i, j) in E)))
  @constraint(gvrp.formulation, deg[i in C], sum(x[e] for e in δ(data, i)) == 2.0)
  @constraint(gvrp.formulation, hotel_deg[i in F₀], sum(x[e] for e in δ(data, i)) == 2*y[i])
  if data.non_consec
    @constraint(gvrp.formulation, no_edge_between_afss[f in F], sum(x[(f, r)] for r in F if (f, r) in E) == 0.0)
  end
  @constraint(gvrp.formulation, y[F₀[1]] >= 1)

  invalidEdges = vcat(get_invalid_edges_1(data), get_invalid_edges_2(data), get_invalid_edges_3(data), get_invalid_edges_4(data))
  [@constraint(gvrp.formulation, x[e] == 0) for e in invalidEdges]


#=
  routes = [
            [0, 6 , 30, 18, 15, 46, 25, 24],
            [0, 43, 3 , 34, 14, 13, 5],
            [0, 2 , 44, 39, 40, 2],
            [0, 9 , 32, 12, 16, 17, 41, 8, 37],
            [0, 13, 48, 36, 7 , 47, 35],
            [0, 19, 50, 22, 21, 45, 23, 29, 4],
            [0, 31, 1 , 27, 13],
            [0, 11, 42, 33, 21, 20, 38],
            [0, 49, 26, 28, 10, 3],
           ]
=#

  if @isdefined routes
    ids = Dict{Int,Int}()
    ids_ = Dict{Int,Int}()
    for route in routes
      for id in route
        if !haskey(ids, id)
          found = false
          for i in 1:length(data.G′.V′)
            v = data.G′.V′[i]
            if v.id_vertex == id + 1
              found = true
              ids[id] = i 
              ids_[i] = id
              break
            end
          end
          !found && error("Node of id $id wasn't found")
        end
      end
    end
    edgesWeights = Dict{Tuple{Int64, Int64},Int64}()
    for route in routes
      println([ids[route[i]] for i in 1:length(route)])
      for i in 2:length(route)
        e = ed(ids[route[i]], ids[route[i - 1]])
        if haskey(edgesWeights, e)
          edgesWeights[e] = edgesWeights[e] + 1
        else
          edgesWeights[e] = 1
        end
      end 
    end
    for e in keys(edgesWeights)
      println(f(data, e))
      println(data.G′.V′[e[1]], " ", data.G′.V′[e[2]])
      println("Edge $e; Edge weight $(edgesWeights[e]); EdgeId $(ed(ids_[e[1]], ids_[e[2]]))")
      @constraint(gvrp.formulation, x[e] == edgesWeights[e])
    end
  end

  # Build the model directed graph G=(V,A)
  function build_graph()

    v_source = v_sink = 1
    L, U = 0, length(C) # max and min number of paths is equal to number of AFSs

    V_ = deepcopy(V)

    #AFSs vertex split
    new_id = length(V_) + 1
    afssIn = Dict{Int, Int}()
    afssOut = Dict{Int, Int}()
    for f in F
      push!(V_, new_id)
      afssIn[f] = new_id
      push!(V_, new_id + 1)
      afssOut[f] = new_id + 1
      new_id = new_id + 2
    end

    # node ids of G from 0 to |V|
    G = VrpGraph(gvrp, V_, v_source, v_sink, (L, U))

    # resourves, R = R_M = {1,2} = {cap_res_id, fuel_res_id}}

    l_i_time, u_i_time = 0.0, T
    time_res_id = add_resource!(G, main=true)
    l_i_fuel, u_i_fuel = 0.0, β
    fuel_res_id = add_resource!(G, main=true)

    for i in V_
      if i in F
        set_resource_bounds!(G, afssIn[i], fuel_res_id, l_i_fuel, u_i_fuel)
        set_resource_bounds!(G, afssIn[i], time_res_id, l_i_time, u_i_time)
        set_resource_bounds!(G, afssOut[i], fuel_res_id, l_i_fuel, u_i_fuel)
        set_resource_bounds!(G, afssOut[i], time_res_id, l_i_time, u_i_time)
      else
        set_resource_bounds!(G, i, fuel_res_id, l_i_fuel, u_i_fuel)
        set_resource_bounds!(G, i, time_res_id, l_i_time, u_i_time)
      end
    end

    #edges
    for (i, j) in E
      if i in F && j in F && data.non_consec
        continue
      end
      # add arcs i - > j
      arc_id = nothing
      if i in F && j in F
        arc_id = add_arc!(G, afssOut[i], afssIn[j])
      elseif j in F
        arc_id = add_arc!(G, i, afssIn[j])
      elseif i in F
        arc_id = add_arc!(G, afssOut[i], j)
      else
        arc_id = add_arc!(G, i, j)
      end
      add_arc_var_mapping!(G, arc_id, x[(i, j)])
      set_arc_consumption!(G, arc_id, fuel_res_id, f(data,(i,j)))
      set_arc_consumption!(G, arc_id, time_res_id, t(data,(i,j)))
      # add arcs j - > i
      arc_id = nothing
      if i in F && j in F
        arc_id = add_arc!(G, afssOut[j], afssIn[i])
      elseif j in F
        arc_id = add_arc!(G, afssOut[j], i)
      elseif i in F
        arc_id = add_arc!(G, j, afssIn[i])
      else
        arc_id = add_arc!(G, j, i)
      end
      add_arc_var_mapping!(G, arc_id, x[(i, j)])
      set_arc_consumption!(G, arc_id, fuel_res_id, f(data,(i,j)))
      set_arc_consumption!(G, arc_id, time_res_id, t(data,(i,j)))
    end
    for f in F
      arc_id = add_arc!(G, afssIn[f], afssOut[f])
      set_arc_consumption!(G, arc_id, fuel_res_id, - β)
      set_arc_consumption!(G, arc_id, time_res_id, 0.0)
    end
    return G
  end

  G = build_graph()
  add_graph!(gvrp, G)

  set_vertex_packing_sets!(gvrp, [[(G, i)] for i in C])

  define_elementarity_sets_distance_matrix!(gvrp, G, [[d(data, ed(i, j)) for i in C] for j in C])

  set_branching_priority!(gvrp, "y", 1)
  set_branching_priority!(gvrp, "x", 1)

  #  function edge_ub_callback()
  #     for (i,j) in E
  #       e = (i,j)
  #        if i in C && j in C && get_value(gvrp.optimizer, x[e]) > 1.001
  #           println("Adding edge ub cut for e = ", e)
  #           add_dynamic_constr!(gvrp.optimizer, [x[e]], [1.0], <=, 1.0, "edge_ub")
  #        end
  #     end
  #  end
  #  add_cut_callback!(gvrp, edge_ub_callback, "edge_ub")

  function maxflow_mincut_callback()
    M = 100000
    # for all routes
    g = SparseMaxFlowMinCut.ArcFlow[]
    for (i, j) in E
      e = (i, j)
      value::Float64 = sum(get_value(gvrp.optimizer, x[e]))
      if value > 0.0001
        flow_::Int = trunc(floor(value, digits=5) * M)
        push!(g, SparseMaxFlowMinCut.ArcFlow(i, j, flow_)) # arc i -> j
        push!(g, SparseMaxFlowMinCut.ArcFlow(j, i, flow_)) # arc j -> i
      end
    end

    added_cuts = []
    s = data.depot_id
    for c in C
      maxFlow, flows, cut = SparseMaxFlowMinCut.find_maxflow_mincut(SparseMaxFlowMinCut.Graph(n, g), s, c)
      if (maxFlow / M) < (2 - 0.001) && !in(cut, added_cuts)
        set1, set2 = [], []
        [cut[i] == 1 ? push!(set1, i) : push!(set2, i) for i in 1:n]
        println(set1, "\\", set2)
        lhs_vars = [x[ed(i, j)] for i in set2 for j in set1 if ed(i, j) in E]
        lhs_coeff = [1.0 for i in set2 for j in set1 if ed(i, j) in E]
        add_dynamic_constr!(gvrp.optimizer, lhs_vars, lhs_coeff, >=, 2.0, "mincut")
        push!(added_cuts, cut)
        # constant LB routes
        setIn = c in set1 ? set1 : set2
        S₀ = vcat([data.depot_id], [i for i in setIn if i in C])
        println("Calculating routes LB for $(S₀)")
        nRoutesLB = calculateGVRP_NRoutesLB(data, S₀)
        println("LB = $nRoutesLB")
#        nRoutesLB = 1
        lhs_vars = [x[ed(i, j)] for i in S₀ for j in V if !in(j, S₀) && ed(i, j) in E]
        lhs_coeff = [1.0 for i in S₀ for j in V if !in(j, S₀) && ed(i, j) in E]
        add_dynamic_constr!(gvrp.optimizer, lhs_vars, lhs_coeff, >=, 2.0 * nRoutesLB, "mincut")
      end
    end
    if length(added_cuts) > 0 
      n_2_path_cuts_gh += length(added_cuts)
      println(">>>>> Add min cuts : ", length(added_cuts), " cut(s) added") 
    end
  end
  add_cut_callback!(gvrp, maxflow_mincut_callback, "mincut")

  function maxflow_mincut_time_callback()
    # solve model
    M = Model(solver = CplexSolver(
                                   CPX_PARAM_MIPDISPLAY=0,
                                   CPX_PARAM_SCRIND=0
                                  ))
    @variable(M, 0 <= y[1:n] <= 1, Int)
    @variable(M, 0 <= w[e in E] <= 1, Int)
    @variable(M, 0 <= z[e in E] <= 1, Int)
    x´ = Dict{Tuple{Int64,Int64},Float64}(e => get_value(gvrp.optimizer, x[e]) for e in E)
    @objective(M, Max, 2 * sum(z[e] * x´[e] * t(data, e) for e in E)/T - sum(w[e] * x´[e] for e in E))
    @constraint(M, update_z_1[e in E], z[e] >= y[e[1]])
    @constraint(M, update_z_2[e in E], z[e] >= y[e[2]])
    @constraint(M, update_z_3[(i, j) in E], y[i] + y[j] >= z[(i, j)])
    @constraint(M, update_w_1[(i, j) in E], w[(i, j)] >= y[i] - y[j])
    @constraint(M, update_w_2[(i, j) in E], w[(i, j)] >= y[j] - y[i])
    @constraint(M, update_w_3[(i, j) in E], 2 - (y[i] + y[j]) >= w[(i, j)])
    @constraint(M, update_w_4[(i, j) in E], y[i] + y[j] >= w[(i, j)])
    @constraint(M, y[1] == 0 )
    @constraint(M, sum(y[i] for i in C) >= 1)
    solve(M)
    # valid set
    if getobjectivevalue(M) > 0.001
      y´ = getvalue(y)
      w´ = getvalue(w)
      z´ = getvalue(z)
      # get bipartition
      setIn, setOut = [], []
      cutEdges = [e for e in E if w´[e] > 0.5]
      [y´[i] > 0.5 ? push!(setIn, i) : push!(setOut, i) for i in 1:n]
      println(setIn, "\\", setOut)
      # constant routes
      setIn = c in set1 ? set1 : set2
      S₀ = vcat([data.depot_id], [i for i in setIn if i in C])
      nRoutesLB = calculateGVRP_NRoutesLB(data, S₀)
      if sum(x´[e] for e in cutEdges) < nRoutesLB
        lhs_vars = [x[ed(i, j)] for i in S₀ for j in V if !in(j, S₀) && ed(i, j) in E]
        lhs_coeff = [1.0 for i in S₀ for j in V if !in(j, S₀) && ed(i, j) in E]
        add_dynamic_constr!(gvrp.optimizer, lhs_vars, lhs_coeff, >=, 2.0 * nRoutesLB, "constant_subcycle")

        lhs_vars = [x[e] for e in cutEdges]
        lhs_coeff = [1.0 for e in cutEdges]
        add_dynamic_constr!(gvrp.optimizer, lhs_vars, lhs_coeff, >=, 2.0 * nRoutesLB, "subcycle_time")
      end
      # variable routes
      lhs_vars = vcat([x[e] for e in E if w´[e] > 0.5], [x[e] for e in E if z´[e] > 0.5])
      lhs_coeff = vcat([1.0 for e in E if w´[e] > 0.5], [- t(data, e) * 2/T for e in E if z´[e] > 0.5])

      add_dynamic_constr!(gvrp.optimizer, lhs_vars, lhs_coeff, >=, 0.0, "mincuttime")
      n_k_path_cuts_time += 1
    end
  end
#  add_cut_callback!(gvrp, maxflow_mincut_time_callback, "mincuttime")
  
  function fuel_callback()
    # solve model
    M = Model(solver = CplexSolver(
                                   CPX_PARAM_MIPDISPLAY=0,
                                   CPX_PARAM_SCRIND=0
                                  ))
    @variable(M, 0 <= y[V] <= 1, Int)
    @variable(M, 0 <= w[e in E] <= 1, Int)
    @variable(M, 0 <= z[e in E] <= 1, Int)
    x´ = Dict{Tuple{Int64,Int64},Float64}(e => get_value(gvrp.optimizer, x[e]) for e in E)
    @objective(M, Max, sum(z[e] * x´[e] * f(data, e) for e in E) - (β/2) * sum(w[e] * x´[e] for e in E))
    @constraint(M, update_z_1[e in E], z[e] >= y[e[1]])
    @constraint(M, update_z_2[e in E], z[e] >= y[e[2]])
    @constraint(M, update_z_3[(i, j) in E], y[i] + y[j] >= z[(i, j)])
    @constraint(M, update_w_1[(i, j) in E], w[(i, j)] >= y[i] - y[j])
    @constraint(M, update_w_2[(i, j) in E], w[(i, j)] >= y[j] - y[i])
    @constraint(M, update_w_3[(i, j) in E], 2 - (y[i] + y[j]) >= w[(i, j)])
    @constraint(M, update_w_4[(i, j) in E], y[i] + y[j] >= w[(i, j)])
    @constraint(M, y[data.depot_id] == 0 )
    @constraint(M, sum(y[i] for i in V if !in(i, C)) == 0)
    @constraint(M, sum(y[i] for i in C) >= 1)
    solve(M)
    # valid set
    if getobjectivevalue(M) > 0.001
      y´ = getvalue(y)
      w´ = getvalue(w)
      z´ = getvalue(z)
      setIn = [y´[i] > 0.5 for i in C]
      cutEdges = [e for e in E if w´[e] > 0.5]
      println(setIn)
      # fuel
      lhs_vars = vcat([x[e] for e in E if z´[e] > 0.5], [x[e] for e in cutEdges])
      lhs_coeff = vcat([f(data, e) for e in E if z´[e] > 0.5], [- β/2 for e in cutEdges])
      add_dynamic_constr!(gvrp.optimizer, lhs_vars, lhs_coeff, <=, 0.0, "fuel")
      # routes
      nRoutesLB = calculateGVRP_NRoutesLB(data, vcat([data.depot_id], [i for i in C if y´[i] > 0.5]))
      if sum(x´[e] for e in cutEdges) < nRoutesLB
        lhs_vars = [x[e] for e in cutEdges]
        lhs_coeff = [1.0 for e in cutEdges]
        add_dynamic_constr!(gvrp.optimizer, lhs_vars, lhs_coeff, >=, 2.0 * nRoutesLB, "subcycle_fuel")
      end
      n_energy_cuts += 1
    end
  end
#  add_cut_callback!(gvrp, fuel_callback, "fuel")

  return (gvrp, x, y)
end
