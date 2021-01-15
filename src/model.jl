using CPLEX

## debug by marcos

function cplex(data::DataGVRP)
  E = edges(data) # set of edges of the input graph G′
  n = nb_vertices(data)
  V = [i for i in 1:n] # set of vertices of the input graph G′
  V′ = [i for i in 0:n] # V ⋃ {0}, where 0 is a dummy vertex

  β = data.β
  C = data.C # Set of customers vertices
  F = data.F # Set of AFSs vertices
  T = data.T # Set of AFSs vertices

  ed(i, j) = i < j ? (i, j) : (j, i)

  gvrp_ = Model(solver = CplexSolver())
  F´ = deepcopy(F) 
  popfirst!(F´)
  @variable(gvrp_, x[e in E] >= 0, Int)
  @variable(gvrp_, y[i in F] >= 0, Int)
  @objective(gvrp_, Min, sum(d(data, e) * x[e] for e in E))
  @constraint(gvrp_, deg[i in C], sum(x[e] for e in δ(data, i)) == 2.0)
  @constraint(gvrp_, hotel_deg[i in F], sum(x[e] for e in δ(data, i)) == 2*y[i])
  @constraint(gvrp_, no_edge_between_afss[f in F´], sum(x[(f, r)] for r in F´ if (f, r) in data.G′.E) == 0.0)
  @constraint(gvrp_, y[F[1]] >= 1)
  time_user_cut_added_cuts = []
  time_lazy_added_cuts = []
  fuel_user_cut_added_cuts = []
  fuel_lazy_added_cuts = []
  function separa(cb, corte)
    x´ = Dict{Tuple{Int64,Int64},Float64}(e => getvalue(x[e]) for e in E)
    # fuel separation
    # dfs
    """
    edges = []
    function dfs(i::Int64)
      for e in δ(data, i) 
        if x´[e] > 0.0001 
          if e[1] in F && e[1] != i
            # get edges
               
          end
          push!(edge, e)
        end
      end
    end
    for f in F
    end
    """
    # time separation
    M = Model(solver = CplexSolver(
                                   CPX_PARAM_MIPDISPLAY=0,
                                   CPX_PARAM_SCRIND=0
                                  ))
    @variable(M, 0 <= y[1:n] <= 1, Int)
    @variable(M, 0 <= w[e in E] <= 1, Int)
    @variable(M, 0 <= z[e in E] <= 1, Int)
    @objective(M, Max, (2/T) * sum(z[e] * x´[e] * t(data, e) for e in E) - sum(w[e] * x´[e] for e in E))
    @constraint(M, update_z_1[e in E], z[e] >= y[e[1]])
    @constraint(M, update_z_2[e in E], z[e] >= y[e[2]])
    @constraint(M, update_z_3[(i, j) in E], y[i] + y[j] >= z[(i, j)])
    @constraint(M, update_w_1[(i, j) in E], w[(i, j)] >= y[i] - y[j])
    @constraint(M, update_w_2[(i, j) in E], w[(i, j)] >= y[j] - y[i])
    @constraint(M, update_w_3[(i, j) in E], 2 - (y[i] + y[j]) >= w[(i, j)])
    @constraint(M, update_w_4[(i, j) in E], y[i] + y[j] >= w[(i, j)])
    @constraint(M, y[F[1]] == 0 )
    @constraint(M, sum(y[i] for i in C) >= 1)
    solve(M)
    # valid set
    if getobjectivevalue(M) > 0.001
      y´ = getvalue(y)
      w´ = getvalue(w)
      z´ = getvalue(z)
      # get bipartition
      setIn, setOut = [], []
      [y´[i] > 0.5 ? push!(setIn, i) : push!(setOut, i) for i in 1:n]

      #       println(y´)
      #       println(w´)
      #       println(z´)
      #       println("S: ", setIn)
      #       println("V\\S: ", setOut)
      #       println([x[e] for e in E if w´[e] > 0.5])
      #      println()
      #      println([x[e] for e in E if z´[e] > 0.5])
      lhs_vars = vcat([x[e] for e in E if w´[e] > 0.5], [x[e] for e in E if z´[e] > 0.5])
      lhs_coeff = vcat([1.0 for e in E if w´[e] > 0.5], [- t(data, e) * 2/T for e in E if z´[e] > 0.5])

      if corte
        if setIn in time_user_cut_added_cuts
          println("=======>REPEATED USER CUT!!!")
        end
        push!(time_user_cut_added_cuts, setIn)
        @usercut(cb, sum(lhs_vars[i] * lhs_coeff[i] for i in 1:length(lhs_vars)) >= 0.0)
        unsafe_store!(cb.userinteraction_p, convert(Cint,2), 1)
      else
        if setIn in time_lazy_added_cuts
          println("=======>REPEATED LAZY CUT!!!")
        end
        push!(time_lazy_added_cuts, setIn)
        @lazyconstraint(cb, sum(lhs_vars[i] * lhs_coeff[i] for i in 1:length(lhs_vars)) >= 0.0)
      end
    end
  end
  function separa_corte(cb)
    separa(cb, true)
  end
  addcutcallback(gvrp_, separa_corte)
  function separa_restr(cb)
    separa(cb, false)
  end
  addlazycallback(gvrp_, separa_restr)
  solve(gvrp_)
  obj = getobjectivevalue(gvrp_)
end

function build_model(data::DataGVRP)
  E = edges(data) # set of edges of the input graph G′
  n = nb_vertices(data)
  V = [i for i in 1:n] # set of vertices of the input graph G′
  V′ = [i for i in 0:n] # V ⋃ {0}, where 0 is a dummy vertex

  β = data.β
  C = data.C # Set of customers vertices
  F = data.F # Set of AFSs vertices
  T = data.T # Set of AFSs vertices
  F´ = deepcopy(F) 
  popfirst!(F´)

  ed(i, j) = i < j ? (i, j) : (j, i)

# println("Customers ")
# for i in data.C
#   println(i, " ", data.G′.V′[i], ", δ($i) = $(δ(data, i))")
# end
# println("AFSs ")
# for i in data.F
#   println(i, " ", data.G′.V′[i], ", δ($i) = $(δ(data, i))")
# end
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
  @variable(gvrp.formulation, x[e in E], Int)
  @variable(gvrp.formulation, 2 * length(C) >= y[i in F] >= 0, Int)
  @objective(gvrp.formulation, Min, sum(d(data, e) * x[e] for e in E))
  @constraint(gvrp.formulation, deg[i in C], sum(x[e] for e in δ(data, i)) == 2.0)
  @constraint(gvrp.formulation, hotel_deg[i in F], sum(x[e] for e in δ(data, i)) == 2*y[i])
  @constraint(gvrp.formulation, no_edge_between_afss[f in F´], sum(x[(f, r)] for r in F´ if (f, r) in E) == 0.0)
  @constraint(gvrp.formulation, y[F[1]] >= 1)

  # Build the model directed graph G=(V,A)
  function build_graph()

    v_source = v_sink = 0
#    L = U = length(F) # max and min number of paths is equal to number of AFSs
    L = length(F´)
    U = length(C) # max and min number of paths is equal to number of AFSs

    # node ids of G from 0 to |V|
    G = VrpGraph(gvrp, V′, v_source, v_sink, (L, U))
    # resourves, R = R_M = {1,2} = {cap_res_id, fuel_res_id}}
    time_res_id = add_resource!(G, main=true)
    fuel_res_id = add_resource!(G, main=true)
    for i in V′
      l_i_fuel, u_i_fuel = 0.0, β
      set_resource_bounds!(G, i, fuel_res_id, l_i_fuel, u_i_fuel)
      l_i_time, u_i_time = 0.0, T
      set_resource_bounds!(G, i, time_res_id, l_i_time, u_i_time)
    end

    # Build set of arcs A from E′ (two arcs for each edge (i,j))
    for f in F # setting the arcs between source, sink, and black vertices
      # source -> i(AFS)
      arc_id = add_arc!(G, v_source, f)
      set_arc_consumption!(G, arc_id, time_res_id, data.G′.V′[f].service_time/2.0)
      set_arc_consumption!(G, arc_id, fuel_res_id, 0.0)
      # i(AFS) -> sink
      arc_id = add_arc!(G, f, v_sink)
      set_arc_consumption!(G, arc_id, time_res_id, data.G′.V′[f].service_time/2.0)
      set_arc_consumption!(G, arc_id, fuel_res_id, 0.0)
    end
    for (i, j) in E
      # add arcs i - > j
      arc_id = add_arc!(G, i, j)
      add_arc_var_mapping!(G, arc_id, x[(i, j)])

      if i in F && j in F && i != F[1] && j != F[1]
        set_arc_consumption!(G, arc_id, fuel_res_id, β + 1)
      elseif i in F && j in F && i == F[1] || j == F[1]
        set_arc_consumption!(G, arc_id, fuel_res_id, - f(data,(i,j)))
      else
        set_arc_consumption!(G, arc_id, fuel_res_id,  f(data,(i,j)))
      end
      set_arc_consumption!(G, arc_id, time_res_id, t(data,(i,j)))
      # add arcs j - > i
      arc_id = add_arc!(G, j, i)
      add_arc_var_mapping!(G, arc_id, x[(i, j)])

      if i in F && j in F && i != F[1] && j != F[1]
        set_arc_consumption!(G, arc_id, fuel_res_id, β + 1)
      elseif i in F && j in F && i == F[1] || j == F[1]
        set_arc_consumption!(G, arc_id, fuel_res_id, - f(data,(i,j)))
      else
        set_arc_consumption!(G, arc_id, fuel_res_id,  f(data,(i,j)))
      end
      set_arc_consumption!(G, arc_id, time_res_id, t(data,(i,j)))
    end
    return G
  end

  G = build_graph()
  add_graph!(gvrp, G)

  set_vertex_packing_sets!(gvrp, [[(G, i)] for i in C])

  set_additional_vertex_elementarity_sets!(gvrp, [(G,[f]) for f in data.F])
  define_elementarity_sets_distance_matrix!(gvrp, G, [[ed(i, j) in data.G′.E ? d(data, ed(i, j)) : 0.0 for i in V] for j in V])
  add_capacity_cut_separator!(gvrp, [ ([(G, i)], 2.0*data.G′.V′[i].service_time) for i in C], 2.0*floor(data.T) )

  set_branching_priority!(gvrp, "x", 1)
  #set_branching_priority!(gvrp, "y", 1)

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
    s = F[1]
    for c in C
      maxFlow, flows, cut = SparseMaxFlowMinCut.find_maxflow_mincut(SparseMaxFlowMinCut.Graph(n, g), s, c)
      if (maxFlow / M) < (2 - 0.001) && !in(cut, added_cuts)
        set1, set2 = [], []
        [cut[i] == 1 ? push!(set1, i) : push!(set2, i) for i in 1:n]
        println(set1, "\\", set2)
        lhs_vars = [x[ed(i, j)] for i in set2 for j in set1 if ed(i, j) in E]
        lhs_coeff = [1.0 for i in set2 for j in set1 if ed(i, j) in E]
        setIn = c in set1 ? set1 : set2
        add_dynamic_constr!(gvrp.optimizer, lhs_vars, lhs_coeff, >=, 2.0, "mincut")

        push!(added_cuts, cut)
      end
    end
    if length(added_cuts) > 0 
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
    @objective(M, Max, (2/T) * sum(z[e] * x´[e] * t(data, e) for e in E) - sum(w[e] * x´[e] for e in E))
    @constraint(M, update_z_1[e in E], z[e] >= y[e[1]])
    @constraint(M, update_z_2[e in E], z[e] >= y[e[2]])
    @constraint(M, update_z_3[(i, j) in E], y[i] + y[j] >= z[(i, j)])
    @constraint(M, update_w_1[(i, j) in E], w[(i, j)] >= y[i] - y[j])
    @constraint(M, update_w_2[(i, j) in E], w[(i, j)] >= y[j] - y[i])
    @constraint(M, update_w_3[(i, j) in E], 2 - (y[i] + y[j]) >= w[(i, j)])
    @constraint(M, update_w_4[(i, j) in E], y[i] + y[j] >= w[(i, j)])
    @constraint(M, y[F[1]] == 0 )
    @constraint(M, sum(y[i] for i in C) >= 1)
    solve(M)
    # valid set
    if getobjectivevalue(M) > 0.001
      y´ = getvalue(y)
      w´ = getvalue(w)
      z´ = getvalue(z)
      # get bipartition
      setIn, setOut = [], []
      [y´[i] > 0.5 ? push!(setIn, i) : push!(setOut, i) for i in 1:n]

#     println(y´)
#     println(w´)
#     println(z´)
     println("S: ", setIn)
     println("V\\S: ", setOut)
#     println([x[e] for e in E if w´[e] > 0.5])
#     println()
#     println([x[e] for e in E if z´[e] > 0.5])

      lhs_vars = vcat([x[e] for e in E if w´[e] > 0.5], [x[e] for e in E if z´[e] > 0.5])
      lhs_coeff = vcat([1.0 for e in E if w´[e] > 0.5], [- t(data, e) * 2/T for e in E if z´[e] > 0.5])

      add_dynamic_constr!(gvrp.optimizer, lhs_vars, lhs_coeff, >=, 0.0, "mincuttime")
      println(">>>>> Add min cuts time: ", 1, " cut(s) added") 
    end
    """
    M = 100000
    # for all routes
    g = SparseMaxFlowMinCut.ArcFlow[]
    for (i, j) in E
      e = (i, j)
      value::Float64 = get_value(gvrp.optimizer, x[e])
      if value > 0.0001
        flow_::Int = trunc(floor(value, digits=5) * M)
        push!(g, SparseMaxFlowMinCut.ArcFlow(i, j, flow_)) # arc i -> j
        push!(g, SparseMaxFlowMinCut.ArcFlow(j, i, flow_)) # arc j -> i
      end
    end

    added_cuts = []
    s = F[1]
    for c in C
      maxFlow, flows, cut = SparseMaxFlowMinCut.find_maxflow_mincut(SparseMaxFlowMinCut.Graph(n, g), s, c)
      set1, set2 = [], []
      [cut[i] == 1 ? push!(set1, i) : push!(set2, i) for i in 1:n]
      setIn = c in set1 ? set1 : set2
      setOut = c in set1 ? set2 : set1
      if (maxFlow/M) + 0.001 < 2 * sum(get_value(gvrp.optimizer, x[(i, j)]) * t(data, (i, j)) for (i, j) in E if i in setIn || j in setIn)/T  && !in(cut, added_cuts)
        println(setIn, " maxflow: ", (maxFlow/M), " cut:", 2 * sum(get_value(gvrp.optimizer, x[(i, j)]) * t(data, (i, j)) for (i, j) in E if i in setIn || j in setIn)/T)
        lhs_vars = vcat([x[(i, j)] for i in set2 for j in set1 if (i, j) in E], [x[(i, j)] for (i, j) in E if i in setIn || j in setIn])
        lhs_coeff = vcat([1.0 for i in set2 for j in set1 if (i, j) in E], [- 2 * t(data, (i, j))/T for (i, j) in E if i in setIn || j in setIn])

        add_dynamic_constr!(gvrp.optimizer, lhs_vars, lhs_coeff, >=, 0.0, "mincuttime")

        push!(added_cuts, cut)
      end
    end
    if length(added_cuts) > 0 
      println(">>>>> Add min cuts : ", length(added_cuts), " cut(s) added") 
    end
    """
  end
  #add_cut_callback!(gvrp, maxflow_mincut_callback, "mincuttime")
  return (gvrp, x)
end
