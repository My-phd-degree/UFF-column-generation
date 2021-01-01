function build_model(data::DataGVRP)

    E = edges(data) # set of edges of the input graph G′
    n = nb_vertices(data)
    V = [i for i in 1:n] # set of vertices of the input graph G′
    V′ = [i for i in 0:n] # V ⋃ {0}, where 0 is a dummy vertex

    β = data.β
    C = data.C # Set of customers vertices
    F = data.F # Set of AFSs vertices
    T = data.T # Set of AFSs vertices
    K = 1:length(C)

    ed(i, j) = i < j ? (i, j) : (j, i)
    # Formulation
    gvrp = VrpModel()
    @variables(gvrp.formulation, begin
                 x[k in K, e in E], Int
                 y[k in K, i in F] >= 0, Int
               end)
    @objective(gvrp.formulation, Min, sum(d(data, e) * x[k, e] for e in E for k in K))
    @constraints(gvrp.formulation, begin
                  deg[i in C], sum(x[k, e] for e in δ(data, i) for k in K) == 2.0
                  afs_deg[i in F, k in K], sum(x[k, e] for e in δ(data, i)) == 2*y[k, i]
                  route_time_limit[k in K], sum(x[k, e] * t(data, e) for e in E) <= T
                  route_used_at_most_once[k in K], sum(x[k, (1, i)] for i in V if i != 1) <= 2
                  no_edge_between_afs[f in F], sum(x[k, ed(f, r)] for r in F for k in K if r != 1 && f != 1 && f != r) == 0
                  sum(y[k, F[1]] for k in K) >= 1
                end)
    #println(gvrp.formulation)

    # Build the model directed graph G=(V,A)
    function build_graph(k::Int)

        v_source = v_sink = 0
        L = U = length(F) # max and min number of paths is equal to number of AFSs

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
            set_arc_consumption!(G, arc_id, time_res_id, data.G′.V′[f].service_time/2)
            set_arc_consumption!(G, arc_id, fuel_res_id, 0.0)
            # i(AFS) -> sink
            arc_id = add_arc!(G, f, v_sink)
            set_arc_consumption!(G, arc_id, time_res_id, data.G′.V′[f].service_time/2)
            set_arc_consumption!(G, arc_id, fuel_res_id, 0.0)
        end
        for (i, j) in E
            if i in F && j in F && i != 1 && j != 1
              continue
            end

            # add arcs i - > j
            arc_id = add_arc!(G, i, j)
            add_arc_var_mapping!(G, arc_id, x[k, (i, j)])
             
            set_arc_consumption!(G, arc_id, fuel_res_id, f(data,(i,j)))
            set_arc_consumption!(G, arc_id, time_res_id, t(data,(i,j)))
            # add arcs j - > i
            arc_id = add_arc!(G, j, i)
            add_arc_var_mapping!(G, arc_id, x[k, (i, j)])

            set_arc_consumption!(G, arc_id, fuel_res_id, f(data,(i,j)))
            set_arc_consumption!(G, arc_id, time_res_id, t(data,(i,j)))
        end
        return G
    end

    Graphs = [build_graph(k) for k in K] 
    #    G = build_graph()
    for G in Graphs
      add_graph!(gvrp, G)
      #println(G)

#      set_vertex_packing_sets!(gvrp, [[(G, i)] for i in C])
      #set_additional_vertex_elementarity_sets!(gvrp, [(G,[i]) for i in B])

    end
    set_vertex_packing_sets!(gvrp, [[(Graphs[k], i) for k in K] for i in C])
    [define_elementarity_sets_distance_matrix!(gvrp, Graphs[k], [[d(data,ed(i, j)) for i in C] for j in C]) for k in K]

"""
#    add_capacity_cut_separator!(gvrp, [ ([(G, i)], 1.0) for i in C], Float64(Q))
    @expression(gvrp.formulation, edge[e in E], sum(x[k, e] for k in K))
    set_branching_priority!(gvrp, edge, "edge", 1)

    @expression(gvrp.formulation, vehNum[k in K], sum(0.5 * x[k, (1,i)] for i in 2:n))
    set_branching_priority!(gvrp, vehNum, "vehNum", 3)

    @expression(gvrp.formulation, assign[k in K, i in 2:n], sum(0.5 * x[k, e] for e in δ(data, i)))
    set_branching_priority!(gvrp, assign, "assign", 1)
"""
    set_branching_priority!(gvrp, "x", 1)
    set_branching_priority!(gvrp, "y", 1)

    function maxflow_mincut_callback()
        M = 100000
        # for all routes
        g = SparseMaxFlowMinCut.ArcFlow[]
        for (i, j) in E
            e = (i, j)
            value::Float64 = sum(get_value(gvrp.optimizer, x[k, e]) for k in K)
            #            value::Float64 = get_value(gvrp.optimizer, x[e] * t(data, ee))
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
#            if (maxFlow / M) > (T - 0.001) && !in(cut, added_cuts)
            if (maxFlow / M) < (2 - 0.001) && !in(cut, added_cuts)
                set1, set2 = [], []
                [cut[i] == 1 ? push!(set1, i) : push!(set2, i) for i in 1:n]

                lhs_vars = [x[k, ed(i, j)] for i in set2 for j in set1 for k in K]
                lhs_coeff = [1.0 for i in set2 for j in set1 for k in K]

#                add_dynamic_constr!(gvrp.optimizer, lhs_vars, lhs_coeff, >=, 2.0 * floor(ceil(sum(data.G′.V′[i].service_time for i in (c in set1 ? set1 : set2) if i in C)/T)), "mincut")
                add_dynamic_constr!(gvrp.optimizer, lhs_vars, lhs_coeff, >=, 2.0, "mincut")

                push!(added_cuts, cut)
            end
        end
        # for each route
        for k in K
          s = F[1]
          visited = [false for i in V]
          for c in C
            if visited[c]
              continue
            end
            visited[c] = true
            comp = [c]
            q = [c]
            while length(q) > 0
              i = pop!(q)
              for j in V 
                if j != i && get_value(gvrp.optimizer, x[k, ed(i, j)]) > 0.0001 && !visited[j]
                  visited[j] = true
                  push!(q, j)
                  push!(comp, j)
                end
              end
            end
            if !(s in comp) && length(comp) > 1
              println("Route $k S: ", comp)
              lhs_vars = [x[k, ed(i, j)] for i in comp for j in V if !(j in comp)]
              lhs_coeff = [1.0 for i in comp for j in V if !(j in comp)]
              for i in comp
                if i in C
                  for j in comp
                    if i != j
                      add_dynamic_constr!(gvrp.optimizer, 
                                          vcat(lhs_vars, [x[k, ed(i, j)]]), 
                                          vcat(lhs_coeff, [-1.0]), 
                                          >=, 0.0, "mincut")

                    end
                  #push!(added_cuts, cut)
                  end
                end
              end
            end
          end
        end
        """
        if length(added_cuts) > 0 
          println(">>>>> Add min cuts : ", length(added_cuts), " cut(s) added") 
        end
        """
    end
    add_cut_callback!(gvrp, maxflow_mincut_callback, "mincut")

    return (gvrp, x)
end
