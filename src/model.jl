function build_model(data::DataGVRP)
    
    E = edges(data) # set of edges of the input graph G′
    n = nb_vertices(data)
    V = [i for i in 1:n] # set of vertices of the input graph G′
    V′ = [i for i in 0:n] # V ⋃ {0}, where 0 is a dummy vertex

    β = data.β
    C = data.C # Set of customers vertices
    F = data.F # Set of AFSs vertices
    F´ = deepcopy(F) 
    popfirst!(F´)

    T = data.T # General time limit
    M = data.M # Set of vehicles
    
    K = M
    
    ed(i, j) = i < j ? (i, j) : (j, i)

    # data.E′
    # Formulation
    gvrp = VrpModel()
    @variables(gvrp.formulation, begin
                 0 <= x[i in V , j in V, k in M] <= 1, Int
                 0 <= e[i in C] <= data.β
               end)

    @objective(gvrp.formulation, Min, sum( data.G′.cost[ed(i, j)] * x[i,j,k] for i in V, j in V, k in M if i != j && !((i, j) in data.E′)  ) )

    @constraints(gvrp.formulation, begin
                  deg_6_2[k in M], sum(x[1,j,k] for j in V if (j!=1 && !((1, j) in data.E′) ) ) <= 1.0

                  deg_6_3[k in M, i in V], sum(x[i,j,k] for j in V if (j!=i && !((i, j) in data.E′) ) ) == sum(x[j,i,k] for j in V if (j!=i && !((j, i) in data.E′) ) )

                  deg_6_4[i in C], sum( x[i,j,k] for j in V, k in M if (j!=i && !((i, j) in data.E′) ) ) == 1.0

                  deg_6_6[ [(i, j) for i in C, j in C if (i, j) in data.E′] , k in M], e[j] <= e[i] - f(data, ed(i, j))*x[i,j,k] + data.β*(1.0 - x[i,j,k]) + f(data, ed(j, i))*x[j,i,k]

                  deg_6_7_1[j in C],  e[j] <= data.β - sum( f(data, ed( ff, j))*x[ff,j,k] for k in M , ff in F´ if !( (ff, j) in data.E′ ) )
                  
                  deg_6_7_2[j in C],  sum( f(data, ed(j, ff))*x[j,ff,k] for k in M , ff in F´ if !( (j, ff) in data.E′ ) ) <= e[j]

                  deg_6_9[k in M], sum(x[i, j, k] * (t(data, ed(i, j)) + data.G′.V′[i].service_time) for i in V, j in V if (i!=j && !((i, j) in data.E′) ) ) <= T
                  
                  #prepro[e in L, k in M], x[e[1], e[2], k] == 0
                end)
    
    #println(gvrp.formulation)

    # Build the model directed graph G=(V,A)
    function build_graph(k::Int64)
        v_source = v_sink = 0
        L = 0 
        U = length(C) # max and min number of paths is equal to number of AFSs

        # node ids of G from 0 to |V|
        G = VrpGraph(gvrp, V′, v_source, v_sink, (L, U))
        # resourves, R = R_M = {1,2} = {cap_res_id, fuel_res_id}}
        time_res_id = add_resource!(G, main=true)
        fuel_res_id = add_resource!(G, main=true)

        for i in V′
            # l_i, u_i = 0.0, Float64(Q) # accumulated resource consumption interval [l_i, u_i] for the vertex i
            # set_resource_bounds!(G, i, cap_res_id, l_i, u_i)

            l_i_fuel, u_i_fuel = 0.0, β 
            set_resource_bounds!(G, i, fuel_res_id, l_i_fuel, u_i_fuel)
            l_i_time, u_i_time = 0.0, T
            set_resource_bounds!(G, i, time_res_id, l_i_time, u_i_time)
        end

        # Build set of arcs A from E′ (two arcs for each edge (i,j))
        for f in F´ # setting the arcs between source, sink, and black vertices
            # source -> i(AFS)
            arc_id = add_arc!(G, v_source, f)
            set_arc_consumption!(G, arc_id, time_res_id, data.G′.V′[f].service_time/2)
            set_arc_consumption!(G, arc_id, fuel_res_id, 0.0)
            # i(AFS) -> sink
            arc_id = add_arc!(G, f, v_sink)
            set_arc_consumption!(G, arc_id, time_res_id, data.G′.V′[f].service_time/2)
            set_arc_consumption!(G, arc_id, fuel_res_id, 0.0)
        end

        for i in V
            for j in V
                if !((i, j) in data.E′)
                    # resource comsuption for the R = 1 
                    # q_1 = 0.5
                    # if (i in W) && (j in W)
                    #     q_1 = 1.0
                    # elseif (i in B) && (j in B)
                    #     q_1 = Q
                    # end

                    # add arcs i - > j
                    arc_id = add_arc!(G, i, j)
                    add_arc_var_mapping!(G, arc_id, x[i, j, k])
                     
                    if i in F && j in F && i != F[1] && j != F[1]
                        set_arc_consumption!(G, arc_id, fuel_res_id, β + 1)
                    elseif i in F && j in F && i == F[1] || j == F[1]
                        set_arc_consumption!(G, arc_id, fuel_res_id,  - f(data,ed(i,j)))
                    else
                        set_arc_consumption!(G, arc_id, fuel_res_id,  f(data,ed(i,j)))
                    end

                    set_arc_consumption!(G, arc_id, time_res_id, ((data.G′.V′[i].service_time + data.G′.V′[j].service_time)/2) + t(data,ed(i,j)))
                    # add arcs j - > i
                    
                    arc_id = add_arc!(G, j, i)
                    add_arc_var_mapping!(G, arc_id, x[i, j, k])

                    if i in F && j in F && i != F[1] && j != F[1]
                        set_arc_consumption!(G, arc_id, fuel_res_id, β + 1)
                    elseif i in F && j in F && i == F[1] || j == F[1]
                        set_arc_consumption!(G, arc_id, fuel_res_id,  - f(data,ed(i,j)))
                    else
                        set_arc_consumption!(G, arc_id, fuel_res_id,  f(data,ed(i,j)))
                    end

                    set_arc_consumption!(G, arc_id, time_res_id, ((data.G′.V′[i].service_time + data.G′.V′[j].service_time)/2) + t(data,ed(i,j)))
                end
            end
        end
        return G
    end

    Graphs = []
    for k in K
        push!(Graphs, build_graph(k))
    end

    for G in Graphs
      add_graph!(gvrp, G)
    end

    set_vertex_packing_sets!(gvrp, [[(Graphs[k], i) for k in K] for i in C])

    [define_elementarity_sets_distance_matrix!(gvrp, Graphs[k], [[d(data,ed(i, j)) for i in C if !((i, j) in data.E′) ] for j in C ] ) for k in K]
    # add_capacity_cut_separator!(gvrp, [ ([(G, i)], 1.0) for i in W], Float64(Q))

    set_branching_priority!(gvrp, "x", 1)

    function maxflow_mincut_callback()
        println("OKAY 7")
        M = 100000
        added_cuts = []

        # for all routes
        """g = SparseMaxFlowMinCut.ArcFlow[]
        for i in V
            for j in V 
                value::Float64 = sum(get_value(gvrp.optimizer, x[i, j, k]) for k in K)
                #            value::Float64 = get_value(gvrp.optimizer, x[e] * t(data, ee))
                if value > 0.0001
                    flow_::Int = trunc(floor(value, digits=5) * M)
                    push!(g, SparseMaxFlowMinCut.ArcFlow(i, j, flow_)) # arc i -> j
                    push!(g, SparseMaxFlowMinCut.ArcFlow(j, i, flow_)) # arc j -> i
                end
            end
        end

        s = F[1]
        for c in C
            maxFlow, flows, cut = SparseMaxFlowMinCut.find_maxflow_mincut(SparseMaxFlowMinCut.Graph(n, g), s, c)
#            if (maxFlow / M) > (T - 0.001) && !in(cut, added_cuts)
            if (maxFlow / M) < (2 - 0.001) && !in(cut, added_cuts)
                set1, set2 = [], []
                [cut[i] == 1 ? push!(set1, i) : push!(set2, i) for i in 1:n]

                lhs_vars = [x[i, j, k] for i in set2 for j in set1 for k in K]
                lhs_coeff = [1.0 for i in set2 for j in set1 for k in K]

#                add_dynamic_constr!(gvrp.optimizer, lhs_vars, lhs_coeff, >=, 2.0 * floor(ceil(sum(data.G′.V′[i].service_time for i in (c in set1 ? set1 : set2) if i in C)/T)), "mincut")
                add_dynamic_constr!(gvrp.optimizer, lhs_vars, lhs_coeff, >=, 2.0, "mincut")

                push!(added_cuts, cut)
            end
        end
        """
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
            #println(q)
            #println("--------------------------------------------------------------------------------------------------------------------------------------------------")
            while length(q) > 0
              i = pop!(q)
              for j in V 
                #if i != j && get_value(gvrp.optimizer, x[i, j, k]) > 0.0001 || 
                #    get_value(gvrp.optimizer, x[j, i, k]) > 0.0001 &&     !visited[j]
                if i != j && !((i, j) in data.E′) && get_value(gvrp.optimizer, x[i, j, k]) > 0.0001 && !visited[j]
                  visited[j] = true
                  push!(q, j)
                  push!(comp, j)
                end
              end
            end
            #println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
            
            if !(s in comp) && length(comp) > 1
              println("Route $k S: ", comp)
              lhs_vars = [x[i, j, k] for i in comp for j in V if !(j in comp)]
              lhs_coeff = [1.0 for i in comp for j in V if !(j in comp)]
              for i in comp
                if i in C
                  for j in comp
                    if i != j && !((i, j) in data.E′)
                      add_dynamic_constr!(gvrp.optimizer, 
                                          vcat(lhs_vars, [x[i, j, k]]), 
                                          vcat(lhs_coeff, [-1.0]), 
                                          >=, 0.0, "mincut")
                      println("add cut")
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
