function build_model(data::DataGVRP)
    
    # os primeiros F vertices são vertices de abastecimento 

    E = edges(data) # set of edges of the input graph G´
    n = nb_vertices(data)
    V = [i for i in 1:n] # set of vertices of the input graph G´
    V´ = [i for i in 0:n] # V ⋃ {0}, where 0 is a dummy vertex

    β = data.β
    C = data.C # Set of customers vertices
    F = data.F # Set of AFSs vertices
    F´ = deepcopy(F) 
    popfirst!(F´)

    T = data.T # General time limit
    M = data.M # Set of vehicles
    K = M
    FM = data.FM

    if length(data.E´) > 0
        println("COM ELIMINACAO DE ARESTA")
        println(data.E´)
    else
        println("SEM ELIMINACAO DE ARESTA")
        println(data.E´)
    end

    println("-------")
    println("AFs:")
    for f in F´ 
        for k in K
            #FM[f][k] = 1
        end
        println(f)
    end
    println("-------")

    ed(i, j) = i < j ? (i, j) : (j, i)

    # Formulation
    gvrp = VrpModel()
    @variables(gvrp.formulation, begin
                 0 <= x[i in V, j in V, k in M] <= 1, Int
                 0 <= e[i in C] <= data.β
               end)

    @objective(gvrp.formulation, Min, sum( data.G´.cost[ed(i, j)] * x[i,j,k] for i in V, j in V, k in M if i != j && !((i, j) in data.E´)  ) )

    @constraints(gvrp.formulation, begin
                  deg_6_2[k in M], sum(x[1,j,k] for j in V if (j!=1 && !((1, j) in data.E´) ) ) <= 1.0

                  deg_6_3[k in M, i in V], sum(x[i,j,k] for j in V if (j!=i && !((i, j) in data.E´) ) ) == sum(x[j,i,k] for j in V if (j!=i && !((j, i) in data.E´) ) )

                  deg_6_4[i in C], sum( x[i,j,k] for j in V, k in M if (j!=i && !((i, j) in data.E´) ) ) == 1.0

                  deg_6_6[ i in C, j in C, k in M], if !( (i, j) in data.E´) e[j] else 0.0 end <=  if !( (i, j) in data.E´) e[i] - f(data, ed(i, j))*x[i,j,k] + data.β*(1.0 - x[i,j,k]) + f(data, ed(j, i))*x[j,i,k] else 0.0 end

                  deg_6_7_1[j in C],  e[j] <= data.β - sum( f(data, ed( ff, j))*x[ff,j,k] for k in M , ff in data.F if !( (ff, j) in data.E´ ) )
                  
                  deg_6_7_2[j in C],  sum( f(data, ed(j, ff))*x[j,ff,k] for k in M , ff in data.F if !( (j, ff) in data.E´ ) ) <= e[j]

                  deg_6_9[k in M], sum(x[i, j, k] * (t(data, ed(i, j)) + data.G´.V´[i].service_time) for i in V, j in V if (i!=j && !((i, j) in data.E´) ) ) <= T

                  #prepro[e in L, k in M], x[e[1], e[2], k] == 0
                end)
    
    #println(gvrp.formulation)

    # Build the model directed graph G=(V,A)
    function build_graph()
        v_source = v_sink = 0
        #L = length(F´) 
        #U = length(V) # max and min number of paths is equal to number of AFSs
        #L = length(C)

        #L = 0 
        #U = length(C)

        L = U = length(F´)
        
        # node ids of G from 0 to |V|
        G = VrpGraph(gvrp, V´, v_source, v_sink, (L, U))
        # resourves, R = R_M = {1,2} = {cap_res_id, fuel_res_id}}
        time_res_id = add_resource!(G, main=true)
        fuel_res_id = add_resource!(G, main=true)

        for i in V´
            # l_i, u_i = 0.0, Float64(Q) # accumulated resource consumption interval [l_i, u_i] for the vertex i
            # set_resource_bounds!(G, i, cap_res_id, l_i, u_i)

            l_i_fuel, u_i_fuel = 0.0, β 
            set_resource_bounds!(G, i, fuel_res_id, l_i_fuel, u_i_fuel)
            l_i_time, u_i_time = 0.0, T
            set_resource_bounds!(G, i, time_res_id, l_i_time, u_i_time)
        end

        # Build set of arcs A from E´ (two arcs for each edge (i,j))
        for f in F´ # setting the arcs between source, sink, and black vertices
            # source -> i(AFS)
            arc_id = add_arc!(G, v_source, f)
            set_arc_consumption!(G, arc_id, time_res_id, data.G´.V´[f].service_time/2.0)
            set_arc_consumption!(G, arc_id, fuel_res_id, 0.0)
            # i(AFS) -> sink
            arc_id = add_arc!(G, f, v_sink)
            set_arc_consumption!(G, arc_id, time_res_id, data.G´.V´[f].service_time/2.0)
            set_arc_consumption!(G, arc_id, fuel_res_id, 0.0)
        end

        for k in K
            for i in V
                for j in V
                    if !((i, j) in data.E´) && i != j
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

                        # (AFs, AFs)
                        if i in F´ && j in F´ && i != V[1] && j != V[1]
                            set_arc_consumption!(G, arc_id, fuel_res_id, β + 1)
                        # (1, AFs) ou (AFs, 1) 
                        elseif i in F´ && j in F´ && i == V[1] || j == V[1]
                            set_arc_consumption!(G, arc_id, fuel_res_id,  - f(data,ed(i,j)) )
                        # (1, AFs)
                        #elseif i in F´ && j in F´ && i == V[1] && j != V[1]
                        #    set_arc_consumption!(G, arc_id, fuel_res_id,  - f(data,ed(i,j)) )
                        # (AFs, 1)
                        #elseif i in F´ && j in F´ && i != V[1] && j == V[1]
                        #    set_arc_consumption!(G, arc_id, fuel_res_id,  f(data,ed(i,j)) )
                        # (j, i)
                        else
                            set_arc_consumption!(G, arc_id, fuel_res_id,  f(data,ed(i,j)))
                        end
                        #set_arc_consumption!(G, arc_id, fuel_res_id,  f(data,ed(i,j)))
                        
                        #set_arc_consumption!(G, arc_id, time_res_id, t(data,ed(i,j)))
                        set_arc_consumption!(G, arc_id, time_res_id, ((data.G´.V´[i].service_time + data.G´.V´[j].service_time)/2) + t(data,ed(i,j)))
                        
                        #######################################################

                        # add arcs j - > i
                        arc_id = add_arc!(G, j, i)
                        add_arc_var_mapping!(G, arc_id, x[i, j, k])

                        # (AFs, AFs)
                        if i in F´ && j in F´ && i != V[1] && j != V[1]
                            set_arc_consumption!(G, arc_id, fuel_res_id, β + 1)
                        # (1, AFs) ou (AFs, 1) 
                        elseif i in F´ && j in F´ && i == V[1] || j == V[1]
                            set_arc_consumption!(G, arc_id, fuel_res_id,  - f(data,ed(i,j)) )
                        # (1, AFs)
                        #elseif i in F´ && j in F´ && i == V[1] && j != V[1]
                        #    set_arc_consumption!(G, arc_id, fuel_res_id,  - f(data,ed(i,j)) )
                        # (AFs, 1)
                        #elseif i in F´ && j in F´ && i != V[1] && j == V[1]
                        #    set_arc_consumption!(G, arc_id, fuel_res_id,  f(data,ed(i,j)) )
                        # (j, i)
                        else
                            set_arc_consumption!(G, arc_id, fuel_res_id,  f(data,ed(i,j)) )
                        end
                        #set_arc_consumption!(G, arc_id, fuel_res_id,  f(data,ed(i,j)))

                        #set_arc_consumption!(G, arc_id, time_res_id, t(data,ed(i,j)))
                        set_arc_consumption!(G, arc_id, time_res_id, ((data.G´.V´[i].service_time + data.G´.V´[j].service_time)/2) + t(data,ed(i,j)))
                    end
                end
            end
        end
        return G
    end

    # Build the model directed graph G=(V,A)
    # still have to test
    function build_graph2(k::Int64)
        v_source = v_sink = 0
        L = length(F´) 
        #U = length(V) # max and min number of paths is equal to number of AFSs
        #L = length(C)

        #L = 0 
        U = length(C)

        # node ids of G from 0 to |V|
        G = VrpGraph(gvrp, V´, v_source, v_sink, (L, U))
        # resources, R = R_M = {1,2} = {cap_res_id, fuel_res_id}}
        time_res_id = add_resource!(G, main=true)
        fuel_res_id = add_resource!(G, main=true)

        for i in V´
            # l_i, u_i = 0.0, Float64(Q) # accumulated resource consumption interval [l_i, u_i] for the vertex i
            # set_resource_bounds!(G, i, cap_res_id, l_i, u_i)

            l_i_fuel, u_i_fuel = 0.0, β 
            set_resource_bounds!(G, i, fuel_res_id, l_i_fuel, u_i_fuel)
            l_i_time, u_i_time = 0.0, T
            set_resource_bounds!(G, i, time_res_id, l_i_time, u_i_time)
        end

        # Build set of arcs A from E´ (two arcs for each edge (i,j))
        for f in F´ # setting the arcs between source, sink, and black vertices
            # source -> i(AFS)
            arc_id = add_arc!(G, v_source, f)
            set_arc_consumption!(G, arc_id, time_res_id, data.G´.V´[f].service_time/2.0)
            set_arc_consumption!(G, arc_id, fuel_res_id, 0.0)
            # i(AFS) -> sink
            arc_id = add_arc!(G, f, v_sink)
            set_arc_consumption!(G, arc_id, time_res_id, data.G´.V´[f].service_time/2.0)
            set_arc_consumption!(G, arc_id, fuel_res_id, 0.0)
        end

        for i in V
            for j in V
                if !((i, j) in data.E´) && i != j
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

                    # (AFs, AFs)
                    if i in F´ && j in F´ && i != V[1] && j != V[1]
                        set_arc_consumption!(G, arc_id, fuel_res_id, β + 1)
                    # (1, AFs) ou (AFs, 1) 
                    elseif i in F´ && j in F´ && i == V[1] || j == V[1]
                        set_arc_consumption!(G, arc_id, fuel_res_id,  - f(data,ed(i,j)) )
                    # (1, AFs)
                    #elseif i in F´ && j in F´ && i == V[1] && j != V[1]
                    #    set_arc_consumption!(G, arc_id, fuel_res_id,  - f(data,ed(i,j)) )
                    # (AFs, 1)
                    #elseif i in F´ && j in F´ && i != V[1] && j == V[1]
                    #    set_arc_consumption!(G, arc_id, fuel_res_id,  f(data,ed(i,j)) )
                    # (j, i)
                    else
                        set_arc_consumption!(G, arc_id, fuel_res_id,  f(data,ed(i,j)))
                    end
                    #set_arc_consumption!(G, arc_id, fuel_res_id,  f(data,ed(i,j)))
                    
                    #set_arc_consumption!(G, arc_id, time_res_id, t(data,ed(i,j)))
                    set_arc_consumption!(G, arc_id, time_res_id, ((data.G´.V´[i].service_time + data.G´.V´[j].service_time)/2) + t(data,ed(i,j)))
                    
                    #######################################################

                    # add arcs j - > i
                    arc_id = add_arc!(G, j, i)
                    add_arc_var_mapping!(G, arc_id, x[i, j, k])

                    # (AFs, AFs)
                    if i in F´ && j in F´ && i != V[1] && j != V[1]
                        set_arc_consumption!(G, arc_id, fuel_res_id, β + 1)
                    # (1, AFs) ou (AFs, 1) 
                    elseif i in F´ && j in F´ && i == V[1] || j == V[1]
                        set_arc_consumption!(G, arc_id, fuel_res_id,  - f(data,ed(i,j)) )
                    # (1, AFs)
                    #elseif i in F´ && j in F´ && i == V[1] && j != V[1]
                    #    set_arc_consumption!(G, arc_id, fuel_res_id,  - f(data,ed(i,j)) )
                    # (AFs, 1)
                    #elseif i in F´ && j in F´ && i != V[1] && j == V[1]
                    #    set_arc_consumption!(G, arc_id, fuel_res_id,  f(data,ed(i,j)) )
                    # (j, i)
                    else
                        set_arc_consumption!(G, arc_id, fuel_res_id,  f(data,ed(i,j)) )
                    end
                    #set_arc_consumption!(G, arc_id, fuel_res_id,  f(data,ed(i,j)))

                    #set_arc_consumption!(G, arc_id, time_res_id, t(data,ed(i,j)))
                    set_arc_consumption!(G, arc_id, time_res_id, ((data.G´.V´[i].service_time + data.G´.V´[j].service_time)/2) + t(data,ed(i,j)))
                end
            end
        end
        return G
    end

    """
    for i in V
        for j in V
            if i < j && !((i, j) in data.E´) )
            end
        end
    end

    if i < j && !((i, j) in data.E´)
    """

    flag = 1
    if flag == 1
        G = build_graph()
        add_graph!(gvrp, G)

        set_vertex_packing_sets!(gvrp, [[(G, i)] for i in C])
        set_additional_vertex_elementarity_sets!(gvrp, [(G,[f]) for f in data.F])
        define_elementarity_sets_distance_matrix!(gvrp, G, [[ d(data,ed(i, j) ) for i in V] for j in V])

        #add_capacity_cut_separator!(gvrp, [ ([(G, i)], 2.0*floor(data.G´.V´[i].service_time) ) for i in C], 2.0*floor(data.T) )

        # 5.64 ,  com capacity_cut
        # 5.81 , 3.75 , 3.83 , 4.39  , sem capacity_cut
        # quantos 1.0 eu devo acumular até chegar no máximo Q
        # infos data.G´.V´[i].service_time e data.T
        # transformar minutos 0.5 em minutos inteiro
    
    # Trabalhos futuros: paralelizar o processo por rota
    elseif flag == 2
        Graphs = []
        for k in K
            push!(Graphs, build_graph2(k))
        end

        for G in Graphs
          add_graph!(gvrp, G)
        end
        
        #set_vertex_packing_sets!(gvrp, [[(Graphs[k], i) for k in K] for i in V])
        set_vertex_packing_sets!(gvrp, [[(Graphs[k], i) for k in K] for i in C])
        #set_vertex_packing_sets!(gvrp, [[(Graphs[k], i) for i in C] for k in K])

        #add_capacity_cut_separator!(gvrp, [ ([(Graphs[k], i) for k in K], 1.0) for i in C], Float64( length(C) )) #β

        # com 2.84
        # sem cap 3.04 
        [set_additional_vertex_elementarity_sets!(gvrp, [(Graphs[k],[f]) for k in K]) for f in F´] 

        [define_elementarity_sets_distance_matrix!(gvrp, Graphs[k], [[d(data,ed(i, j)) for i in V] for j in V]) for k in K]
    end

    set_branching_priority!(gvrp, "x", 1)
    # "Solution not found" p/ prioridade de branch na var e
    #set_branching_priority!(gvrp, "e", 2)

    function maxflow_mincut_kPathCut_callback()
        println("find sets ...")
        M = 100000
        added_cuts = []

        # for all routes
        g = SparseMaxFlowMinCut.ArcFlow[]
        for i in V
            for j in V 
                if !((i, j) in data.E´)
                    value::Float64 = sum(get_value(gvrp.optimizer, x[i, j, k]) for k in K)
                    #            value::Float64 = get_value(gvrp.optimizer, x[e] * t(data, ee))
                    if value > 0.0001
                        flow_::Int = trunc(floor(value, digits=5) * M)
                        push!(g, SparseMaxFlowMinCut.ArcFlow(i, j, flow_)) # arc i -> j
                        push!(g, SparseMaxFlowMinCut.ArcFlow(j, i, flow_)) # arc j -> i
                    end
                end
            end
        end

        s = V[1]
        for i in V
            for j in V
                maxFlow, flows, cut = SparseMaxFlowMinCut.find_maxflow_mincut(SparseMaxFlowMinCut.Graph(n, g), s, j)
                #if (maxFlow / M) > (T - 0.001) && !in(cut, added_cuts)
                if !((i, j) in data.E´) && (maxFlow / M) < (2 - 0.001) && !in(cut, added_cuts)
                    set1, set2 = [], []
                    [cut[i] == 1 ? push!(set1, i) : push!(set2, i) for i in 1:n]

                    lhs_vars = [x[i, j, k] for i in set2 for j in set1 for k in K]
                    lhs_coeff = [1.0 for i in set2 for j in set1 for k in K]

                    #add_dynamic_constr!(gvrp.optimizer, lhs_vars, lhs_coeff, >=, 2.0 * floor(ceil(sum(data.G´.V´[i].service_time for i in (c in set1 ? set1 : set2) if i in C)/T)), "mincut")
                    add_dynamic_constr!(gvrp.optimizer, lhs_vars, lhs_coeff, >=, 2.0, "mincut")

                    push!(added_cuts, cut)
                end
            end
        end

        # for each route
        for k in K
          s = V[1]
          visited = [false for i in V]
          for c in V
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
                if i != j && !((i, j) in data.E´) && ( get_value(gvrp.optimizer, x[i, j, k]) > 0.0001 || get_value(gvrp.optimizer, x[j, i, k]) > 0.0001 ) && !visited[j]
                  visited[j] = true
                  push!(q, j)
                  push!(comp, j)
                end
              end
            end
            #println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
            
            if !(s in comp) && length(comp) > 1
              println("Route $k S: ", comp)
              lhs_vars = [x[i, j, k] for i in comp for j in V if !(j in comp && !((i, j) in data.E´) )]
              lhs_coeff = [1.0 for i in comp for j in V if !(j in comp && !((i, j) in data.E´) )]
              for i in comp
                if i in V
                  for j in comp
                    if i != j && !((i, j) in data.E´)
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
        
        if length(added_cuts) > 0 
          println(">>>>> Add min cuts : ", length(added_cuts), " cut(s) added") 
        end
    end
    add_cut_callback!(gvrp, maxflow_mincut_kPathCut_callback, "mincut")
    return (gvrp, x)
end
