function build_model(data::DataGVRP)

    E = edges(data) # set of edges of the input graph G′
    n = nb_vertices(data)
    V = [i for i in 1:n] # set of vertices of the input graph G′
    V′ = [i for i in 0:n] # V ⋃ {0}, where 0 is a dummy vertex

    β = data.β
    C = data.C # Set of customers vertices
    F = data.F # Set of AFSs vertices
    T = data.T # Set of AFSs vertices

    ed(i, j) = i < j ? (i, j) : (j, i)

    # Formulation
    gvrp = VrpModel()
    @variables(gvrp.formulation, begin
                 x[e in E], Int
                 y[i in F] >= 0, Int
               end)
    @objective(gvrp.formulation, Min, sum(d(data, e) * x[e] for e in E))
    @constraints(gvrp.formulation, begin
                  deg[i in C], sum(x[e] for e in δ(data, i)) == 2.0
                  afs_deg[i in F], sum(x[e] for e in δ(data, i)) == 2*y[i]
                  y[F[1]] >= 1
                end)
    #println(gvrp.formulation)

    # Build the model directed graph G=(V,A)
    function build_graph()

        v_source = v_sink = 0
        L = U = length(F) # max and min number of paths is equal to number of AFSs

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
            # resource comsuption for the R = 1 
            # q_1 = 0.5
            # if (i in W) && (j in W)
            #     q_1 = 1.0
            # elseif (i in B) && (j in B)
            #     q_1 = Q
            # end

            # add arcs i - > j
            arc_id = add_arc!(G, i, j)
            add_arc_var_mapping!(G, arc_id, x[(i, j)])
             
            set_arc_consumption!(G, arc_id, fuel_res_id, f(data,(i,j)))
            set_arc_consumption!(G, arc_id, time_res_id, ((data.G′.V′[i].service_time + data.G′.V′[j].service_time)/2) + t(data,(i,j)))
            # add arcs j - > i
            arc_id = add_arc!(G, j, i)
            add_arc_var_mapping!(G, arc_id, x[(i, j)])

            set_arc_consumption!(G, arc_id, fuel_res_id, f(data,(i,j)))
            set_arc_consumption!(G, arc_id, time_res_id, ((data.G′.V′[i].service_time + data.G′.V′[j].service_time)/2) + t(data,(i,j)))
        end
        return G
    end

    G = build_graph()
    add_graph!(gvrp, G)
    #println(G)

    set_vertex_packing_sets!(gvrp, [[(G, i)] for i in C])
    # set_additional_vertex_elementarity_sets!(gvrp, [(G,[i]) for i in B])

    define_elementarity_sets_distance_matrix!(gvrp, G, [[d(data,ed(i, j)) for i in C] for j in C])

    # add_capacity_cut_separator!(gvrp, [ ([(G, i)], 1.0) for i in W], Float64(Q))

    set_branching_priority!(gvrp, "x", 1)
    set_branching_priority!(gvrp, "y", 1)

    function maxflow_mincut_callback()
        M = 100000
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
        for t in C
            maxFlow, flows, cut = SparseMaxFlowMinCut.find_maxflow_mincut(SparseMaxFlowMinCut.Graph(n, g), s, t)
            if (maxFlow / M) < (2 - 0.001) && !in(cut, added_cuts)
                set1, set2 = [], []
                [cut[i] == 1 ? push!(set1, i) : push!(set2, i) for i in 1:n]

                setIn = s in set1 ? set2 : set1
                setOut = s in set1 ? set1 : set2

                lhs_vars = vcat([x[ed(i, j)] for i in setIn for j in setOut], [x[e] for e in E if e[1] in setIn && e[2] in setIn])
                lhs_coeff = vcat([1.0 for i in setIn for j in setOut], [-2.0 * (d(data, e)/data.ε)/T for e in E if e[1] in setIn && e[2] in setIn])

                rhs = sum(data.G′.V′[i].service_time for i in setIn)/T

                for i in setIn
                  print(i, ", ")
                end
                print(": $rhs \n")

                add_dynamic_constr!(gvrp.optimizer, lhs_vars, lhs_coeff, >=, 2.0 * rhs, "mincut")

                push!(added_cuts, cut)
            end
        end
        if length(added_cuts) > 0 
          println(">>>>> Add min cuts : ", length(added_cuts), " cut(s) added") 
        end
    end
    add_cut_callback!(gvrp, maxflow_mincut_callback, "mincut")

    return (gvrp, x)
end
