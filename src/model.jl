function build_model(data::DataBWTSP)

    E = edges(data) # set of edges of the input graph G′
    n = nb_vertices(data)
    V = [i for i in 1:n] # set of vertices of the input graph G′
    V′ = [i for i in 0:n] # V ⋃ {0}, where 0 is a dummy vertex

    Q = data.Q #  Maximum number of white between two consecutive black
    D = data.D
    W = data.White # Set of white vertices
    B = data.Black # Set of black vertices
    q = data.q
    service_time = data.service_time

    ed(i, j) = i < j ? (i, j) : (j, i)

    # Formulation
    bwtsp = VrpModel()
    @variable(bwtsp.formulation, x[e in E], Int)
    @variable(bwtsp.formulation, y[i in B] >= 0, Int)
    @objective(bwtsp.formulation, Min, sum(c(data, e) * x[e] for e in E))
    @constraint(bwtsp.formulation, deg[i in W], sum(x[e] for e in δ(data, i)) == 2.0)
    @constraint(bwtsp.formulation, hotel_deg[i in B], sum(x[e] for e in δ(data, i)) == 2*y[i])
    @constraint(bwtsp.formulation, y[B[1]] >= 1)
    @constraint(bwtsp.formulation, sum(y[i] for i in B) == q)
    #println(bwtsp.formulation)

    # Build the model directed graph G=(V,A)
    function build_graph()

        v_source = v_sink = 0
        L = U = q # max and min number of paths is equal to number of black nodes

        # node ids of G from 0 to |V|
        G = VrpGraph(bwtsp, V′, v_source, v_sink, (L, U))
        # resourves, R = R_M = {1,2} = {cap_res_id, dist_res_id}}
        # cap_res_id = add_resource!(G, main=true)
        dist_res_id = add_resource!(G, main=true)
        for i in V′
            # l_i, u_i = 0.0, Float64(Q) # accumulated resource consumption interval [l_i, u_i] for the vertex i
            # set_resource_bounds!(G, i, cap_res_id, l_i, u_i)

            l_i, u_i = 0.0, Float64(D)
            set_resource_bounds!(G, i, dist_res_id, l_i, u_i)
        end

        # Build set of arcs A from E′ (two arcs for each edge (i,j))
        for i in B # setting the arcs between source, sink, and black vertices
            # source -> i(black)
            arc_id = add_arc!(G, v_source, i)
            # set_arc_consumption!(G, arc_id, cap_res_id, 0.0)
            set_arc_consumption!(G, arc_id, dist_res_id, 0.0)
            # i(black) -> sink
            arc_id = add_arc!(G, i, v_sink)
            # set_arc_consumption!(G, arc_id, cap_res_id, 0.0)
            set_arc_consumption!(G, arc_id, dist_res_id, 0.0)
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
            # set_arc_consumption!(G, arc_id, cap_res_id, q_1)
            st = 0
            if i in W && j in W
               st = service_time
            elseif i in W || j in W
               st = service_time/2
            end
             
            #st = j in W ? service_time : 0.0
            set_arc_consumption!(G, arc_id, dist_res_id, d(data,(i,j)) + st)
             # add arcs j - > i
            arc_id = add_arc!(G, j, i)
            add_arc_var_mapping!(G, arc_id, x[(i, j)])
            # set_arc_consumption!(G, arc_id, cap_res_id, q_1)
            #st = i in W ? service_time : 0.0
            set_arc_consumption!(G, arc_id, dist_res_id, d(data,(i,j)) + st)
        end
        return G
    end

    G = build_graph()
    add_graph!(bwtsp, G)
    #println(G)

    set_vertex_packing_sets!(bwtsp, [[(G, i)] for i in W])
    # set_additional_vertex_elementarity_sets!(bwtsp, [(G,[i]) for i in B])

    define_elementarity_sets_distance_matrix!(bwtsp, G, [[d(data,ed(i, j)) for i in W] for j in W])

    # add_capacity_cut_separator!(bwtsp, [ ([(G, i)], 1.0) for i in W], Float64(Q))

    set_branching_priority!(bwtsp, "x", 1)
    set_branching_priority!(bwtsp, "y", 1)

    function maxflow_mincut_callback()
        M = 100000
        g = SparseMaxFlowMinCut.ArcFlow[]
        for (i, j) in E
            e = (i, j)
            value::Float64 = get_value(bwtsp.optimizer, x[e])
            if value > 0.0001
                flow_::Int = trunc(floor(value, digits=5) * M)
                push!(g, SparseMaxFlowMinCut.ArcFlow(i, j, flow_)) # arc i -> j
                push!(g, SparseMaxFlowMinCut.ArcFlow(j, i, flow_)) # arc j -> i
            end
        end

        added_cuts = []
        s = B[1]
        for t in W
            maxFlow, flows, cut = SparseMaxFlowMinCut.find_maxflow_mincut(SparseMaxFlowMinCut.Graph(n, g), s, t)
            if (maxFlow / M) < (2 - 0.001) && !in(cut, added_cuts)
                set1, set2 = [], []
                [cut[i] == 1 ? push!(set1, i) : push!(set2, i) for i in 1:n]

                add_dynamic_constr!(bwtsp.optimizer, [x[ed(i, j)] for i in set1 for j in set2], [1.0 for i in set1 for j in set2], >=, 2.0, "mincut")
                push!(added_cuts, cut)
            end
        end
        if length(added_cuts) > 0 println(">>>>> Add min cuts : ", length(added_cuts), " cut(s) added") end
    end
    add_cut_callback!(bwtsp, maxflow_mincut_callback, "mincut")


    return (bwtsp, x)
end
