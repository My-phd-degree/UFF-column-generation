include("gvrp_afs_tree.jl")

using CPLEX

ed(i, j) = i < j ? (i, j) : (j, i)

function build_model_compact_with_arcs(data::DataGVRP)
  # uniderected GVRP instance
  E = edges(data) # set of edges of the input graph G′
  n = nb_vertices(data)
  V = [i for i in 1:n] # set of vertices of the input graph G′
  V′ = [i for i in 0:n] # V ⋃ {0}, where 0 is a dummy vertex
  depot_id = data.depot_id
  β = data.β
  C = data.C # Set of customers vertices
  C₀ = vcat([depot_id], C) # Set of customers with depot
  F = data.F # Set of AFSs vertices
  F₀ = vcat([depot_id], F) # Set of AFSs with depot
  T = data.T # Set of AFSs vertices
  gvrp_afs_tree = data.gvrp_afs_tree
  fuel = f
  # directed GVRP instance
  D = InputGraph(
    deepcopy(data.G′.V′), 
    vcat([(i, j) for (i, j) in E if !in(i, F) || !in(j, F)], [(j, i) for (i, j) in E if !in(i, F) || !in(j, F)]), 
    merge(Dict{Tuple{Int64, Int64},Float64}((i, j) => v for ((i, j), v) in data.G′.cost if !in(i, F) || !in(j, F)), Dict{Tuple{Int64, Int64},Float64}((j, i) => v for ((i, j), v) in data.G′.cost if !in(i, F) || !in(j, F)))
  )
  directedData = DataGVRP(D, depot_id, data.coord, data.round, true, deepcopy(F), C, β, T, data.ρ, data.ε, data.reduced_graph, data.gvrp_afs_tree)
  A = D.E
  #create compact AFSs paths
  afss_pairs = Dict{Int64, Tuple{Int64, Int64}}()
  for f in F 
    for r in F 
      if f < r &&
        !in(depot_id, getAFSTreePath(f, r, gvrp_afs_tree))
        mayBeUsed = false 
        for i in C₀
          for j in C₀
            for f´ in F₀
              for r´ in F₀
                if (f´, i) in A && 
                  (i, f) in A &&
                  (r, j) in A &&
                  (j, r´) in A && 
                  gvrp_afs_tree.times[f´] + t(directedData, (f´, i)) + t(directedData, (i, f)) + gvrp_afs_tree.pairTimes[(f, r)] + t(directedData, (r, j)) + t(directedData, (j, r´)) + gvrp_afs_tree.times[r´] <= T && 
                  fuel(directedData, (f´, i)) + fuel(directedData, (i, f)) <= β && 
                  fuel(directedData, (r, j)) + fuel(directedData, (j, r´)) <= β
                  mayBeUsed = true
                  break
                end
              end
              mayBeUsed && break
            end
            mayBeUsed && break
          end
          mayBeUsed && break
        end
        if mayBeUsed
          n = n + 1
          afss_pairs[n] = (f, r) 
          push!(directedData.G′.V′, Vertex(n, 0.0, 0.0, sum([directedData.G′.V′[s].service_time for s in getAFSTreePath(f, r, gvrp_afs_tree)]), gvrp_afs_tree.pairCosts[(f, r)]))
          push!(directedData.F, n)
          n = n + 1
          afss_pairs[n] = (r, f) 
          push!(directedData.G′.V′, Vertex(n, 0.0, 0.0, sum([directedData.G′.V′[s].service_time for s in getAFSTreePath(r, f, gvrp_afs_tree)]), gvrp_afs_tree.pairCosts[(r, f)]))
          push!(directedData.F, n)
          #new edges
          for i in C₀
            if (i, f) in A
              push!(A, (i, n - 1))
              push!(A, (n, i))
              directedData.G′.cost[(i, n - 1)] = directedData.G′.cost[(n, i)] = EUC_dist(directedData.G′.V′[i], directedData.G′.V′[f])
            end
            if (i, r) in A
              push!(A, (i, n))
              push!(A, (n - 1, i))
              directedData.G′.cost[(i, n)] = directedData.G′.cost[(n - 1, i)] = EUC_dist(directedData.G′.V′[i], directedData.G′.V′[r])
            end
          end
        end
      end
    end
  end
  F = directedData.F
  V = [i for i in 1:n] # set of vertices of the input graph G′
  F´ = copy(F) 
  pushfirst!(F´, depot_id)
  
  # Formulation
  gvrp = VrpModel()
  @variable(gvrp.formulation, 1 >= x[a in A] >= 0, Int)
#  @variable(gvrp.formulation, 2 * length(C) >= y[i in F´] >= 0, Int)
@objective(gvrp.formulation, Min, sum((d(directedData, (i, j)) + directedData.G′.V′[i].weight) * x[(i, j)] for (i, j) in A))
  @constraint(gvrp.formulation, deg[i in V], sum(x[a] for a in δ⁻(directedData, i)) == sum(x[a] for a in δ⁺(directedData, i)))
  @constraint(gvrp.formulation, degCustomers[i in C], sum(x[a] for a in δ⁻(directedData, i)) == 1)
#  @constraint(gvrp.formulation, hotel_deg[i in F´], sum(x[e] for e in δ⁻(directedData, i)) == y[i])
  @constraint(gvrp.formulation, no_edge_between_afss[f in F], sum(x[(f, r)] for r in F if (f, r) in A) == 0.0)
#  @constraint(gvrp.formulation, y[F´[1]] >= 1)

#  routes = [
#           [0, 12, 5, 10, 39, 10, 49, 9, 50, 38, 50, 21, 2, 0],
#           [0, 6, 23, 48, 23, 7, 23, 24, 23, 0],
#           [0, 32, 31, 8, 31, 28, 31, 26, 31, 1, 32, 27, 0],
#           [0, 18, 13, 41, 40, 19, 42, 19, 4, 18, 0],
#           [0, 32, 22, 2, 35, 3, 35, 36, 35, 20, 2, 29, 2, 0, ],
#           [0, 12, 37, 12, 17, 12, 47, 18, 25, 18, 14, 18, 0],
#           [0, 46, 12, 45, 15, 44, 45, 33, 10, 30, 34, 50, 16, 11, 32, 0],
#          ]

if @isdefined routes
  print
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
    for (i, j) in A
      if i in F && j in F
         continue
       end
      # add arcs i - > j
      arc_id = nothing
      if j in F
        arc_id = add_arc!(G, i, afssIn[j])
      elseif i in F
        arc_id = add_arc!(G, afssOut[i], j)
      else
        arc_id = add_arc!(G, i, j)
      end
      add_arc_var_mapping!(G, arc_id, x[(i, j)])
      set_arc_consumption!(G, arc_id, fuel_res_id, f(directedData,(i,j)))
      set_arc_consumption!(G, arc_id, time_res_id, t(directedData,(i,j)))
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

#  set_branching_priority!(gvrp, "y", 1)
  set_branching_priority!(gvrp, "x", 1)

  function maxflow_mincut_callback()
    M = 100000
    # for all routes
    g = SparseMaxFlowMinCut.ArcFlow[]
    for (i, j) in A
      a = (i, j)
      value::Float64 = sum(get_value(gvrp.optimizer, x[a]))
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
        lhs_vars = [x[(i, j)] for i in set2 for j in set1 if (i, j) in A]
        lhs_coeff = [1.0 for i in set2 for j in set1 if (i, j) in A]
        setIn = c in set1 ? set1 : set2
        add_dynamic_constr!(gvrp.optimizer, lhs_vars, lhs_coeff, >=, 1.0, "mincut")

        push!(added_cuts, cut)
      end
    end
    if length(added_cuts) > 0 
      println(">>>>> Add min cuts : ", length(added_cuts), " cut(s) added") 
    end
  end
  add_cut_callback!(gvrp, maxflow_mincut_callback, "mincut")

#  return (gvrp, x, y)
  return (directedData, gvrp, x, afss_pairs)
end
