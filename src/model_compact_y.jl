include("gvrp_afs_tree.jl")

using CPLEX

ed(i, j) = i < j ? (i, j) : (j, i)

function build_model_compact_y(data::DataGVRP)

  E = edges(data) # set of edges of the input graph G′
  n = nb_vertices(data)
  V = [i for i in 1:n] # set of vertices of the input graph G′
  V′ = [i for i in 0:n] # V ⋃ {0}, where 0 is a dummy vertex
  depot_id = data.depot_id

  β = data.β
  C = data.C # Set of customers vertices
  C₀ = vcat(data.depot_id, C) # Set of customers and depot
  EC₀ = [(i, j) for (i, j) in E if i in C₀ && j in C₀]
  E₀ = [(i, j) for (i, j) in E if i == depot_id && j in C]
  EC = [(i, j) for (i, j) in E if i in C && j in C]
  F = data.F # Set of AFSs vertices
  F₀ = vcat([depot_id], data.F) # Set of AFSs vertices with depot
  F′ = vcat([depot_id], F)
  T = data.T # Set of AFSs vertices
  gvrp_afs_tree = data.gvrp_afs_tree
  fuel = f
  P = Dict{Tuple{Int64, Int64, Int64, Int64},Float64}()
  function δ′(P::Dict{Tuple{Int64, Int64, Int64, Int64},Float64}, i::Integer)
    return [p for p in keys(P) if p[1] == i || p[4] == i]
  end
  function δ′′(P::Dict{Tuple{Int64, Int64, Int64, Int64},Float64}, i::Integer, j::Integer)
    return [p for p in keys(P) if (p[1] == i && p[4] == j) || (p[1] == j && p[4] == i)]
  end
  # calculating P
  for i in C₀
    for j in C₀
      if i < j
        for f in F
          (!(ed(i, f) in E) || fuel(data, ed(i, f)) > β) && continue
          for r in F
            (!(ed(r, j) in E) || fuel(data, ed(r, j)) > β) && continue

            for f′ in F₀
              ((!(ed(f′, i) in E) && f′ != depot_id) || fuel(data, ed(f′, i)) + fuel(data, ed(i, f)) > β) && continue
              mayBeUsed = false
              for r′ in F₀
                if ed(j, r′) in E &&
                  !(depot_id in getAFSTreePath(f, r, gvrp_afs_tree)) &&
                  fuel(data, ed(r, j)) + fuel(data, ed(j, r′)) <= β &&
                  gvrp_afs_tree.times[f′] + t(data, ed(f′, i)) + t(data, ed(i, f)) + gvrp_afs_tree.pairTimes[(f, r)] + t(data, ed(r, j)) + t(data, ed(j, r′)) + gvrp_afs_tree.times[r′] <= T
                  P[(i, f, r, j)] = d(data, ed(i, f)) + gvrp_afs_tree.pairCosts[(f, r)] + d(data, ed(r, j))
                  mayBeUsed = true
                  break
                end
              end
              mayBeUsed && break
            end
          end
        end
      end
    end
  end
  # Formulation
  gvrp = VrpModel()
  @variable(gvrp.formulation, x[e in EC₀], Int)
  @variable(gvrp.formulation, y[p in keys(P)] >= 0, Int)
  @objective(gvrp.formulation, Min, sum((d(data, e) * x[e] for e in EC₀)) + sum(y[p] * c for (p, c) in P))
  @constraint(gvrp.formulation, deg[i in C], sum(x[(j, k)] for (j, k) in δ(data, i) if j in C₀ && k in C₀) + sum(y[p] for p in δ′(P, i)) == 2.0)
  @constraint(gvrp.formulation, customers_edges[e in EC], x[e] <= 1)
  @constraint(gvrp.formulation, depot_edges[e in E₀], x[e] <= 2)

# routes = [
#           [0, 2, 7, 5, 21, 20, 21, 14, 17, 21, 9, 21, 13, 2, 0],
#           [0, 16, 10, 2, 8, 18, 3, 18, 19, 18, 2, 6, 0],
#           [0, 16, 15, 12, 15, 4, 15, 11, 15, 1, 16, 0],
#          ]

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
    xWeights = Dict{Tuple{Int64, Int64},Int64}()
    yWeights = Dict{Tuple{Int64, Int64, Int64, Int64},Int64}()
    for route in routes
      i = 2
      while i <= length(route)
        j, k = (ids[route[i - 1]], ids[route[i]])
        e = ed(j, k)
        if e in EC₀
          #x
          xWeights[e] = haskey(xWeights, e) ? xWeights[e] + 1 : 1
        else
          #y
          while true 
            i = i + 1
            ids[route[i]] in C₀ && break
          end
          p = (j, k, ids[route[i - 1]], ids[route[i]])
          p′ = (ids_[p[1]], ids_[p[2]], ids_[p[3]], ids_[p[4]])
          if !(p in keys(P)) 
            p = reverse(p)
          end
          !(p in keys(P)) && error("Variable y[$p = $(p′)] does not exist")
          if haskey(yWeights, p)
            yWeights[p] = yWeights[p] + 1
          else
            yWeights[p] = 1
          end
        end
        i = i + 1
      end 
    end
    for e in keys(xWeights)
      println(data.G′.V′[e[1]].id_vertex, " ", data.G′.V′[e[2]].id_vertex)
      println("Edge $e; Edge weight $(xWeights[e]); EdgeId $(ed(ids_[e[1]], ids_[e[2]]))")
      @constraint(gvrp.formulation, x[e] == xWeights[e])
    end
    for p in keys(yWeights)
      println(data.G′.V′[p[1]].id_vertex, " ", data.G′.V′[p[2]].id_vertex, " ", data.G′.V′[p[3]].id_vertex, " ", data.G′.V′[p[4]].id_vertex)
      println("Path $p; Path weight $(yWeights[p]); PathId $((ids_[p[1]], ids_[p[2]], ids_[p[3]], ids_[p[4]]))")
      @constraint(gvrp.formulation, y[p] == yWeights[p])
    end
  end
  # remove dominated paths
  for i in C₀
    for j in C₀
      if i < j
        P′ = δ′′(P, i, j)
        for p in P′
          (i, f, r, j) = p[1] == i ? p : reverse(p)
          for p′ in P′ 
            (i, f′, r′, j) = p′[1] == i ? p′ : reverse(p′)
            # if p′ dominates p
            if p != p′ && 
              fuel(data, ed(i, f′)) <= fuel(data, ed(i, f)) && 
              fuel(data, ed(r′, j)) <= fuel(data, ed(r, j)) && 
              gvrp_afs_tree.pairCosts[(f′, r′)] <= gvrp_afs_tree.pairCosts[(f, r)]
              delete!(P, p)
              break
            end
          end
        end
      end
    end
  end
  # Build the model directed graph G=(V,A)
  function build_graph()
    v_source = v_sink = depot_id
    L, U = 0, length(C) # max and min number of paths is equal to number of AFSs

    V_ = deepcopy(C₀)

    #AFSs vertex split
    new_id = length(V) + 1
    #going
    firstAfssPaths = Dict{Tuple{Int64, Int64, Int64, Int64}, Int64}()
    lastAfssPaths = Dict{Tuple{Int64, Int64, Int64, Int64}, Int64}()
    #returning
    firstAfssPaths′ = Dict{Tuple{Int64, Int64, Int64, Int64}, Int64}()
    lastAfssPaths′ = Dict{Tuple{Int64, Int64, Int64, Int64}, Int64}()
    for p in keys(P)
      #going
      push!(V_, new_id)
      firstAfssPaths[p] = new_id
      push!(V_, new_id + 1)
      lastAfssPaths[p] = new_id + 1
      push!(V_, new_id + 2)
      #returning
      firstAfssPaths′[p] = new_id + 2
      push!(V_, new_id + 3)
      lastAfssPaths′[p] = new_id + 3
      new_id = new_id + 4
    end
    # node ids of G from 0 to |V|
    G = VrpGraph(gvrp, V_, v_source, v_sink, (L, U))
    # resourves, R = R_M = {1,2} = {cap_res_id, fuel_res_id}}
    l_i_time, u_i_time = 0.0, T
    time_res_id = add_resource!(G, main=true)
    l_i_fuel, u_i_fuel = 0.0, β
    fuel_res_id = add_resource!(G, main=true)
    for p in keys(P)
      #going
      set_resource_bounds!(G, firstAfssPaths[p], fuel_res_id, l_i_fuel, u_i_fuel)
      set_resource_bounds!(G, firstAfssPaths[p], time_res_id, l_i_time, u_i_time)
      set_resource_bounds!(G, lastAfssPaths[p], fuel_res_id, l_i_fuel, u_i_fuel)
      set_resource_bounds!(G, lastAfssPaths[p], time_res_id, l_i_time, u_i_time)
      #returning
      set_resource_bounds!(G, firstAfssPaths′[p], fuel_res_id, l_i_fuel, u_i_fuel)
      set_resource_bounds!(G, firstAfssPaths′[p], time_res_id, l_i_time, u_i_time)
      set_resource_bounds!(G, lastAfssPaths′[p], fuel_res_id, l_i_fuel, u_i_fuel)
      set_resource_bounds!(G, lastAfssPaths′[p], time_res_id, l_i_time, u_i_time)
    end
    for i in C₀
      set_resource_bounds!(G, i, fuel_res_id, l_i_fuel, u_i_fuel)
      set_resource_bounds!(G, i, time_res_id, l_i_time, u_i_time)
    end
    #edges
    for (i, j) in EC₀
      # add arcs i - > j
      arc_id = add_arc!(G, i, j)
      add_arc_var_mapping!(G, arc_id, x[(i, j)])
      set_arc_consumption!(G, arc_id, fuel_res_id, f(data,(i,j)))
      set_arc_consumption!(G, arc_id, time_res_id, t(data,(i,j)))
      # add arcs j - > i
      arc_id = add_arc!(G, j, i)
      add_arc_var_mapping!(G, arc_id, x[(i, j)])
      set_arc_consumption!(G, arc_id, fuel_res_id, f(data,(i,j)))
      set_arc_consumption!(G, arc_id, time_res_id, t(data,(i,j)))
    end
    for (i, f, r, j) in keys(P)
      p = (i, f, r, j)
      # (i, f, r, j)
      # add arcs i - > f
      arc_id = add_arc!(G, i, firstAfssPaths[p])
      set_arc_consumption!(G, arc_id, fuel_res_id, fuel(data,ed(i,f)))
      set_arc_consumption!(G, arc_id, time_res_id, t(data,ed(i,f)))
      # add arcs f - > ... - > r
      p = (i, f, r, j)
      arc_id = add_arc!(G, firstAfssPaths[p], lastAfssPaths[p])
      add_arc_var_mapping!(G, arc_id, y[p])
      set_arc_consumption!(G, arc_id, fuel_res_id, - β)
      set_arc_consumption!(G, arc_id, time_res_id, gvrp_afs_tree.pairTimes[(f, r)])
      # add arcs r - > j
      arc_id = add_arc!(G, lastAfssPaths[p], j)
      set_arc_consumption!(G, arc_id, fuel_res_id, fuel(data,ed(r,j)))
      set_arc_consumption!(G, arc_id, time_res_id, t(data,ed(r,j)))
      # (j, r, f, i)
      # add arcs j - > r
      arc_id = add_arc!(G, j, firstAfssPaths′[p])
      set_arc_consumption!(G, arc_id, fuel_res_id, fuel(data,ed(j,r)))
      set_arc_consumption!(G, arc_id, time_res_id, t(data,ed(j,r)))
      # add arcs r - > ... - > f
      p = (i, f, r, j)
      arc_id = add_arc!(G, firstAfssPaths′[p], lastAfssPaths′[p])
      add_arc_var_mapping!(G, arc_id, y[p])
      set_arc_consumption!(G, arc_id, fuel_res_id, - β)
      set_arc_consumption!(G, arc_id, time_res_id, gvrp_afs_tree.pairTimes[(r, f)])
      # add arcs f - > i
      arc_id = add_arc!(G, lastAfssPaths′[p], i)
      set_arc_consumption!(G, arc_id, fuel_res_id, fuel(data,ed(f,i)))
      set_arc_consumption!(G, arc_id, time_res_id, t(data,ed(f,i)))
    end
    return G
  end

  G = build_graph()
  add_graph!(gvrp, G)

  set_vertex_packing_sets!(gvrp, [[(G, i)] for i in C])

  define_elementarity_sets_distance_matrix!(gvrp, G, [[d(data, ed(i, j)) for i in C] for j in C])

  set_branching_priority!(gvrp, "y", 1)
  set_branching_priority!(gvrp, "x", 1)

  function maxflow_mincut_callback()
    M = 100000
    # for all routes
    g = SparseMaxFlowMinCut.ArcFlow[]
    for i in C₀
      for j in C₀
        if i < j
          value::Float64 = sum(get_value(gvrp.optimizer, y[p]) for p in δ′′(P, i, j)) + (ed(i, j) in EC₀ ? get_value(gvrp.optimizer, x[ed(i, j)]) : 0.0)
          if value > 0.0001
            flow_::Int = trunc(floor(value, digits=5) * M)
            push!(g, SparseMaxFlowMinCut.ArcFlow(i, j, flow_)) # arc i -> j
            push!(g, SparseMaxFlowMinCut.ArcFlow(j, i, flow_)) # arc j -> i
          end
        end
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
        E′ = [ed(i, j) for i in set2 for j in set1 if ed(i, j) in EC₀] 
        P′ = [p for i in set2 for j in set1 for p in δ′′(P, i, j)]
        lhs_vars = vcat([x[e] for e in E′], [y[p] for p in P′])
        lhs_coeff = [1.0 for i in 1:(length(E′) + length(P′))]
        add_dynamic_constr!(gvrp.optimizer, lhs_vars, lhs_coeff, >=, 2.0, "mincut")
        push!(added_cuts, cut)
      end
    end
    if length(added_cuts) > 0 
      println(">>>>> Add min cuts : ", length(added_cuts), " cut(s) added") 
    end
  end
  add_cut_callback!(gvrp, maxflow_mincut_callback, "mincut")


  return (gvrp, P, x, y)
end
