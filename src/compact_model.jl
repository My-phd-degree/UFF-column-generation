function compact_model(data::DataGVRP)
  E = edges(data) # set of edges of the input graph G′
  β = data.β
  C = data.C # Set of customers vertices
  depot_id = data.depot_id
  C₀ = vcat([depot_id], C) # customers and depot
  F = deepcopy(data.F) # Set of AFSs vertices
  T = data.T # Time limit 
  F₀ = deepcopy(vcat([depot_id], F)) # AFSs and depot
  n = nb_vertices(data)
  fuel = f
  gvrp_afs_tree = data.gvrp_afs_tree
  #updating graph
  V′ = data.G′.V′
    #new afss
  afss_pairs = Dict{Int64, Tuple{Int64, Int64}}()
  for f in F 
    for r in F 
      if f < r &&
        !(depot_id in getAFSTreePath(f, r, gvrp_afs_tree))
        mayBeUsed = false 
        for i in C₀
          for j in C₀
            for f´ in F₀
              for r´ in F₀
                if ed(f´, i) in E && 
                  ed(i, f) in E &&
                  ed(r, j) in E &&
                  ed(j, r´) in E && 
                  gvrp_afs_tree.times[f´] + t(data, ed(f´, i)) + t(data, ed(i, f)) + gvrp_afs_tree.pairTimes[(f, r)] + t(data, ed(r, j)) + t(data, ed(j, r´)) + gvrp_afs_tree.times[r´] <= T && 
                  fuel(data, ed(f´, i)) + fuel(data, ed(i, f)) <= β && 
                  fuel(data, ed(r, j)) + fuel(data, ed(j, r´)) <= β
                  mayBeUsed = true
                  break
                end
              end
              if mayBeUsed
                break
              end
            end
            if mayBeUsed
              break
            end
          end
          if mayBeUsed
            break
          end
        end
        if mayBeUsed
          n = n + 1
          afss_pairs[n] = (f, r) 
          push!(data.G′.V′, Vertex(n, 0.0, 0.0, sum([V′[s].service_time for s in getAFSTreePath(f, r, gvrp_afs_tree)]), gvrp_afs_tree.pairCosts[(f, r)]))
          push!(data.F, n)
          n = n + 1
          afss_pairs[n] = (r, f) 
          push!(data.G′.V′, Vertex(n, 0.0, 0.0, sum([V′[s].service_time for s in getAFSTreePath(r, f, gvrp_afs_tree)]), gvrp_afs_tree.pairCosts[(r, f)]))
          push!(data.F, n)
          #new edges
          for i in C₀
            if ed(i, f) in E
              push!(data.G′.E, (i, n - 1))
              push!(data.G′.E, (n, i))
              data.G′.cost[(i, n - 1)] = data.G′.cost[(n, i)] = EUC_dist(V′[i], V′[f])
            end
            if ed(i, r) in E
              push!(data.G′.E, (i, n))
              push!(data.G′.E, (n - 1, i))
              data.G′.cost[(i, n)] = data.G′.cost[(n - 1, i)] = EUC_dist(V′[i], V′[r])
            end
          end
        end
      end
    end
  end
  data.non_consec = true
  return afss_pairs 
end

