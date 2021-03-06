include("gvrp_afs_tree.jl")

function calculateGVRPReducedGraphTime(data::DataGVRP, gvrp_afs_tree::GVRP_AFS_Tree) 
  times = Dict{Tuple{Int64,Int64}, Float64}()
  depot = data.depot_id
  fuel = f
  C₀ = vcat([depot], data.C)
  F₀ = vcat([depot], data.F)
  for i in C₀
    for j in C₀
      if i < j
        minTime = typemax(Float64)
        for f in F₀
          if !(ed(f, i) in data.G′.E) && i != depot
            continue
          end
          for r in F₀
            if !(ed(j, r) in data.G′.E)
              continue
            end
            # P^{G[F_0]}_{0f} \cup (i, j) \cup P^{G[F_0]}_{r0}
            time = t(data, ed(i, j))
            if ed(i, j) in data.G′.E && 
              fuel(data, ed(f, i)) + fuel(data, ed(i, j)) + fuel(data, ed(j, r)) <= data.β && 
              gvrp_afs_tree.times[f] + t(data, ed(f, i)) + time + t(data, ed(j, r)) + gvrp_afs_tree.times[r] <= data.T
              minTime = time < minTime ? time : minTime
            else
              for f′ in F₀
                if !(ed(i, f′) in data.G′.E) && i != depot
                  continue
                end
                time = t(data, ed(i, f′)) + t(data, ed(f′, j))
                if fuel(data, ed(f, i)) + fuel(data, ed(i, f′)) <= data.β && 
                  fuel(data, ed(f′, j)) + fuel(data, ed(j, r)) <= data.β && 
                  gvrp_afs_tree.times[f] + t(data, ed(f, i)) + time + t(data, ed(j, r)) + gvrp_afs_tree.times[r] <= data.T
                  minTime = time < minTime ? time : minTime
                else
                  for r′ in F₀
                    if !(ed(r′, j) in data.G′.E)
                      continue
                    end
                    time = t(data, ed(i, f′)) + gvrp_afs_tree.pairTimes[(f′, r′)] + t(data, ed(r′, j))
                    if fuel(data, ed(f, i)) + fuel(data, ed(i, f′)) <= data.β && 
                      fuel(data, ed(r′, j)) + fuel(data, ed(j, r)) <= data.β 
                      if gvrp_afs_tree.times[f] + t(data, ed(f, i)) + time + t(data, ed(j, r)) + gvrp_afs_tree.times[r] <= data.T
                        minTime = time < minTime ? time : minTime
                      else
                        time = t(data, ed(i, f′)) + gvrp_afs_tree.times[f′] + gvrp_afs_tree.times[r′] + t(data, ed(r′, j))
                        minTime = time < minTime ? time : minTime
                      end
                    end
                  end
                end
              end
            end
          end
        end
        times[ed(i, j)] = minTime
      end
    end
  end
  return times
end
