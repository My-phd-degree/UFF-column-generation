function calculateGVRPReducedGraphTime(data::DataGVRP) 
  depot = data.depot_id
  times, C₀, F₀, T, β = Dict{Tuple{Int64,Int64}, Float64}(), vcat([depot], data.C), vcat([depot], data.F), data.T, data.β
  fuel(i, j) = f(data, ed(i, j))
  t′(i, j) = t(data, ed(i, j))
#  nondominatedAFSs = Dict{Int64, Array{Int64}}(i => [f for f in F₀ if all([fuel(i, f) <= fuel(i, r) || data.gvrp_afs_tree.times[f] + t′(i, f) <= data.gvrp_afs_tree.times[r] + t′(i, r) for r in F₀]) && fuel(i, f) <= β && t′(depot, i) + t′(i, f) + data.gvrp_afs_tree.times[f] <= T] for i in C₀)
  for i in C₀
    for j in C₀
      i >= j && continue
#      println("Customers pair $i and $j")
      bestFound, minTime, time = false, typemax(Float64), t′(i, j)
#      for f in nondominatedAFSs[i]
      for f in F₀
        bestFound && break
        fuel′, time′ = fuel(f, i), data.gvrp_afs_tree.times[f] + t′(f, i)
#        for r in nondominatedAFSs[j]
        for r in F₀
#          println("\tExternal AFSs pair $f and $r")
          (fuel′ + fuel(i, j) + fuel(j, r) <= β &&  time′ + time + t′(j, r) + data.gvrp_afs_tree.times[r] <= T) && (minTime = time; bestFound = true; break)
#          for f′ in nondominatedAFSs[i]
          for f′ in F₀
            time″ = time′ + t′(i, f′)
            fuel′ + fuel(i, f′) > β && continue
            (f′ == depot && (time″ > T || t′(f′, j) + t′(j, r) + data.gvrp_afs_tree.times[r] > T)) && continue
            (f′ != depot && time″ + t′(f′, j) + t′(j, r) + data.gvrp_afs_tree.times[r] > T) && continue
            (fuel(f′, j) + fuel(j, r) <= β) && (minTime = min(t′(i, f′) + t′(f′, j), minTime); continue)
            f′ == depot && continue
#            for r′ in nondominatedAFSs[j]
            for r′ in F₀
              fuel(r′, j) + fuel(j, r) > β && continue
              time‴ = time″ + data.gvrp_afs_tree.pairTimes[(f′, r′)] + t′(r′, j)
              time‴ + t′(j, r) + data.gvrp_afs_tree.times[r] <= T && (minTime = min(t′(i, f′) + data.gvrp_afs_tree.pairTimes[(f′, r′)] + t′(r′, j),  minTime))
              (time″ + data.gvrp_afs_tree.times[f′] <= T && data.gvrp_afs_tree.times[r′] + t′(r′, j) + t′(j, r) + data.gvrp_afs_tree.times[r] <= T) && (minTime = min(t′(i, f′) + data.gvrp_afs_tree.times[f′] + data.gvrp_afs_tree.times[r′] + t′(r′, j), minTime))
            end
          end
        end
      end
      times[ed(i, j)] = minTime
    end
  end
  for p in keys(times)
#  println(p, " ", times[p])
  end
  return times
end
