mutable struct GVRP_AFS_Tree
  F0::Array{Int64} 
  times::Dict{Int64,Float64}
  pred::Dict{Int64,Int64}
  pairCosts::Dict{Tuple{Int64,Int64},Float64}
  pairTimes::Dict{Tuple{Int64,Int64},Float64}
  pairPreds::Dict{Tuple{Int64,Int64},Int64}
end

function calculateGVRP_AFS_Tree(data::DataGVRP)
  F0 = vcat([data.depot_id], data.F)
  pairTimes = Dict{Tuple{Int64,Int64},Float64}((f, r) => typemax(Float64) for f in F0 for r in F0);
  pairCosts = Dict{Tuple{Int64,Int64},Float64}((f, r) => typemax(Float64) for f in F0 for r in F0);
  pairPreds = Dict{Tuple{Int64,Int64},Int64}()
  #dijkstra
  fuel = f
  for f in F0
    for r in F0
      pairPreds[(f, r)] = r
    end
    pairTimes[(f, f)] = 0.0;
    pairCosts[(f, f)] = 0.0;
    q = [f]
    while !isempty(q)
      curr = popfirst!(q);
      for r in F0
        if ed(curr, r) in data.G′.E
          cost = d(data, ed(curr, r))
          time = t(data, ed(curr, r))
          if fuel(data, ed(curr, r)) <= data.β && 
            pairCosts[(f, curr)] + cost < pairCosts[(f, r)] && 
            pairTimes[(f, curr)] + time < pairTimes[(f, r)] 
            pairCosts[(f, r)] = pairCosts[(f, curr)] + cost
            pairTimes[(f, r)] = pairTimes[(f, curr)] + time
            pairPreds[(f, r)] = curr
            push!(q, r)
          end
        end
      end
    end
  end
  times = Dict{Int64,Float64}(i => pairTimes[(data.depot_id, i)] for i in F0)
  pred = Dict{Int64,Int64}(i => pairPreds[(data.depot_id, i)] for i in F0)
  return GVRP_AFS_Tree(F0, times, pred, pairCosts, pairTimes, pairPreds)
end
