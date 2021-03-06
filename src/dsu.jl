mutable struct DSU
  pred::Array{Int64}
  rank::Array{Int64}
  n::Int64
end

function createDSU(n::Int64)
  return DSU([i for i in 1:n], [0 for i in 1:n], n)
end

function joinDSU(dsu::DSU, i::Int64, j::Int64)
  if i != j && i >= 1 && j >= 1 && i <= dsu.n && j <= dsu.n
    a, b = findSetDSU(dsu, i), findSetDSU(dsu, j)
    if dsu.rank[a] < dsu.rank[b]
      a, b = b, a
    end
    dsu.pred[b] = a
    if dsu.rank[a] == dsu.rank[b]
      dsu.rank[a] = dsu.rank[a] + 1
    end
  end
end

function findSetDSU(dsu::DSU, i::Int64)
  if dsu.pred[i] == i
    return i
  end
  dsu.pred[i] = dsu.pred[findSetDSU(dsu, dsu.pred[i])]
  return dsu.pred[i]
end

function cleanDSU(dsu::DSU)
  for i in 1:dsu.n
    dsu.pred[i] = i
    dsu.rank[i] = 0
  end
end
