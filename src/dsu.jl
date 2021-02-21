mutable struct DSU
  pred::Array{Int64}
  rank::Array{Int64}
  n::Int64
end

function createDSU(n::Int64)
  return DSU([i for i in 1:n], [0 for i in 1:n], n)
end

function join(dsu::DSU, i::Int64, j::Int64)
  if i != j && i < n && j < n
    a, b = findSet(dsu, i), findSet(dsu, j)
    if dsu.rank[a] < dsu.rank[b]
      a, b = b, a
    end
    dsu.pred[b] = a
    if dsu.rank[a] == dsu.rank[b]
      dsu.rank[a] = dsu.rank[a] + 1
    end
  end
end

function findSet(dsu::DSU, i::Int64)
  if (dsu.pred[i] == i)
    return i
  end
  dsu.pred[i] = dsu.pred[dsu.findSet(dsu.pred[i])]
  return dsu.pred[i]
end

function cleanDSU(duse::DSU)
  for i in 1:n
    dsu.pred[i] = i
    dsu.rank[i] = 0
  end
end
