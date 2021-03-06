include("bpp.jl")
include("dsu.jl")

function calculateClosestsCustomers(data::DataGVRP, S₀::Array{Int64}, reduced_graph::Dict{Tuple{Int64, Int64}, Float64})
  η = Dict{Int64, Float64}()
  pi = Dict{Int64, Float64}()
  for i in S₀
    closest, secondClosest = typemax(Float64), typemax(Float64)
    for j in S₀
      if i != j
        if reduced_graph[ed(i, j)] - (data.G′.V′[i].service_time + data.G′.V′[j].service_time)/2.0 < closest
          secondClosest = closest
          closest = reduced_graph[ed(i, j)] - (data.G′.V′[i].service_time + data.G′.V′[j].service_time)/2.0
        elseif reduced_graph[ed(i, j)]  - (data.G′.V′[i].service_time + data.G′.V′[j].service_time)/2.0 < secondClosest
          secondClosest = reduced_graph[ed(i, j)] - (data.G′.V′[i].service_time + data.G′.V′[j].service_time)/2.0
        end
      end
    end
    η[i] = closest
    pi[i] = secondClosest
  end
  return η, pi
end

function calculateGvrpLBByImprovedMST(data::DataGVRP, S₀::Array{Int64}, η::Dict{Int64,Float64}, pi::Dict{Int64,Float64}, gvrpReducedGraph::Dict{Tuple{Int64,Int64}, Float64}) 
  #boruvka
  bestLB = 0.0
  adjMatrix = Set{Tuple{Int64, Int64}}()
  dsu = createDSU(max(S₀...));
  nTrees = length(S₀)
  MSTCost = 0.0
  bestEdge = Dict{Int64, Tuple{Int64, Int64}}()
  bestEdgeCost = Dict{Int64, Float64}(i => typemax(Float64) for i in S₀)
  while nTrees > 1 
    #calculate best cuts
    for i in S₀
      setI = findSetDSU(dsu, i)
      for j in S₀
        e = ed(i, j)
        if setI != findSetDSU(dsu, j) && i != j && e in data.G′.E
          cost = gvrpReducedGraph[e] - (data.G′.V′[i].service_time + data.G′.V′[j].service_time)/2.0
          if cost < bestEdgeCost[setI] 
            bestEdgeCost[setI] = cost
            bestEdge[setI] = e
          end
        end
      end
    end
    #insert best edges
    for i in S₀
      setI = findSetDSU(dsu, i)
      #if edge exist
      if setI in keys(bestEdge)
        (j, k) = bestEdge[setI]
        setJ, setK = findSetDSU(dsu, j), findSetDSU(dsu, k)
        if setJ != setK
          push!(adjMatrix, ed(j, k))
          joinDSU(dsu, j, k)
          nTrees = nTrees - 1
          MSTCost = MSTCost + bestEdgeCost[setI]
        end
      end
    end
    bestEdge = Dict{Int64, Tuple{Int64, Int64}}()
    bestEdgeCost = Dict{Int64, Float64}(i => typemax(Float64) for i in S₀)
  end
  MSTCost_ = 0.0
  for e in adjMatrix
    MSTCost_ = MSTCost_ + gvrpReducedGraph[e] - (data.G′.V′[e[1]].service_time + data.G′.V′[e[2]].service_time)/2.0
#    println("(", data.G′.V′[e[1]].id_vertex, ", ", data.G′.V′[e[2]].id_vertex, "): ", gvrpReducedGraph[e] - (data.G′.V′[e[1]].service_time + data.G′.V′[e[2]].service_time)/2.0)
  end
#  println(MSTCost)
#  println(MSTCost_)
#  println()
  #greedy boruvka
  for i in S₀
    cleanDSU(dsu)
#    println("MST without ", data.G′.V′[i].id_vertex)
#    print("\tRemoved edges: ")
    #remove i from the graph   
    for j in S₀
      e = ed(i, j)
      if e in adjMatrix
#        print("(", data.G′.V′[e[1]].id_vertex, ", ", data.G′.V′[e[2]].id_vertex, "): ", (gvrpReducedGraph[e] - (data.G′.V′[e[1]].service_time + data.G′.V′[e[2]].service_time)/2.0), ", ")
        MSTCost = MSTCost - (gvrpReducedGraph[e] - (data.G′.V′[e[1]].service_time + data.G′.V′[e[2]].service_time)/2.0)
        delete!(adjMatrix, e)
      end
    end
#    println()
#    println("\tnew MST cost ", MSTCost)
#    print("\tnew edges: ")
    #populate dsu
    for (j, k) in adjMatrix
      joinDSU(dsu, j, k)
#      print("(", data.G′.V′[j].id_vertex, ", ", data.G′.V′[k].id_vertex, ")")
    end
#    println()
    sets = Set{Int64}(findSetDSU(dsu, i) for i in S₀)
    #get number of components
    nTrees = length(sets)
    #borukva
    while nTrees > 2
#      println("\t# trees ", nTrees)
#      println("\tselecting edges: ")
      #calculate best cuts
      for j in S₀
        if i != j
          setJ = findSetDSU(dsu, j)
          for k in S₀
            e = ed(j, k)
            if k != i && k != j && e in data.G′.E
              cost = gvrpReducedGraph[e] - (data.G′.V′[j].service_time + data.G′.V′[k].service_time)/2.0
              if setJ != findSetDSU(dsu, k) && cost < bestEdgeCost[setJ] 
#                println("(", data.G′.V′[e[1]].id_vertex, ", ", data.G′.V′[e[2]].id_vertex, "): ", cost, "<=", bestEdgeCost[setJ])
                bestEdgeCost[setJ] = cost
                bestEdge[setJ] = e 
              end
            end
          end
        end
      end
#      print("\tinserted edges: ")
      #insert best edges
      for j in S₀
        setJ = findSetDSU(dsu, j)
        #if edge exist
        if setJ in keys(bestEdge)
          (k, l) = bestEdge[setJ]
          setK = findSetDSU(dsu, k)
          setL = findSetDSU(dsu, l)
          if setK != setL
#            print("(", data.G′.V′[k].id_vertex, ", ", data.G′.V′[l].id_vertex, "): ", bestEdgeCost[setJ], ", ")
            push!(adjMatrix, ed(k, l))
            joinDSU(dsu, setK, setL)
            nTrees = nTrees - 1
            MSTCost = MSTCost + bestEdgeCost[setJ]
          end 
        end
      end
#      println()
#      println("\tcurrent MST cost ", MSTCost)
      bestEdge = Dict{Int64, Tuple{Int64, Int64}}()
      bestEdgeCost = Dict{Int64, Float64}(i => typemax(Float64) for i in S₀)
      bestLB = max(η[i] + pi[i] + MSTCost, bestLB)
    end
  end
  return bestLB + sum([data.G′.V′[i].service_time for i in data.C])
end

function calculateGVRP_BPP_NRoutesLB(data::DataGVRP, S₀::Array{Int64}, η::Dict{Int64,Float64}, pi::Dict{Int64,Float64}) 
  dataBPP = DataBPP([data.G′.V′[i].service_time + (η[i] + pi[i])/2.0 for i in S₀], data.T)
  return solveBPP(dataBPP);
end

function calculateGVRP_NRoutesLB(data::DataGVRP, S₀::Array{Int64}, reduced_graph::Dict{Tuple{Int64, Int64}, Float64})
  η, pi  = calculateClosestsCustomers(data, S₀, reduced_graph)
  return max(calculateGVRP_BPP_NRoutesLB(data, S₀, η, pi), ceil(calculateGvrpLBByImprovedMST(data, S₀, η, pi, reduced_graph)/data.T))
end

