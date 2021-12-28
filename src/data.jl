import Unicode
using Distances

mutable struct Vertex
    id_vertex::Int
    pos_x::Float64
    pos_y::Float64
    service_time::Float64
    weight::Float64
end

# undirected graph
mutable struct InputGraph
    V′::Array{Vertex} # set of vertices
    E::Array{Tuple{Int64,Int64}} # set of edges
    cost::Dict{Tuple{Int64,Int64},Float64} # cost for each edge
end

mutable struct DataGVRP
    G′::InputGraph
    depot_id ::Int
    coord::Bool # instance with NODE_COORD_SECTION
    round::Bool # Is the distance matrix rounded?
    non_consec::Bool # if the instance must be solved without using edges between AFSs
    F::Array{Int64} # AFSs nodes
    C::Array{Int64} # Customers nodes
    β::Float64 # Total distance between two consecutive black vertices
    T::Float64 # Route time limit
    ρ::Float64 # Vehicle fuel comsumption rate
    ε::Float64 # Vehicle average speed
    reduced_graph::Dict{Tuple{Int64, Int64}, Float64} # reduced graph
    gvrp_afs_tree::Any
end

vertices(data::DataGVRP) = [i.id_vertex for i in data.G′.V′[1:end]] # return set of vertices

function distance(data::DataGVRP, e::Tuple{Int64,Int64})
  if haskey(data.G′.cost, e) # use already calculated value
    return data.G′.cost[e]
  end
  u, v = e 
  vertices = data.G′.V′
  # array <vertices> is indexed from 1 (depot is vertices[1], customer 1 is vertices[2], and so on)
  return haversine((vertices[v].pos_x, vertices[v].pos_y), (vertices[u].pos_x, vertices[u].pos_y), 4182.44949)
end

# EUC_2D distance
function EUC_dist(u::Vertex, v::Vertex)
    x_sq = (v.pos_x - u.pos_x)^2
    y_sq = (v.pos_y - u.pos_y)^2
    return floor(sqrt(x_sq + y_sq) + 0.5)
end

contains(p, s) = findnext(s, p, 1) != nothing

function read_Andelmin_Bartolini_Instance(app::Dict{String,Any})
    data = DataGVRP(InputGraph([], [], Dict{Tuple{Int64, Int64},Float64}()), 1, true, false, false, [], [], 0.0, 0.0, 0.0, 0.0, Dict{Tuple{Int64, Int64}, Float64}(), nothing)
    haskey(app, "non-consec") && app["non-consec"] && error("The instances of Andelmin and Bartolini requires edges between AFSs")
    sepChar = ' '
    FVertices = []
    CVertices = []
    open(app["instance"]) do f
      # get vehicle data
      line = readline(f)
      aux = split(line, [' ']; limit=0, keepempty=false)
      data.T = parse(Float64, aux[4])
      data.β = parse(Float64, aux[5])
      data.ε = parse(Float64, aux[6])
      data.ρ = data.ε > 1.7 ? 0.2137 : 0.2 
      CustomerServiceTime = parse(Float64, aux[7])
      FServiceTime = parse(Float64, aux[8])
      while !eof(f)
        # get vertex
        line = readline(f)
        aux = split(line, [' ']; limit=0, keepempty=false)
        length(aux) == 0 && continue
        if strip(line) == "Infeasible customers"
          !eof(f) && [filter!(v->v.id_vertex ≠ parse(Int64, id), CVertices) for id in split(readline(f), [' ']; limit=0, keepempty=false)]
        elseif aux[2] == "d"
          push!(data.G′.V′, Vertex(parse(Int64, aux[1]), parse(Float64, aux[3]), parse(Float64, aux[4]), 0.0, 0.0))
        elseif aux[2] == "f"
          push!(FVertices, Vertex(parse(Int64, aux[1]), parse(Float64, aux[3]), parse(Float64, aux[4]), FServiceTime, 0.0))
        else
          push!(CVertices, Vertex(parse(Int64, aux[1]), parse(Float64, aux[3]), parse(Float64, aux[4]), CustomerServiceTime, 0.0))
        end
      end
    end
    data.β = data.β * data.ρ
    data.G′.V′[1].id_vertex = 1
    i = 2
    for f in FVertices
      f.id_vertex = i
      push!(data.G′.V′, f)
      push!(data.F, i)
      i = i + 1
    end
    for c in CVertices
      c.id_vertex = i
      push!(data.G′.V′, c)
      push!(data.C, i)
      i = i + 1
    end
    # create edges
    for i in vertices(data)
      for j in vertices(data) # add arcs between vertices
        if i < j 
          e = (i, j)
          push!(data.G′.E, e) # add edge e
          data.G′.cost[e] = round(distance(data, e), digits=6)
        end
      end
    end
#    data.gvrp_afs_tree = calculateGVRP_AFS_Tree(data)
#    data.reduced_graph = calculateGVRPReducedGraphTime(data)

#   invalidEdges = vcat(get_invalid_edges_1(data), get_invalid_edges_2(data), get_invalid_edges_3(data), get_invalid_edges_4(data))
#   data.G′.E = setdiff(data.G′.E, invalidEdges)
#   for e in invalidEdges
#     if e in keys(data.G′.cost)
#       delete!(data.G′.cost, e)
#     end
#   end
    return data
end

function readEMHInstance(app::Dict{String,Any})
    G′ = InputGraph([], [], Dict())
    data = DataGVRP(G′, 1, true, false, true, [], [], 0.0, 0.0, 0.0, 0.0, Dict{Tuple{Int64, Int64}, Float64}(), nothing)
    open(app["instance"]) do f
      # Ignore header
      readline(f)
      # Read vertices
      i = 1
      while !eof(f)
        line = readline(f)
        aux = split(line, ['\t']; limit=0, keepempty=false)
        if length(aux) == 0
          break
        end
        v = Vertex(i, 0.0, 0.0, 0.0, 0.0)
        v.pos_x = parse(Float64, aux[3])
        v.pos_y = parse(Float64, aux[4])
        push!(G′.V′, v) 
        if aux[2] == "d"
          v.service_time = 0
        elseif aux[2] == "f"
          # Get AFS
          v.service_time = aux[2] == "f" ? 15 : 0
          push!(data.F, v.id_vertex)
        elseif aux[2] == "c"
          # Get customer
          v.service_time = 30
          push!(data.C, v.id_vertex)
        end
        i += 1
      end
      # Read params
      # Get beta
      line = readline(f)
      data.β = parse(Float64, split(line, ['/']; limit=0, keepempty=false)[2])
      # Get vehicle fuel consumptin rate
      line = readline(f)
      data.ρ = parse(Float64, split(line, ['/']; limit=0, keepempty=false)[2])
      # Get vehicle time limit
      line = readline(f)
#      data.T = 10.75 
      data.T = 645 
      # Get vehicle average speed
      line = readline(f)
      data.ε = parse(Float64, split(line, ['/']; limit=0, keepempty=false)[2])/60.0
    end

    #read preprocessings
    invalidEdges = []
    if haskey(app, "preprocessings") && app["preprocessings"] != nothing
      open(app["preprocessings"]) do f
        while !eof(f)
          # read edge
          line = readline(f)
          edge = split(line, [' ', ',']; limit=0, keepempty=false)
          push!(invalidEdges, (parse(Int, edge[1]) + 1, parse(Int, edge[2]) + 1))
        end
      end
    end

    # create edges
    for i in vertices(data)
      for j in vertices(data) # add arcs between vertices
        if i < j
          e = (i, j)
          push!(G′.E, e) # add edge e
          data.G′.cost[e] = distance(data, e)
        end
      end
    end
    #remove infeasible customers
    data.gvrp_afs_tree = calculateGVRP_AFS_Tree(data)
    fuel = f
    infeasibleCustomers = Set{Int64}()
    for i in data.C
      feasible = false
      for f in data.gvrp_afs_tree.F0
        for r in data.gvrp_afs_tree.F0
          if f <= r &&
            fuel(data, ed(f, i)) + fuel(data, ed(i, r)) <= data.β && 
            data.gvrp_afs_tree.times[f] + t(data, ed(f, i)) + t(data, ed(i, r)) + data.gvrp_afs_tree.times[r] <= data.T
            feasible = true
            break 
          end
        end
        feasible && break
      end
      !feasible && push!(infeasibleCustomers, i)
    end
    #create new graph
    if !isempty(infeasibleCustomers)
      G′ = InputGraph([data.G′.V′[1]], [], Dict{Tuple{Int64,Int64},Float64}())
      C = []
      F = []
      #insert customers
      id = 2
      for i in data.C
        if !in(i, infeasibleCustomers)
          data.G′.V′[i].id_vertex = id
          push!(C, id)
          push!(G′.V′, data.G′.V′[i])
          id = id + 1
        end
      end
      #insert AFSs
      for i in data.F
        data.G′.V′[i].id_vertex = id
        push!(F, id)
        push!(G′.V′, data.G′.V′[i])
        id = id + 1
      end
      data.G′ = G′
      data.C = C
      data.F = F
      # create edges
      for i in vertices(data)
        for j in vertices(data) # add arcs between vertices
          if i < j
            e = (i, j)
            push!(G′.E, e) # add edge e
            data.G′.cost[e] = distance(data, e)
          end
        end
      end
    end
    data.reduced_graph = calculateGVRPReducedGraphTime(data)
    return data
end

function readMatheusInstance(app::Dict{String,Any})
    G′ = InputGraph([], [], Dict())
    data = DataGVRP(G′, 1, true, false, false, [], [], 0.0, 0.0, 0.0, 0.0, Dict{Tuple{Int64, Int64}, Float64}(), nothing)
    data.non_consec = haskey(app, "non-consec") && app["non-consec"]
    sepChar = ';'
    open(app["instance"]) do f
      # vehicle data
      # ignore header
      readline(f)
      # get vehicle average speed
      line = readline(f)
      data.ε = parse(Float64, split(line, [sepChar]; limit=0, keepempty=false)[2])
      # get vehicle time limit
      line = readline(f)
      data.T = parse(Float64, split(line, [sepChar]; limit=0, keepempty=false)[2])
      # get vehicle fuel consumptin rate
      line = readline(f)
      data.ρ = parse(Float64, split(line, [sepChar]; limit=0, keepempty=false)[2])
      # get beta
      line = readline(f)
      data.β = parse(Float64, split(line, [sepChar]; limit=0, keepempty=false)[2])
      # depot
      # ignore headers
      readline(f)
      readline(f)
      line = readline(f)
      aux = split(line, [sepChar]; limit=0, keepempty=false)
      v = Vertex(parse(Int, aux[1]) + 1, parse(Float64, aux[2]), parse(Float64, aux[3]), parse(Float64, aux[4]), 0.0)
      push!(G′.V′, v) 
      i = 2
      # get customers
      # ignore headers
      line = readline(f)
      line = readline(f)
      while true 
        line = readline(f)
        aux = split(line, [sepChar]; limit=0, keepempty=false)
        if length(aux) == 1
          break
        end
        v = Vertex(parse(Int, aux[1]) + 1, parse(Float64, aux[2]), parse(Float64, aux[3]), parse(Float64, aux[4]), 0.0)
        push!(data.C, i)
        push!(G′.V′, v) 
        i = i + 1
      end
      # get AFSs
      # ignore headers
      line = readline(f)
      while !eof(f) 
        line = readline(f)
        aux = split(line, [sepChar]; limit=0, keepempty=false)
        if length(aux) == 1
          break
        end
        v = Vertex(parse(Int, aux[1]) + 1, parse(Float64, aux[2]), parse(Float64, aux[3]), parse(Float64, aux[4]), 0.0)
        push!(data.F, i)
        push!(G′.V′, v) 
        i = i + 1
      end
    end

    #read preprocessings
    invalidEdges = []
    if haskey(app, "preprocessings") && app["preprocessings"] != nothing
      open(app["preprocessings"]) do f
        while !eof(f)
          # read edge
          line = readline(f)
          edge = split(line, [' ', ',']; limit=0, keepempty=false)
          push!(invalidEdges, (parse(Int, edge[1]) + 1, parse(Int, edge[2]) + 1))
        end
      end
    end

    for i in vertices(data)
      for j in vertices(data) # add arcs between vertices
        e = (i, j)
        a, b = data.G′.V′[i], data.G′.V′[j]
        if i < j && !((a.id_vertex, b.id_vertex) in invalidEdges)
          push!(G′.E, e) # add edge e
          data.G′.cost[e] = EUC_dist(a, b)
        end
      end
    end
    data.gvrp_afs_tree = calculateGVRP_AFS_Tree(data)
    data.reduced_graph = calculateGVRPReducedGraphTime(data)
    return data
end

edges(data::DataGVRP) = data.G′.E # return set of arcs
d(data,e) = (e[1] != e[2]) ? (e in data.G′.E ? data.G′.cost[e] : typemax(Float64)) : 0.0 # cost of the edge e
f(data, e) = d(data, e) * data.ρ # fuel of the arc a 
t(data, e) = (d(data, e) / data.ε) + (data.G′.V′[e[1]].service_time + data.G′.V′[e[2]].service_time)/2.0 # time of the arc a 
dimension(data::DataGVRP) = length(data.G′.V′) # return number of vertices
nb_vertices(data::DataGVRP) = length(vertices(data))

# return incident edges of i
function δ(data::DataGVRP, i::Integer)
#    return vcat([(j, i) for j in 1:i - 1 if (j, i) in data.G′.E], [(i, j) for j in i + 1:(length(data.G′.V′)) if (i, j) in data.G′.E])
  return [(j, k) for (j, k) in data.G′.E if j == i || k == i]
end

# return endering arcs of i
function δ⁻(data::DataGVRP, i::Integer)
  return [(j, k) for (j, k) in data.G′.E if k == i]
end

# return endering arcs of i
function δ⁺(data::DataGVRP, i::Integer)
  return [(j, k) for (j, k) in data.G′.E if j == i]
end

function removeInfeasibleCustomers(data::DataGVRP)
end
