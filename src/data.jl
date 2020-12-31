import Unicode
using Distances

mutable struct Vertex
    id_vertex::Int
    pos_x::Float64
    pos_y::Float64
    service_time::Float64
end

# Directed graph
mutable struct InputGraph
    V′::Array{Vertex} # set of vertices
    E::Array{Tuple{Int64,Int64}} # set of edges
    cost::Dict{Tuple{Int64,Int64},Float64} # cost for each edge
end

mutable struct DataGVRP
    G′::InputGraph
    F::Array{Int64} # AFSs nodes
    C::Array{Int64} # Customers nodes
    β::Float64 # Total distance between two consecutive black vertices
    T::Float64 # Route time limit
    ρ::Float64 # Vehicle fuel comsumption rate
    ε::Float64 # Vehicle average speed
    m::Float64 # Qtd of vehicles
end

vertices(data::DataGVRP) = [i.id_vertex for i in data.G′.V′[1:end]] # return set of vertices

function distance(data::DataGVRP, e::Tuple{Int64,Int64})
  if haskey(data.G′.cost, e) # use already calculated value
    return data.G′.cost[e]
  end
  u, v = e 
  vertices = data.G′.V′
  # array <vertices> is indexed from 1 (depot is vertices[1], customer 1 is vertices[2], and so on)
  return GEO_dist(vertices[u], vertices[v])
end

# GEO distance
function GEO_dist(u::Vertex, v::Vertex)
    #    return haversine((v.pos_x, v.pos_y), (u.pos_x, u.pos_y), 4182.44949)
    x_sq = (v.pos_x - u.pos_x)^2
    y_sq = (v.pos_y - u.pos_y)^2
    return floor(sqrt(x_sq + y_sq) + 0.5)
end

contains(p, s) = findnext(s, p, 1) != nothing

function readEMHInstance(app::Dict{String,Any})
    G′ = InputGraph([], [], Dict())
    data = DataGVRP(G′, [], [], 0.0, 0.0, 0.0, 0.0,0.0)

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
        v = Vertex(i, 0.0, 0.0, 0.0)
        v.pos_x = parse(Float64, aux[3])
        v.pos_y = parse(Float64, aux[4])
        push!(G′.V′, v) 
        if aux[2] == "f" || aux[2] == "d"
          # Get AFS
          v.service_time = aux[2] == "f" ? 0.25 : 0
          push!(data.F, v.id_vertex)
        elseif aux[2] == "c"
          # Get customer
          v.service_time = 0.5
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
#      data.T = parse(Float64, split(line, ['/']; limit=0, keepempty=false)[2])
      data.T = 10.75 
      # Get vehicle average speed
      line = readline(f)
      data.ε = parse(Float64, split(line, ['/']; limit=0, keepempty=false)[2])
      
      # Get amount of vehicle
      #line = readline(f)
      #data.m = parse(Float64, split(line, ['/']; limit=0, keepempty=false)[2])
    end

    for i in vertices(data)
      for j in vertices(data) # add arcs between vertices
        if i < j
          e = (i, j)
          push!(G′.E, e) # add edge e
          data.G′.cost[e] = distance(data, e)
        end
      end
    end

    return data
end

edges(data::DataGVRP) = data.G′.E # return set of arcs
d(data,e) = (e[1] != e[2]) ? data.G′.cost[e] : 0.0 # cost of the edge e
f(data, e) = d(data, e) * data.ρ # fuel of the arc a 
t(data, e) = d(data, e) / data.ε # time of the arc a 
dimension(data::DataGVRP) = length(data.G′.V′) # return number of vertices
nb_vertices(data::DataGVRP) = length(vertices(data))

# return incident edges of i
function δ(data::DataGVRP, i::Integer)
    incident_edges = Vector{Tuple}()
    for j in 1:i - 1 push!(incident_edges, (j, i)) end
    for j in i + 1:(length(data.G′.V′)) push!(incident_edges, (i, j)) end
    return incident_edges
end
