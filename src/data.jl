import Unicode

mutable struct Vertex
    id_vertex::Int
    pos_x::Float64
    pos_y::Float64
    service_time::Float64
end

# Directed graph
mutable struct InputGraph
    V′::Array{Vertex} # set of vertices
    A::Array{Tuple{Int64,Int64}} # set of edges
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
end

vertices(data::DataGVRP) = [i.id_vertex for i in data.G′.V′[1:end]] # return set of vertices

function distance(data::DataGVRP, a::Tuple{Int64,Int64})
  if haskey(data.G′.cost, a) # use already calculated value
    return data.G′.cost[a]
  end
  u, v = a
  vertices = data.G′.V′
  # array <vertices> is indexed from 1 (depot is vertices[1], customer 1 is vertices[2], and so on)
  return GEO_dist(vertices[u], vertices[v])
end

# GEO distance
function GEO_dist(u::Vertex, v::Vertex)
    PI = 3.141592
    RRR = 6378.388

    function coordenate(pos)
        deg = trunc(pos) # round(pos)
        min = pos - deg
        return PI * (deg + 5.0 * min / 3.0 ) / 180.0;
    end

    lat_i = coordenate(v.pos_x)
    long_i = coordenate(v.pos_y)
    lat_j = coordenate(u.pos_x)
    long_j = coordenate(u.pos_y)

    q1 = cos(long_i - long_j)
    q2 = cos(lat_i - lat_j)
    q3 = cos(lat_i + lat_j)

    return floor(RRR * acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0)
end

contains(p, s) = findnext(s, p, 1) != nothing

function readEMHIntance(app::Dict{String,Any})
    G′ = InputGraph([], [], Dict())
    data = DataGVRP(G′, [], [], 0.0, 0.0, 0.0, 0.0)

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
        if aux[2] == "f"
          # Get AFS
          push!(data.F, v.id_vertex)
        elseif aux[2] == "c"
          # Get customer
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
      data.T = parse(Float64, split(line, ['/']; limit=0, keepempty=false)[2])
      # Get vehicle average speed
      line = readline(f)
      data.ε = parse(Float64, split(line, ['/']; limit=0, keepempty=false)[2])
    end

    for i in vertices(data)
      for j in vertices(data) # add arcs between vertices
        a = (i, j)
        push!(G′.A, a) # add edge e
        data.G′.cost[a] = distance(data, a)
      end
    end

    return data
end

arcs(data::DataGVRP) = data.G′.A # return set of arcs
d(data,a) = data.G′.cost[a] # distance of the arc a 
f(data,a) = d(data, a) * data.ρ # fuel of the arc a 
t(data,a) = d(data, a) / data.ε # time of the arc a 
dimension(data::DataGVRP) = length(data.G′.V′) # return number of vertices
nb_vertices(data::DataGVRP) = length(vertices(data))

# return incident edges of i
function δ(data::DataGVRP, i::Integer)
    incident_edges = Vector{Tuple}()
    for j in 1:i - 1 push!(incident_edges, (j, i)) end
    for j in i + 1:(length(data.G′.V′)) push!(incident_edges, (i, j)) end
    return incident_edges
end
