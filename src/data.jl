import Unicode

mutable struct Vertex
    id_vertex::Int
    pos_x::Float64
    pos_y::Float64
end

# Undirected graph
mutable struct InputGraph
    V′::Array{Vertex} # set of vertices
    E::Array{Tuple{Int64,Int64}} # set of edges
    cost::Dict{Tuple{Int64,Int64},Float64} # cost for each edge
end

mutable struct DataBWTSP
    G′::InputGraph
    Black::Array{Int64} # block nodes
    White::Array{Int64} # white nodes
    Q::Int64 # Maximum number of white vertices between two consecutive black vertices
    D::Float64 # Total distance between two consecutive black vertices
    q::Int64
    service_time::Float64 # service time
    insType::String # TSP instance "EDGE_WEIGHT_TYPE"
    coord::Bool # instance with NODE_COORD_SECTION
end

vertices(data::DataBWTSP) = [i.id_vertex for i in data.G′.V′[1:end]] # return set of vertices

function distance(data::DataBWTSP, arc::Tuple{Int64,Int64})
    e = (arc[1] < arc[2]) ? arc : (arc[2], arc[1])
    if haskey(data.G′.cost, e) # use already calculated value
        return data.G′.cost[e]
    elseif data.coord
        u, v = arc
        vertices = data.G′.V′
      # array <vertices> is indexed from 1 (depot is vertices[1], customer 1 is vertices[2], and so on)
        if data.insType == "GEO"
            return GEO_dist(vertices[u], vertices[v])
        elseif data.insType == "EUC_2D"
            return EUC_dist(vertices[u], vertices[v])
        else
            error("$(data.insType) type not found!")
        end
    else
        return 0.0
    end
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

# EUC_2D distance
function EUC_dist(u::Vertex, v::Vertex)
    x_sq = (v.pos_x - u.pos_x)^2
    y_sq = (v.pos_y - u.pos_y)^2
    return floor(sqrt(x_sq + y_sq) + 0.5)
end

contains(p, s) = findnext(s, p, 1) != nothing

function readBWTSPData(app::Dict{String,Any})

    str = Unicode.normalize(read(app["instance"], String); stripcc=true)
    breaks_in = [' '; ':'; '\n']
    aux = split(str, breaks_in; limit=0, keepempty=false)

    G′ = InputGraph([], [], Dict())
    data = DataBWTSP(G′, [], [], 0, 0.0, 0, 0.0, "none", false)

    cost_aux, insType_aux = [], " "
    dim = 0

    for i in 1:length(aux)
        if contains(aux[i], "DIMENSION")
            dim = parse(Int, aux[i + 1])
        elseif contains(aux[i], "EDGE_WEIGHT_TYPE")
            data.insType = aux[i + 1]
        elseif contains(aux[i], "EDGE_WEIGHT_FORMAT")
            insType_aux = aux[i + 1]
        elseif contains(aux[i], "NODE_COORD_SECTION")
            data.coord = true
            j = i + 1
            while aux[j] != "EOF"
                v = Vertex(0, 0, 0)
                v.id_vertex = parse(Int, aux[j])
                v.pos_x = parse(Float64, aux[j + 1])
                v.pos_y = parse(Float64, aux[j + 2])
                push!(G′.V′, v) # add v in the vertex array
                j += 3
            end
        elseif contains(aux[i], "EDGE_WEIGHT_SECTION")
            data.insType = insType_aux
            for k in 1:dim
                v = Vertex(0, 0, 0)
                v.id_vertex = k
                push!(G′.V′, v)
            end
            j = i + 1
            while aux[j] != "EOF"
                push!(cost_aux, parse(Float64, aux[j]))
                j += 1
                if j > length(aux)
                    break
                end
            end
        end
    end

    if data.coord
        for i in vertices(data)
            for j in vertices(data) # add edges between vertices
                if i < j
                    e = (i, j)
                    push!(G′.E, e) # add edge e
                    data.G′.cost[e] = distance(data, e)
                end
            end
        end
    elseif data.insType == "LOWER_DIAG_ROW"
        k = 1
        for i in 1:dim
            for j in 1:dim # add edges between vertices
                if j < i
                    e = (j, i)
                    push!(G′.E, e) # add edge e
                    data.G′.cost[e] = cost_aux[k]
                    k += 1
                end
            end
            k += 1
        end
    elseif data.insType == "FULL_MATRIX"
        k = 1
        for i in 1:dim
            for j in 1:dim # add edges between vertices
                if j > i
                    e = (i, j)
                    push!(G′.E, e) # add edge e
                    data.G′.cost[e] = cost_aux[k]
                end
                k += 1
            end
        end
    end

    nb_black = app["nblack"]
    data.Q = app["maxwhite"]
    data.D = app["maxdist"]
    data.Black = [i.id_vertex for i in G′.V′[1:nb_black]]
    data.White = [i.id_vertex for i in G′.V′[(nb_black+1):end]]
    data.q = app["nb_days"]
    data.service_time = app["service"]

    return data
end

function readBWTSPData2(app::Dict{String,Any})

   str = Unicode.normalize(read(app["instance"], String); stripcc=true)
   breaks_in = [' '; ':'; '\n']
   aux = split(str, breaks_in; limit=0, keepempty=false)

   G′ = InputGraph([],[],Dict())
   data = DataBWTSP(G′, [], [], 0, 0.0, 0, 0.0, "EUC_2D", true)

   nodes = []
   h, c, l = 0, 0, 0

   for i in 1:length(aux)
      h = parse(Int, aux[i])
      c = parse(Int, aux[i+1])
      l = parse(Int, aux[i+2])

      j = i + 3
      last = 1

      while last <= h # hotels vertices
         v = Vertex(0, 0, 0)
         v.id_vertex = last
         v.pos_x = parse(Float64, aux[j+1])
         v.pos_y = parse(Float64, aux[j+2])
         #v.s_time = 0.0
         push!(nodes, v) # add v in the vertex array
         last+=1
         j+=3
      end

      while last <= (h + c) # customers vertices
         v = Vertex(0, 0, 0)
         v.id_vertex = last
         v.pos_x = parse(Float64, aux[j+1])
         v.pos_y = parse(Float64, aux[j+2])
         #v.s_time = parse(Float64, aux[j+3])
         push!(nodes, v) # add v in the vertex array
         last+=1
         j+=4
      end

      if j > length(aux)
         break
      end
   end

   data.Black = [i.id_vertex for i in nodes[1:h]]
   data.White = [i.id_vertex for i in nodes[h + 1:end]]
   data.D = l
   G′.V′ = nodes # add vertices to graph G′
   # println(G′.V′)

   function dist(arc)
       u = data.G′.V′[arc[1]]
       v = data.G′.V′[arc[2]]
       x_sq = (v.pos_x - u.pos_x)^2
       y_sq = (v.pos_y - u.pos_y)^2
       return floor(sqrt(x_sq + y_sq), digits=1)
   end


   for i in vertices(data)
      for j in vertices(data) # add edges between customers and hotels
         if i < j
            e = (i,j)
            push!(G′.E, e) # add edge e
            data.G′.cost[e] = dist(e) # cost edge e
         end
      end
   end
   # println(data.G′.V′)
   # println(data.G′.cost)
   data.q = app["nb_days"]
   data.service_time = app["service"]

   return data
end



edges(data::DataBWTSP) = data.G′.E # return set of edges
c(data,e) = data.G′.cost[e] # cost of the edge e
d(data,e) = (e[1] != e[2]) ? data.G′.cost[e] : 0.0 # cost of the edge e
dimension(data::DataBWTSP) = length(data.G′.V′) # return number of vertices
nb_vertices(data::DataBWTSP) = length(vertices(data))

# return incident edges of i
function δ(data::DataBWTSP, i::Integer)
    incident_edges = Vector{Tuple}()
    for j in 1:i - 1 push!(incident_edges, (j, i)) end
    for j in i + 1:(length(data.G′.V′)) push!(incident_edges, (i, j)) end
    return incident_edges
end
