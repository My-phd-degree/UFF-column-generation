import Unicode
#using Distances

mutable struct Vertex
    id_vertex::Int64
    pos_x::Float64
    pos_y::Float64
    service_time::Float64
end

# Directed graph
mutable struct InputGraph
    V´::Array{Vertex} # set of vertices
    E::Array{Tuple{Int64,Int64}} # set of edges
    cost::Dict{Tuple{Int64,Int64},Float64} # cost for each edge
end

mutable struct DataGVRP
    G´::InputGraph
    F::Array{Int64} # AFSs nodes
    C::Array{Int64} # Customers nodes
    M::Array{Int64} # Vehicles IDs
    E´::Array{Tuple{Int64,Int64}}
    FM::Array{Tuple{Int64,Int64}}
    β::Float64 # Total distance between two consecutive black vertices
    T::Float64 # Route time limit
    ρ::Float64 # Vehicle fuel comsumption rate
    ε::Float64 # Vehicle average speed
    m::Int64 # Qtd of vehicles
    max_ed::Float64
    min_ed::Float64
    LB_E::Array{Float64}
end

vertices(data::DataGVRP) = [i.id_vertex for i in data.G´.V´[1:end]] # return set of vertices

function distance(data::DataGVRP, e::Tuple{Int64,Int64})
  if haskey(data.G´.cost, e) # use already calculated value
    return data.G´.cost[e]
  end
  u, v = e 
  vertices = data.G´.V´
  flag = 1
  # array <vertices> is indexed from 1 (depot is vertices[1], customer 1 is vertices[2], and so on)
  return (flag == 1) ? GEO_dist(vertices[u], vertices[v]) : EUC_dist(vertices[u], vertices[v])  
end

# GEO distance
function GEO_dist(u::Vertex, v::Vertex)
    #if (u.id_vertex == 2 && v.id_vertex == 6)
    #  println("............................................")
    #  println(v.pos_x," ",u.pos_x)
    #  println(v.pos_y," ",u.pos_y)
    #end
    #    return haversine((v.pos_x, v.pos_y), (u.pos_x, u.pos_y), 4182.44949)
    x_sq = (v.pos_x - u.pos_x)^2
    y_sq = (v.pos_y - u.pos_y)^2
    return floor(sqrt(x_sq + y_sq) + 0.5)
    #return sqrt(x_sq + y_sq) + 0.5
end

# EUC_2D distance
function EUC_dist(u::Vertex, v::Vertex)
    x_sq = (v.pos_x - u.pos_x)^2
    y_sq = (v.pos_y - u.pos_y)^2
    return floor(sqrt(x_sq + y_sq) + 0.5)
end

contains(p, s) = findnext(s, p, 1) != nothing

function readEMHInstance(app::Dict{String,Any})
    G´ = InputGraph([], [], Dict())
    data = DataGVRP(G´, [], [], [], [], [], 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0,[])

    open(app["instance"]) do f
      # Ignore header
      readline(f)
      # Read vertices
      i = 1
      data.max_ed = 0.0
      data.min_ed = 999999999.9999
      while !eof(f)
        line = readline(f)
        aux = split(line, ['\t']; limit=0, keepempty=false)
        if length(aux) == 0
          break
        end
        v = Vertex(i, 0.0, 0.0, 0.0)
        v.pos_x = parse(Float64, aux[3])
        v.pos_y = parse(Float64, aux[4])
        push!(G´.V´, v) 
        
        if aux[2] == "f" || aux[2] == "d"
          # Get AFS
          v.service_time = aux[2] == "f" ? 0.25 : 0.0
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
      # data.T = parse(Float64, split(line, ['/']; limit=0, keepempty=false)[2])
      data.T = 10.75 
      # Get vehicle average speed
      line = readline(f)
      data.ε = parse(Float64, split(line, ['/']; limit=0, keepempty=false)[2])
      # Get amount of vehicle
      line = readline(f)
      #data.m = parse(Int64, split(line, ['/']; limit=0, keepempty=false)[2])
      data.m = length(data.C)
      
      for k in 1:data.m 
        push!(data.M, k)
        #push!(data.FM, data.m)
        #for f in 0:length(F)-1
        #  push!(data.FM[k], length(F)-1 )
        #end
      end
    end

    #read preprocessings
    data.E´ = []
    if haskey(app, "preprocessings") && app["preprocessings"] != nothing
      open(app["preprocessings"]) do f
        while !eof(f)
          # read edge
          line = readline(f)
          edge = split(line, [' ', ',']; limit=0, keepempty=false)
          push!(data.E´, (parse(Int64, edge[1]) , parse(Int64, edge[2]) ))
        end
      end
    end

    for i in vertices(data)
      for j in vertices(data) # add arcs between vertices
        e = (i, j)
        if haskey(app, "preprocessings") && app["preprocessings"] != nothing
          push!(G´.E, e) # add edge e
          if i < j && !((i, j) in data.E´)
            data.G´.cost[e] = distance(data, e)
            data.max_ed = ( data.max_ed < f(data, ed(i, j)) ) ? f(data, ed(i, j)) : 0.0
            data.min_ed = ( data.min_ed > f(data, ed(i, j)) ) ? f(data, ed(i, j)) : 999999999.9999
          else
            data.G´.cost[e] = 999999999.9999
          end
        elseif i < j
            push!(G´.E, e) # add edge e
            data.G´.cost[e] = distance(data, e)
            #println( data.max_ed, " <= " , f(data, ed(i, j)) , " ", (data.max_ed <= f(data, ed(i, j)) ) ? true : false )
            data.max_ed = ( data.max_ed < f(data, ed(i, j)) ) ? f(data, ed(i, j)) : data.max_ed
            data.min_ed = ( data.min_ed > f(data, ed(i, j)) ) ? f(data, ed(i, j)) : 999999999.9999
        end
      end
    end
    min_LB_E_j(data)
    return data
end

# fix bugs
function readMatheusInstance(app::Dict{String,Any})
    G´ = InputGraph([], [], Dict())
    data = DataGVRP(G´, [], [], [], [], [], 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, [])
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
      v = Vertex(parse(Int64, aux[1]) + 1, parse(Float64, aux[2]), parse(Float64, aux[3]), parse(Float64, aux[4]))
      i = 1
      push!(data.F, i)
      push!(G´.V´, v) 
      i = i + 1
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
        v = Vertex(parse(Int64, aux[1]) + 1, parse(Float64, aux[2]), parse(Float64, aux[3]), parse(Float64, aux[4]))
        push!(data.C, i)
        push!(G´.V´, v) 
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
        v = Vertex(parse(Int64, aux[1]) + 1, parse(Float64, aux[2]), parse(Float64, aux[3]), parse(Float64, aux[4]))
        push!(data.F, i)
        push!(G´.V´, v) 
        i = i + 1
      end
      data.m = length(data.C)
      for k in 1:data.m push!(data.M, k) end
      #data.M = [k for k in 0:data.m]
    end

    #read preprocessings
    data.E´ = []
    if haskey(app, "preprocessings") && app["preprocessings"] != nothing
      open(app["preprocessings"]) do f
        while !eof(f)
          # read edge
          line = readline(f)
          edge = split(line, [' ', ',']; limit=0, keepempty=false)
          push!(data.E´, (parse(Int64, edge[1]) + 1, parse(Int64, edge[2]) + 1))
        end
      end
    end

    for i in vertices(data)
      for j in vertices(data) # add arcs between vertices
        e = (i, j)
        a, b = data.G´.V´[i], data.G´.V´[j]
        if i < j && !((a.id_vertex, b.id_vertex) in data.E´)
          vertices = data.G´.V´
          push!(G´.E, e) # add edge e
          data.G´.cost[e] = distance(data,e)
          #data.G´.cost[e] = EUC_dist(vertices[e[1]], vertices[e[2]])
        end
      end
    end

    return data
end

#&& !((e[1], e[2]) in data.E´) 
edges(data::DataGVRP) = data.G´.E # return set of arcs
d(data,e) = (e[1] != e[2] && !((e[1], e[2]) in data.E´) ) ? data.G´.cost[e] : 999999999.9999 # cost of the edge e
f(data,e) = (e[1] != e[2] && !((e[1], e[2]) in data.E´) ) ? d(data, e) * data.ρ : 999999999.9999 # fuel of the arc a 
t(data,e) = (e[1] != e[2] && !((e[1], e[2]) in data.E´) ) ? d(data, e) / data.ε : 999999999.9999 # time of the arc a 
dimension(data::DataGVRP) = length(data.G´.V´) # return number of vertices
nb_vertices(data::DataGVRP) = length(vertices(data))

function min_LB_E_j(data::DataGVRP)
  C´ = Array{Int64}
  C´ = []
  push!(C´, 1 )
  for j in data.C
    push!(C´, j )
  end
  e_jf = Array{Int64}
  e_jr = Array{Int64}
  e_jf = []
  e_jr = []
  push!(e_jf,0) 
  push!(e_jf,0)
  push!(e_jr,0) 
  push!(e_jr,0)

  # se vc se refere ao  t_{f}^{'} da tese, eles são calculados pelo 
  # caminho minimo de f ao deposito no grapfo induzido de AFSs com arestas de peso <= beta
  
  # bem se ele tem limite de tempo <= T e todas as arestas tem combustivel <= beta, 
  # então ele é viaável
  for j in C´
    for _f in data.F
      for _r in data.F
        t_f = t_r = 0
        for i in 1:1:nb_vertices(data)
          t_f += t( data, ed( i, _f ) )
          t_r += t( data, ed( i, _r ) )
        end
        if t_f + t( data, ed( _f, j ) ) + t( data, ed( j, _r ) ) + t_r <= data.T && f( data, ed( _f, j ) ) + f( data, ed( j, _r ) ) <= data.β
          e_jf[1], e_jf[2] = j, _f
          e_jr[1], e_jr[2] = j, _r
        end
        #LB_E[i] = f(data, ed(i, j))
      end
    end
    push!( data.LB_E, min( f( data, ed( e_jf[1], e_jf[2] ) ) , f( data, ed( e_jr[1], e_jr[2] ) ) ) )
  end
end

function lowerBoundNbVehicles(data::DataGVRP) 
   sum_demand = 0
   alpha = nb_vertices(data)
   println("lowerBoundNbVehicles")
   for i in data.G´.V´
      println(i)
   end

   for i in vertices(data)
      println(i)
      sum_demand += ceil(alpha*data.G´.V´[i].service_time)
   end
   return Int( ceil( sum_demand / (ceil(data.T) ) ) )
end

# return incident edges of i
function δ(data::DataGVRP, i::Integer)
    incident_edges = Vector{Tuple}()
    for j in 1:i - 1 push!(incident_edges, (j, i)) end
    for j in i + 1:(length(data.G´.V´)) push!(incident_edges, (i, j)) end
    return incident_edges
end