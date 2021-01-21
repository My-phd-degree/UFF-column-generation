import Unicode
#using Distances
using Random
using Crayons
using Crayons.Box

MAX_INT     = 999999999
MAX_DOUBLE  = 999999999

MIN_DOUBLE  = -9.999999999999e8

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
    max_d::Float64
    min_d::Float64
    max_f::Float64
    min_f::Float64
    max_t::Float64
    min_t::Float64
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
  return (flag == 1) ? EUC_dist(vertices[u], vertices[v]) : GEO_dist(vertices[u], vertices[v])   
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
    #return floor(sqrt(x_sq + y_sq) + 0.5)
    #return ceil(sqrt(x_sq + y_sq) + 0.5)
    return sqrt(x_sq + y_sq) + 0.5
end

# EUC_2D distance
function EUC_dist(u::Vertex, v::Vertex)
    x_sq = (v.pos_x - u.pos_x)^2
    y_sq = (v.pos_y - u.pos_y)^2
    return floor(sqrt(x_sq + y_sq) + 0.5)
    #return ceil(sqrt(x_sq + y_sq) + 0.5)
    #return sqrt(x_sq + y_sq) + 0.5
end

contains(p, s) = findnext(s, p, 1) != nothing

function readEMHInstance(app::Dict{String,Any})
    G´ = InputGraph([], [], Dict())
    data = DataGVRP(G´, [], [], [], [], [], 0.0,0.0,0.0,0.0,0,0.0,0.0,0.0,0.0,0.0,0.0,[])
    
    #graph = Vector{Tuple{Int64,Float64,Float64}}
    #graph = [(0.0, 0.0, 0.0)]
    #graph = [(1, 2, 7.0),  (1, 3, 9.0),  (1, 6, 14.0), (2, 3, 10.0),
    #         (2, 4, 15.0), (3, 4, 11.0), (3, 6, 2.0),  (4, 5, 6.0),
    #         (5, 6, 9.0)]

    graph = [(1, 1, 999999999.9),  (2, 2, 999999999.9),  (3, 3, 999999999.9), (4, 4, 999999999.9),
             (5, 5, 999999999.9), (6, 6, 999999999.9), (7, 7, 999999999.9),  (8, 8, 999999999.9),
             (9, 9, 999999999.9)]

    open(app["instance"]) do f
      # Ignore header
      readline(f)
      # Read vertices
      i = 1
      data.max_d = data.max_f = data.max_t = 0.0
      data.min_d = data.min_f = data.min_t = 999999999.9
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
    #graph = []

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
          cost = 999999999.9
          if i < j && !((i, j) in data.E´)
            data.G´.cost[e] = distance(data, e)
            data.max_d  = ( data.max_d < d(data, ed(i, j)) ) ? d(data, ed(i, j)) : data.max_d
            data.max_f  = ( data.max_f < f(data, ed(i, j)) ) ? f(data, ed(i, j)) : data.max_f
            data.max_t  = ( data.max_t < t(data, ed(i, j)) ) ? t(data, ed(i, j)) : data.max_t

            data.min_d = ( data.min_d > d(data, ed(i, j)) ) ? d(data, ed(i, j)) : data.min_d
            data.min_f = ( data.min_f > f(data, ed(i, j)) ) ? f(data, ed(i, j)) : data.min_f
            data.min_t = ( data.min_t > t(data, ed(i, j)) ) ? t(data, ed(i, j)) : data.min_t
            cost = convert(Float64, t(data, ed(i, j)) ) + data.G´.V´[j].service_time
            push!(graph, (i, j, cost))
          else
            data.G´.cost[e] = 999999999.9
            push!(graph, (i, j, cost))
          end
          #append!( graph, [ (i, j , cost ) ] )
        elseif i < j
            push!(G´.E, e) # add edge e
            data.G´.cost[e] = distance(data, e)
            #cost = ceil( t(data, ed(i, j)) )
            #println(typeof( t(data, ed(i, j)) ))
            cost = convert(Float64, t(data, ed(i, j)) ) + data.G´.V´[j].service_time
            #println(typeof(cost))
            
            #cost2 = 9.0#(data, ed(i, j))
            
            push!(graph, (i, j, cost))
            #append!( graph, [ (i, j , cost ) ] )

            data.max_d  = ( data.max_d < d(data, ed(i, j)) ) ? d(data, ed(i, j)) : data.max_d
            data.max_f  = ( data.max_f < f(data, ed(i, j)) ) ? f(data, ed(i, j)) : data.max_f
            data.max_t  = ( data.max_t < t(data, ed(i, j)) ) ? t(data, ed(i, j)) : data.max_t

            data.min_d = ( data.min_d > d(data, ed(i, j)) ) ? d(data, ed(i, j)) : data.min_d
            data.min_f = ( data.min_f > f(data, ed(i, j)) ) ? f(data, ed(i, j)) : data.min_f
            data.min_t = ( data.min_t > t(data, ed(i, j)) ) ? t(data, ed(i, j)) : data.min_t
        end
      end
    end
    
    #Shortest path from a to e: a → c → d → e (cost 26)
    g = Digraph(graph)
    min_LB_E_j(data, g)
    
    #println("BFS ...")
    #println( "cost_BFS_path (1, 2): " , cost_BFS_path(data,1,2) )
    #println( "cost_BFS_path (1, 3): " , cost_BFS_path(data,1,3) )
    #println( "cost_BFS_path (1, 4): " , cost_BFS_path(data,1,4) )
    
    """
    C´ = Array{Int64}
    C´ = []
    push!(C´, 1 )
    for j in data.C
      push!(C´, j )
    end
    for j in C´
    push!( data.LB_E, 0.0 )
    end
    """

    println("DijkstraPath ...")
    src, dst = 4, 9
    path, cost = dijkstrapath(g, src, dst)
    println("Shortest path from $src to $dst: ", isempty(path) ? "no possible path" : join(path, " → "), " (cost $cost)")

    print_matrix(data, "time_cost")
    #print_matrix(data, "fuel_cost")
    #print_matrix(data, "distance_cost")
    return data
end

# fix bugs
function readMatheusInstance(app::Dict{String,Any})
    G´ = InputGraph([], [], Dict())
    data = DataGVRP(G´, [], [], [], [], [], 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0,0.0,0.0,0.0,0.0,[])
    
    graph = [(1, 1, 999999999.9),  (2, 2, 999999999.9),  (3, 3, 999999999.9), (4, 4, 999999999.9),
             (5, 5, 999999999.9), (6, 6, 999999999.9), (7, 7, 999999999.9),  (8, 8, 999999999.9),
             (9, 9, 999999999.9)]

    sepChar = ';'
    open(app["instance"]) do f
      data.max_d = data.max_f = data.max_t = 0.0
      data.min_d = data.min_f = data.min_t = 999999999.9
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

    """
    original read
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
    """

        for i in vertices(data)
      for j in vertices(data) # add arcs between vertices
        e = (i, j)
        if haskey(app, "preprocessings") && app["preprocessings"] != nothing
          push!(G´.E, e) # add edge e
          cost = 999999999.9
          if i < j && !((i, j) in data.E´)
            data.G´.cost[e] = distance(data, e)
            data.max_d  = ( data.max_d < d(data, ed(i, j)) ) ? d(data, ed(i, j)) : data.max_d
            data.max_f  = ( data.max_f < f(data, ed(i, j)) ) ? f(data, ed(i, j)) : data.max_f
            data.max_t  = ( data.max_t < t(data, ed(i, j)) ) ? t(data, ed(i, j)) : data.max_t

            data.min_d = ( data.min_d > d(data, ed(i, j)) ) ? d(data, ed(i, j)) : data.min_d
            data.min_f = ( data.min_f > f(data, ed(i, j)) ) ? f(data, ed(i, j)) : data.min_f
            data.min_t = ( data.min_t > t(data, ed(i, j)) ) ? t(data, ed(i, j)) : data.min_t
            cost = convert(Float64, t(data, ed(i, j)) ) + data.G´.V´[j].service_time
            push!(graph, (i, j, cost))
          else
            data.G´.cost[e] = 999999999.9
            push!(graph, (i, j, cost))
          end
          #append!( graph, [ (i, j , cost ) ] )
        elseif i < j
            push!(G´.E, e) # add edge e
            data.G´.cost[e] = distance(data, e)
            #cost = ceil( t(data, ed(i, j)) )
            #println(typeof( t(data, ed(i, j)) ))
            cost = convert(Float64, t(data, ed(i, j)) ) + data.G´.V´[j].service_time
            #println(typeof(cost))
            
            #cost2 = 9.0#(data, ed(i, j))
            
            push!(graph, (i, j, cost))
            #append!( graph, [ (i, j , cost ) ] )

            data.max_d  = ( data.max_d < d(data, ed(i, j)) ) ? d(data, ed(i, j)) : data.max_d
            data.max_f  = ( data.max_f < f(data, ed(i, j)) ) ? f(data, ed(i, j)) : data.max_f
            data.max_t  = ( data.max_t < t(data, ed(i, j)) ) ? t(data, ed(i, j)) : data.max_t

            data.min_d = ( data.min_d > d(data, ed(i, j)) ) ? d(data, ed(i, j)) : data.min_d
            data.min_f = ( data.min_f > f(data, ed(i, j)) ) ? f(data, ed(i, j)) : data.min_f
            data.min_t = ( data.min_t > t(data, ed(i, j)) ) ? t(data, ed(i, j)) : data.min_t
        end
      end
    end
    
    #Shortest path from a to e: a → c → d → e (cost 26)
    g = Digraph(graph)
    min_LB_E_j(data, g)
    
    #println("BFS ...")
    #println( "cost_BFS_path (1, 2): " , cost_BFS_path(data,1,2) )
    #println( "cost_BFS_path (1, 3): " , cost_BFS_path(data,1,3) )
    #println( "cost_BFS_path (1, 4): " , cost_BFS_path(data,1,4) )
    
    """
    C´ = Array{Int64}
    C´ = []
    push!(C´, 1 )
    for j in data.C
      push!(C´, j )
    end
    for j in C´
    push!( data.LB_E, 0.0 )
    end
    """

    println("DijkstraPath ...")
    src, dst = 4, 9
    path, cost = dijkstrapath(g, src, dst)
    println("Shortest path from $src to $dst: ", isempty(path) ? "no possible path" : join(path, " → "), " (cost $cost)")

    print_matrix(data, "time_cost")
    #print_matrix(data, "fuel_cost")
    #print_matrix(data, "distance_cost")
    return data
end

function read_Andelmin_Bartolini_Instance(app::Dict{String,Any})
    G′ = InputGraph([], [], Dict())
    data = DataGVRP(G′, 1, [], [], 0.0, 0.0, 0.0, 0.0)
    sepChar = ' '
    open(app["instance"]) do f
      # Ignore header and get paper default values
      readline(f)
      Customers     = 50
      Stations      = 22
      MaxTime       = 660
      MaxDistance   = 300
      Speed         = 0.66666667
      ServiceTime   = 30
      RefuelTime    = 15
      i = 1
      while !eof(f)
        # Read vertices
        line = readline(f)
        aux = split(line, [' ']; limit=0, keepempty=false)
        if length(aux) == 0
          break
        end
        
        v = Vertex(i, 0.0, 0.0, 0.0)
        v.pos_x = parse(Float64, aux[3])
        v.pos_y = parse(Float64, aux[4])
        push!(G′.V′, v)
        if aux[2] == "d"
          v.service_time = 0
        elseif aux[2] == "f"
          # Get AFS
          v.service_time = aux[2] == "f" ? RefuelTime : 0
          push!(data.F, v.id_vertex)
        elseif aux[2] == "c"
          # Get customer
          v.service_time = ServiceTime
          push!(data.C, v.id_vertex)
        end
        i += 1
      end
      # Read params
      
      #Vehicle capacity
      data.β = 60
      # Get vehicle fuel consumptin rate
      data.ρ = 0.2
      # Get vehicle time limit
      data.T = MaxTime
      # Get vehicle average speed
      data.ε = Speed

      line = readline(f)#
      #line = readline(f)#Infeasible customers 
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
        if i < j && !((i, j) in invalidEdges)
          e = (i, j)
          push!(G′.E, e) # add edge e
          data.G′.cost[e] = distance(data, e)
        end
      end
    end

    return data
end

#&& !((e[1], e[2]) in data.E´)
edges(data::DataGVRP) = data.G´.E # return set of arcs
d(data,e) = (e[1] != e[2] && !((e[1], e[2]) in data.E´) ) ? data.G´.cost[e] : 999999999.9 # cost of the edge e
f(data,e) = (e[1] != e[2] && !((e[1], e[2]) in data.E´) ) ? d(data, e) * data.ρ : 999999999.9 # fuel of the arc a 
t(data,e) = (e[1] != e[2] && !((e[1], e[2]) in data.E´) ) ? d(data, e) / data.ε : 999999999.9 # time of the arc a 
dimension(data::DataGVRP) = length(data.G´.V´) # return number of vertices
nb_vertices(data::DataGVRP) = length(vertices(data))

function get_LB_E(data::DataGVRP, index::Integer)
  j = 1
  value = 0.0
  for i in data.C 
    #value = data.LB_E[i]
    if i == index
      value = data.LB_E[j]
    end
    j += 1
    #for body
  end
  return value
end

function fill_with(x, value)
  x[:] .= value
end

# fix "bug"
# https://stackoverflow.com/questions/62550582/error-loaderror-boundserror-attempt-to-access-0-element-arraycandidate-1-a
#                           5 -> 4             2 -> 5
function dijsktra(data::DataGVRP, source::Int64, target::Int64)
  println("***********************************************************************")
  E = edges(data) # set of edges of the input graph G´
  n = nb_vertices(data)
  V = [i for i in 1:n] # set of vertices of the input graph G´
  V´= [i for i in 2:n]

  dist = Array{Float64}(undef, n)   # int dist[N]
  prev = Array{Float64}(undef, n)   # int prev[N]
  selected = zeros(Int64, n)        # selected[N]={0}
  fill_with(selected, 0)
  #selected = []
  #push!(selected, 0)
  path = Array{Int64}(undef, n)

  #i = m = min = start = d = j = 1   # index in Julia starts (for this implementation) in 1 not in 0 (like C/C++)
  #i::Int64      = 1 
  #m::Int64      = 1
  #min::Float64  = 1
  #start::Int64  = 1
  #d::Float64    = 1
  #j::Int64      = 1

  for i in V´
    dist[i] = 999999999
    prev[i] = -1
  end
  start = source
  selected[start] = 2               # second item
  dist[start] = 0                   # first item
  
  println(target," -> ", selected)
  println(start," -> ", selected[target])
  if selected[target] == 1 println("entra while") end
  while selected[target] == 0
    min = 999999999
    m = 1
    for i in V´
      c = 999999999
      if !( (start+1, i) in data.E´ || start+1 == i )
        c = data.G´.cost[ ed(start+1, i) ]
      end
      d = dist[start+1] + c

      if d < dist[i] && selected[i] == 1
        dist[i] = d
        prev[i] = start+1
      elseif min > dist[i] && selected[i] == 1
        min = dist[i]
        m = i
      end
    end
    start = m
    #  selected[start] = 1;
    println("*-> ", selected)
    selected[start] = 2
    #break
  end
  start = target
  j = 1
  while start != -1
        #j++
        path[j] = start#+65;
        j += 1
        start = prev[start]#;
  end

  #path[j]='\0';
  #strrev(path);
  #printf("%s", path);
  return dist[target];

  """
    int dist[N],prev[N],selected[N]={0},i,m,min,start,d,j;
    char path[N];
    for(i=1;i< N;i++)
    {
        dist[i] = IN;
        prev[i] = -1;
    }
    start = source;
    selected[start]=1;
    dist[start] = 0;
    while(selected[target] ==0)
    {
        min = IN;
        m = 0;
        for(i=1;i< N;i++)
        {
            d = dist[start] +cost[start][i];
            if(d< dist[i]&&selected[i]==0)
            {
                dist[i] = d;
                prev[i] = start;
            }
            if(min>dist[i] && selected[i]==0)
            {
                min = dist[i];
                m = i;
            }
        }
        start = m;
        selected[start] = 1;
    }
    start = target;
    j = 0;
    while(start != -1)
    {
        path[j++] = start+65;
        start = prev[start];
    }
    path[j]=' 0';
    strrev(path);
    printf("s", path);
    return dist[target];
  """
  #return 0.0
end

# min path cost i to j (satisfying problem restrictions)
function cost_BFS_path(data::DataGVRP, ii::Integer, jj::Integer)
  
  #(i -> j)
  if ii == jj
    return 0
  end

  """
  time_cost = t( data, ed( ii, jj ) )
  if time_cost <= data.T && f( data, ed( ii, jj ) ) <= data.β
    return time_cost
  elseif ii == jj 
    return 0 
  end
  """

  println(ii, " - ", jj)
  n = nb_vertices(data)
  V = data.F#[i for i in 1:n]
  E = edges(data)
  
  #---------------------------------------------------------
  #heurística p/ embaralhar a ordem de inicio da BFS ou algoritmo posterior
  s1 = V[ ii ]
  s2 = deepcopy(V)
  #popfirst!(s2) 
  filter!(e->e≠ii,s2)
  filter!(e->e≠jj,s2) 
  
  s3 = s2[randperm(length(s2))]
  #s3 = deepcopy(s2)                    # preserva a ordem 
  
  V´ = [ s1; [jj]; s3 ]
  print(BLUE_FG, ii, " ", jj, " ")
  stack = CrayonStack()
  print(stack, "[\t ")
  #print(stack, " ")
  for i in V´ 
    stack = CrayonStack()
    if i == ii || i == jj || i in data.F || i ==length(V´)+1
      if (i == ii || i == jj)
        print(BLUE_FG, "[", i, "]\t")
      elseif i in data.F && !(i == ii || i == jj)
        print(GREEN_FG, "[", i, "]\t")
      else
        print(stack, "[", i, "]\t")
      end
    else
      print(stack, i, "\t")
    end
  end
  stack = CrayonStack()
  println(stack, "] β: ",data.β," ","T: ",data.T)
  #---------------------------------------------------------
  
  s = V[1]                                                  #V[i] começa a BFS apartir do i até 1
  visited = [false for i in V]

  time_cost           = data.G´.V´[ ii ].service_time                  # começa zerado OU atente o i e depois começa a rota
  residual_fuel_cost  = data.β
  total_fuel_cost     = data.β                                  # veículo já sai cheio do depósito (otimizar: veículo sai com o suficiente p/ sair)
  total_distance_cost = 0.0

  #"""
  for c in V
    if visited[c]
      continue
    end
    visited[c] = true
    comp = [c]
    q = [c]
    #println(q)
    while length(q) > 0
      i = pop!(q)
      for j in V
        """
        # viável apenas em tempo
        if i != j && !((i, j) in data.E´) && !( time_cost > data.T)
          visited[j] = true
          push!(q, j)
          push!(comp, j)
          time_cost += t( data, ed( i, j ) ) + data.G´.V´[i].service_time # é um lower_bound, add o último e posteriorente viabilidade em combustível
          total_fuel_cost += f( data, ed( i, j ) )    # lower_bound, f( data, ed( i, j ) ) + residual
          obj_cost += data.G´.cost[ed(i, j)]
        end
        """
        #"""
        # viável em tempo && combustível
        if i != j && !((i, j) in data.E´) && ( time_cost < data.T && f( data, ed( i, j ) ) < data.β) || ( f( data, ed( i, j ) ) > residual_fuel_cost ) && !visited[j]
          visited[j] = true
          push!(q, j)
          push!(comp, j)
          time_cost += t( data, ed( i, j ) ) + data.G´.V´[j].service_time
          if i in data.F residual_fuel_cost = data.β end
          residual_fuel_cost  = residual_fuel_cost - f( data, ed( i, j ) )
          total_fuel_cost     += f( data, ed( i, j ) ) + residual_fuel_cost
          total_distance_cost += data.G´.cost[ed(i, j)]
        end
        #"""
      end
    end
    if (s in comp) && length(comp) > 1 
      println("min Route S' $ii -> $jj: ", comp) 
    end
  end
  #"""
  println("total_distance_cost: $total_distance_cost total_fuel_cost: $total_fuel_cost residual_fuel_cost: $residual_fuel_cost time_cost: $time_cost")
  #return total_fuel_cost
  #return residual_fuel_cost
  return time_cost
end

function min_LB_E_j(data::DataGVRP, g)
  #LB_E = Array{Float64}
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

  min_cost_e_jf1 = 999999999
  min_cost_e_jr1 = 999999999
  # se vc se refere ao  t_{f}^{'} da tese, eles são calculados pelo 
  # caminho minimo de f ao deposito no grafo induzido de AFSs com arestas de peso <= beta
  
  # bem se ele tem limite de tempo <= T e todas as arestas tem combustivel <= beta, 
  # então ele é viaável
  for j in C´
    for _f in data.F
      for _r in data.F
        #if _f < _r && ( _f!=1 && _r!=1 ) && !(( _f, _r ) in data.E´)
        if _f < _r && !(( _f, _r ) in data.E´)
        #if !(( _f, _r ) in data.E´)
          t_f = t_r = 0
          jj = 1

          src, dst = 1, _f
          path, cost_t = dijkstrapath(g, src, dst)
          #println("\n",isempty(path) ? "no possible path" : join(path, " → ")," len: ",length(path))
          #!( cost_t <= data.T && length(path)<1 )
          t_f = !( cost_t <= data.T) ? 999999999 : cost_t
          cost_fuel = (cost_t * data.ε * data.ρ ) / length(path)
          if t_f != 999999999
            t_f = !( cost_fuel <= data.β) ? 999999999 : cost_t 
            #t_f = cost_BFS_path(data, _f, jj)
          end

          src, dst = 1, _r
          path, cost_t = dijkstrapath(g, src, dst)
          #println("\n",isempty(path) ? "no possible path" : join(path, " → ")," len: ",length(path),"\n")
          #!( cost_t <= data.T && length(path)<1 )
          t_r = !( cost_t <= data.T) ? 999999999 : cost_t
          cost_fuel = ( cost_t * data.ε * data.ρ ) / length(path)
          if t_r != 999999999
            t_r = !( cost_fuel <= data.β) ? 999999999 : cost_t
            #t_r = cost_BFS_path(data, _r, jj)
          end

          if t_f + t( data, ed( _f, j ) ) + t( data, ed( j, _r ) ) + t_r <= data.T && f( data, ed( _f, j ) ) + f( data, ed( j, _r ) ) <= data.β
            if min_cost_e_jf1 > t_f
              min_cost_e_jf1 = t_f
              e_jf[1], e_jf[2] = j, _f
            elseif min_cost_e_jr1 > t_r
              min_cost_e_jr1 = t_r
              e_jr[1], e_jr[2] = j, _r
            end
          end
          #LB_E[i] = f(data, ed(i, j))
        end
      end
    end
    #push!( data.LB_E, 0.0 )
    #push!( data.LB_E, data.min_f )
    push!( data.LB_E, min( f( data, ed( e_jf[1], e_jf[2] ) ) , f( data, ed( e_jr[1], e_jr[2] ) ) ) )
    #push!( data.LB_E, min( min_cost_e_jf1 , min_cost_e_jr1 ) )
  end
  println(RED_FG, data.LB_E)
  stack = CrayonStack()
  println(stack,"...")
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

#--------------------------------------------------------------------------------------
struct Digraph{T <: Real,U}
    edges::Dict{Tuple{U,U},T}
    verts::Set{U}
end
 
function Digraph(edges::Vector{Tuple{U,U,T}}) where {T <: Real,U}
    vnames = Set{U}(v for edge in edges for v in edge[1:2])
    adjmat = Dict((edge[1], edge[2]) => edge[3] for edge in edges)
    return Digraph(adjmat, vnames)
end
 
vertices(g::Digraph) = g.verts
edges(g::Digraph)    = g.edges
 
neighbours(g::Digraph, v) = Set((b, c) for ((a, b), c) in edges(g) if a == v)
 
function dijkstrapath(g::Digraph{T,U}, source::U, dest::U) where {T, U}
    @assert source ∈ vertices(g) "$source is not a vertex in the graph"
 
    # Easy case
    if source == dest return [source], 0 end
    # Initialize variables
    inf  = typemax(T)
    dist = Dict(v => inf for v in vertices(g))
    prev = Dict(v => v   for v in vertices(g))
    dist[source] = 0
    Q = copy(vertices(g))
    neigh = Dict(v => neighbours(g, v) for v in vertices(g))
 
    # Main loop
    while !isempty(Q)
        u = reduce((x, y) -> dist[x] < dist[y] ? x : y, Q)
        pop!(Q, u)
        if dist[u] == inf || u == dest break end
        for (v, cost) in neigh[u]
            alt = dist[u] + cost
            if alt < dist[v]
                dist[v] = alt
                prev[v] = u
            end
        end
    end
 
    # Return path
    rst, cost = U[], dist[dest]
    if prev[dest] == dest
        return rst, cost
    else
        while dest != source
            pushfirst!(rst, dest)
            dest = prev[dest]
        end
        pushfirst!(rst, dest)
        return rst, cost
    end
end
#--------------------------------------------------------------------------------------

#testgraph = [("a", "b", 7),  ("a", "c", 9),  ("a", "f", 14), ("b", "c", 10),
#             ("b", "d", 15), ("c", "d", 11), ("c", "f", 2),  ("d", "e", 6),
#             ("e", "f", 9)]

# Print all possible paths
# @printf("\n%4s | %3s | %s\n", "src", "dst", "path")
# @printf("----------------\n")
# for src in vertices(g), dst in vertices(g)
#    path, cost = dijkstrapath(g, src, dst)
#    @printf("%4s | %3s | %s\n", src, dst, isempty(path) ? "no possible path" : join(path, " → ") * " ($cost)")
# end