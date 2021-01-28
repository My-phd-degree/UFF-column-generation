mutable struct Solution
  cost::Union{Int,Float64}
  routes::Array{Array{Int}}
end

ed(i, j) = i < j ? (i, j) : (j, i)

# build Solution from the variables x
function getsolution(data::DataGVRP, optimizer::VrpOptimizer, x, objval, app::Dict{String,Any})
  objval = 0.0
  routes = []
  E, dim = edges(data), dimension(data)

#  println("Num of paths $(get_number_of_positive_paths(optimizer))")
  for i in 1:get_number_of_positive_paths(optimizer)
    λ = get_value(optimizer, i)
    if λ >= 0.001 && λ <= 1 - 0.001
      println("Frac λ with $λ")
    end

    #valid λ
    if λ > 1 - 0.001
      #get λ used edges
      vars = [x[e] for e in E]
      E_vals = get_values(optimizer, vars, i)
      E_weights = Dict{Tuple{Int64, Float64}, Int64}()
      #get route
      for i in 1:length(E)
        e = E[i]
        e_val = E_vals[i]
        #        if e_val > 0.001
        #          println("$(e) : $(e_val)")
        #        end
        E_weights[e] = e_val
      end

      #get route
      function getRoute(i::Int64, E_weights::Dict{Tuple{Int64, Float64}, Int64}, route::Array{Int64})
        if i == data.depot_id
          feasible = true
          for (e, weight) in E_weights
            if weight > 0.001
              feasible = false
            end
          end  
          return feasible
        end
        for e in δ(data, i)
          if E_weights[e] > 0.001
            j = e[1] == i ? e[2] : e[1]
            E_weights[e] = E_weights[e] - 1
            push!(route, j)
            if getRoute(j, E_weights, route)
              return true
            end
            E_weights[e] = E_weights[e] + 1
            pop!(route)
          end
        end
        return false
      end

      route = [data.depot_id]
      for e in δ(data, data.depot_id)
        if E_weights[e] > 0.001
          j = e[1] == data.depot_id ? e[2] : e[1]
          E_weights[e] = E_weights[e] - 1
          push!(route, j)
          getRoute(j, E_weights, route)
          break
        end
      end
      push!(routes, route)
    end
  end
  
  for route in routes
    n = length(route)
    for i in 2:n
      objval = objval + d(data, ed(route[i - 1], route[i]))
    end
  end
  return Solution(objval, routes)
end

function print_routes(data::DataGVRP, solution::Solution)
  for (i, r) in enumerate(solution.routes)
    print("Route #$i: ")
    for j in r
      print("$(data.G′.V′[j].id_vertex) ")
    end
    println()
  end
end

# checks the feasiblity of a solution

function checksolution(data::DataGVRP, solution)
  nTimesCustomersVisited = Dict{Int64, Int64}([i => 0 for i in data.C])
  routes = solution.routes
  for route in routes
    n = length(route)
    n == 0 && error("Route is empty")
    route[1] != 1 && error("Route does not begins at the depot")
    route[n] != 1 && error("Route does not ends at the depot")
    consumedFuel = consumedTime = 0
    for i in 2:n
      a, b = route[i - 1], route[i]
      #update resources
      consumedFuel = consumedFuel + f(data, ed(a, b))
      consumedTime = consumedTime + t(data, ed(a, b))
      #error checking
      consumedFuel > data.β && error("No fuel in $b: $consumedFuel")
      consumedTime > data.T && error("Time exceeded in $b: $consumedTime")
      # restore fuel
      if b in data.F
        consumedFuel = 0
      elseif b in data.C
        nTimesCustomersVisited[b] = nTimesCustomersVisited[b] + 1
        nTimesCustomersVisited[b] > 1 && error("Customer $b visited more than once")
      end
    end
  end
end

# read solution from file (CVRPLIB format)
function readsolution(app::Dict{String,Any})
  str = read(app["sol"], String)
  breaks_in = [' '; ':'; '\n';'\t';'\r']
  aux = split(str, breaks_in; limit=0, keepempty=false)
  sol = Solution(0, [])
  j = 3
  while j <= length(aux)
    r = []
    while j <= length(aux)
      push!(r, parse(Int, aux[j]))
      j += 1
      if contains(lowercase(aux[j]), "cost") || contains(lowercase(aux[j]), "route")
        break
      end
    end
    push!(sol.routes, r)
    if contains(lowercase(aux[j]), "cost")
      if app["noround"]
        sol.cost = parse(Float64, aux[j + 1])
      else
        sol.cost = parse(Int, aux[j + 1])
      end
      return sol
    end
    j += 2 # skip "Route" and "#j:" elements
  end
  error("The solution file was not read successfully.")
  return sol
end

# write solution in a file
function writesolution(solpath::String, data::DataGVRP, solution::Solution)
  open(solpath, "w") do f
    for (i, r) in enumerate(solution.routes)
      write(f, "Route #$i: ")
      for j in r
        write(f, "$(data.G′.V′[j].id_vertex) ")
      end
      write(f, "\n")
    end
    write(f, "Cost $(solution.cost)\n")
  end
end

function drawsolution(tikzpath::String, data::DataGVRP, routes::Array{Array{Int64}})

  open(tikzpath, "w") do f
    write(f,"\\documentclass[crop,tikz]{standalone}\n\\begin{document}\n")
    write(f,"\\usetikzlibrary{arrows,positioning,automata,shadows,fit,shapes,calc,shapes.geometric}
          \\tikzset{triangle_black/.style={regular polygon, regular polygon sides=3, minimum size=0.3cm, fill=black}}
          \\tikzset{triangle/.style={regular polygon, regular polygon sides=3, minimum size=0.3cm, fill=white}}
          \\tikzset{square/.style={regular polygon, regular polygon sides=4, minimum size=0.3cm, fill=white}}
          \\tikzset{circufe/.style={circle,draw, minimum size=0.2cm, fill=white}}
          \\tikzset{raio/.style={circle,draw, minimum size=1.45cm, fill=white}} 
          \\tikzset{circle_new/.style={circle,draw, minimum size=0.2cm, fill=white}}   
          ")
    # get limits to draw
    pos_x_vals = [i.pos_x for i in data.G′.V′]
    pos_y_vals = [i.pos_y for i in data.G′.V′]
    scale_fac = 1/(max(maximum(pos_x_vals),maximum(pos_y_vals))/10)
    # draw
    write(f,"\\begin{tikzpicture}\n\\tikzstyle{every node}=[circle, draw, fill=black!50, inner sep=0pt, minimum width=4pt]")
    # nodes
    # depot
    for i in 1:length(data.G′.V′)
      type = ""
      id_vertex = data.G′.V′[i].id_vertex
      if i == data.depot_id
        type = "triangle_black"
      elseif i in data.C
        type = "square"
      else
        type = "triangle"
      end
      write(f, "\t\\node[$(type)] (v$(id_vertex)) at ($(scale_fac*pos_x_vals[i]),$(scale_fac*pos_y_vals[i])) {\\footnotesize $(id_vertex)};\n")
    end
    # define colors
    colors = [
              "black", 
              "red", 
              "green", 
              "blue", 
              "cyan", 
              "magenta", 
              "yellow",
              "orange",
              "purple",
              "black!50!gray", 
              "red!50!gray", 
              "green!50!gray", 
              "blue!50!gray", 
              "cyan!50!gray", 
              "magenta!50!gray", 
              "yellow!50!gray",
              "orange!50!gray",
              "purple!50!gray",
             ]
    color_index = 1
    # get edges in F_0 visited multiple times
    F_0_edges = Dict{Tuple{Int64,Int64}, Int64}()
    F_0 = vcat(data.F, [data.depot_id])
#    println("F_0:")
#   [print("$(data.G′.V′[f].id_vertex) ") for f in F_0]
#   println()
#   [print("$(f) ") for f in F_0]
#   println("\nC:")
#   [print("$(data.G′.V′[f].id_vertex) ") for f in data.C]
#   println()
#   [print("$(f) ") for f in data.C]
#   println()
    for (i, j) in data.G′.E
      if i in F_0 && j in F_0
        # get number of times this edge is visited
        count = 0
        for route in routes
          for k in 2:length(route)
            if ed(route[k - 1], route[k]) == (i, j)
              count = count + 1
            end
          end
        end
        #draw edge
        color_index = 1
        for route in routes
          for k in 2:length(route)
            if ed(route[k - 1], route[k]) == (i, j)
              write(f, "\t\\path [$(colors[color_index]), draw,-latex] (v$(data.G′.V′[route[k - 1]].id_vertex)) edge [bend left = $(count*20)] (v$(data.G′.V′[route[k]].id_vertex));\n")
#              println("($i, $j)) $(data.G′.V′[route[k - 1]].id_vertex) $(data.G′.V′[route[k]].id_vertex)")
              count = count - 1
            end
          end
          color_index = color_index + 1
        end
      end
    end
    # get edges
    color_index = 1
    for (r, route) in enumerate(routes)
      n = length(route)
      for i in 2:n
        e = (data.G′.V′[route[i - 1]].id_vertex, data.G′.V′[route[i]].id_vertex)
#        println(e, " ($(route[i - 1]), $(route[i]))")
        if !(route[i - 1] in F_0 && route[i] in F_0)
          write(f, "\t\\path [$(colors[color_index]), draw,-latex] (v$(e[1])) edge (v$(e[2]));\n")
        end
      end
      color_index = color_index + 1
    end
    write(f, "\\end{tikzpicture}\n")
    write(f, "\\end{document}\n")
  end   
end
