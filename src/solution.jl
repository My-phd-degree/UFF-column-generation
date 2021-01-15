mutable struct Solution
  cost::Union{Int,Float64}
  routes::Array{Array{Int}}
end

function get_routes(adj_matrix::Array, data::DataGVRP)
end

# build Solution from the variables x
function getsolution(data::DataGVRP, optimizer::VrpOptimizer, x, objval, app::Dict{String,Any})
  # adj list
  E, dim = edges(data), dimension(data)
  adj_list = [[] for i in 1:dim]
  for e in E 
    val = get_value(optimizer, x[e])
    if val > 0.5
      push!(adj_list[e[1]], e[2])
      push!(adj_list[e[2]], e[1])
      if val > 1.5
        push!(adj_list[e[1]], e[2])
        push!(adj_list[e[2]], e[1])
      end
    end
  end
  # DFS
  visited, routes = [false for i in 1:dim], []
  """
  function dfs(f::Int64, adj_list::Array{Array{Int64}})
    visited[f] = true
  end
  for f in data.F
    dfs(f, adj_list)
  end
  """

  print(adj_list)
  for i in 1:length(adj_list)
    print(i, ": ", adj_list[i], "\n")
  end
  # for j in adj_list
  i = 1 #j[1]
  first = i
  if !visited[i]
    r, prev = [], first
    push!(r, i)
    visited[i] = true
    length(adj_list[i]) != 2 && i in data.C && error("Problem trying to recover the route from the x values. " *
                                                     "Customer $i has $(length(adj_list[i])) incident edges.")
    next, prev = (adj_list[i][1] == prev) ? adj_list[i][2] : adj_list[i][1], i
    maxit, it = dim, 0
    print("($prev, $next), ")
    while next != first && it < maxit
      length(adj_list[next]) != 2 && next in data.C && error("Problem trying to recover the route from the x values. " *
                                                             "Customer $next has $(length(adj_list[next])) incident edges.")
      push!(r, next)
      visited[next] = true
      aux = next
      next, prev = (adj_list[next][1] == prev) ? adj_list[next][2] : adj_list[next][1], aux
      it += 1
      print("($prev, $next), ")
    end
    print("\n")
    push!(routes, r)
  end
  # end
  # if !app["noround"]
  # objval = trunc(Int, round(objval))
  # end
  return Solution(objval, routes)
end

function print_routes(solution)
  for (i, r) in enumerate(solution.routes)
    print("Route #$i: ")
    for j in r
      print("$j ")
    end
    println()
  end
end

# checks the feasiblity of a solution
ed(i, j) = i < j ? (i, j) : (j, i)

function checksolution(data::DataGVRP, solution)
  dim, T, β = dimension(data), data.T, data.β
  visits = [0 for i in 1:dim]
  sum_cost = 0.0
  for (i, r) in enumerate(solution.routes)
    sum_time, sum_fuel, prev = 0.0, 0.0, r[1]
    visits[r[1]] += 1
    println(r)
    for j in r[2:end]
      visits[j] += 1
      j in data.C && visits[j] == 2 && error("Customer $j was visited more than once")
      sum_cost += d(data, ed(prev, j))
      sum_time += t(data, ed(prev, j))
      sum_fuel = (prev in data.F) ? 0.0 : sum_fuel
      sum_fuel += f(data, ed(prev, j))
      println(j)
      (sum_time > T) && error("Route $r is violating the limit T. Total time spent is at least $(sum_time) and T is $T")
      (sum_fuel > β) && error("Route is violating the limit β. Total fuel spent is at least $(sum_fuel) and β is $β")
      prev = j
    end
    if prev != r[1] 
      sum_cost += d(data, ed(prev, r[1])) 
      sum_fuel += f(data, ed(prev, r[1]))
    end
    print("Route $i (time: $sum_time, fuel: $sum_fuel)\n")
  end
  !isempty(filter(a -> a == 0, visits)) && error("The following vertices were not visited: $(filter(a -> a == 0, visits))")
  (abs(solution.cost - sum_cost) > 0.001) && error("Cost calculated from the routes ($sum_cost) is different from that passed as" *
    " argument ($(solution.cost)).")
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
function writesolution(solpath, solution)
  open(solpath, "w") do f
    for (i, r) in enumerate(solution.routes)
      write(f, "Route #$i: ")
      for j in r
        write(f, "$j ")
      end
      write(f, "\n")
    end
    write(f, "Cost $(solution.cost)\n")
  end
end

# write solution as TikZ figure (.tex)
function drawsolution(tikzpath, data, solution)
  open(tikzpath, "w") do f
    write(f, "\\documentclass[crop,tikz]{standalone}\n\\begin{document}\n")
    # get limits to draw
    pos_x_vals = [i.pos_x for i in data.G′.V′]
    pos_y_vals = [i.pos_y for i in data.G′.V′]

    if data.insType == "EUC_2D" || data.insType == "GEO"

      if data.insType == "EUC_2D"
        scale_fac = 1 / (max(maximum(pos_x_vals), maximum(pos_y_vals)) / 10)
      elseif data.insType == "GEO"
        scale_fac = 1 / (max(maximum(pos_x_vals), maximum(pos_y_vals)) / 10)
      end

      write(f, "\\begin{tikzpicture}[thick, scale=1, every node/.style={scale=0.3}]\n")
      for i in data.G′.V′
        x_plot = scale_fac * i.pos_x
        y_plot = scale_fac * i.pos_y
        if i.id_vertex in data.F # plot balck vertices
          write(f, "\t\\node[draw, line width=0.1mm, circle, fill=black, inner sep=0.05cm, text=white] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {\\footnotesize $(i.id_vertex)};\n")
        else
          write(f, "\t\\node[draw, line width=0.1mm, circle, fill=white, inner sep=0.05cm] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {\\footnotesize $(i.id_vertex)};\n")
        end
      end

      for r in solution.routes
        prev = r[1]
        first = r[1]
        edge_style = "-,line width=0.8pt"
        for i in r[2:end]
          e = (prev, i)
          write(f, "\t\\draw[$(edge_style)] (v$(e[1])) -- (v$(e[2]));\n")
          prev = i
        end
        write(f, "\t\\draw[$(edge_style)] (v$(first)) -- (v$(prev));\n")
      end

      write(f, "\\end{tikzpicture}\n")
      write(f, "\\end{document}\n")
    else
      println("Draw funciton for $(data.insType) not defined yet!")
    end
  end
end
