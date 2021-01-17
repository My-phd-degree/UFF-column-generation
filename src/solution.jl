#import Pkg; 
#Pkg.add("Distances")
#Pkg.add("Crayons")
#Pkg.add("Plots")
#Pkg.add("GraphRecipes")
#Pkg.add("LightGraphs")

using Crayons                   # https://github.com/KristofferC/Crayons.jl
using Crayons.Box               # ...
using GraphRecipes, Plots       # https://docs.juliaplots.org/latest/tutorial/#tutorial
using LightGraphs               # https://docs.juliaplots.org/latest/graphrecipes/examples/

mutable struct Solution
  cost::Union{Int,Float64}
  routes::Array{Array{Int}}
end

# build Solution from the variables x
function getsolution(data::DataGVRP, optimizer::VrpOptimizer, x, e, objval, app::Dict{String,Any})
  E, dim = edges(data), dimension(data)
  n = nb_vertices(data)
  V = [i for i in 1:n]
  T = data.T
  C = data.C
  K = data.M
  log = 1

  F´ = deepcopy(data.F)
  popfirst!(F´)

  adj_list = [[] for i in 1:dim]
  

  println("\nMatrix Costs:\n")
  
  print("\t")
  for j in V
    stack = CrayonStack()
    if j in F´
      print(GREEN_FG, j, "\t")
    else
      print(stack, j, "\t")
    end
  end
  println("\n")
  for i in V
    stack = CrayonStack()
    if i in F´
      print(GREEN_FG, i, "\t")
    else
      print(stack, i, "\t")
    end
    for j in V
      if i != j && !((i, j) in data.E´) 
        stack = CrayonStack()
        print(stack, data.G´.cost[ed(i, j)], "\t")
      else
        print(RED_FG, "***", "\t")
        stack = CrayonStack()
      end
    end
    println("")
  end
  println("")

  println(data.LB_E)
  for k in K
    stack = CrayonStack()
    print(stack, "Route $k S: ", k)
    visited = [false for i in V]
    #print("[ ")  
    #for lb in data.LB_E
    #  print(lb," ")
    #end
    
    #print("]")

    for i in V
      for j in V
        #if visited[j] continue end
        if i != j && !((i, j) in data.E´) #&& !visited[j]
          val = get_value(optimizer, x[i, j, k])
          #####################################
          if val > 0.5
            push!(adj_list[i], j)
            push!(adj_list[j], i)
            if val > 1.5
              push!(adj_list[i], j)
              push!(adj_list[j], i)
            end
            if i in F´
              print("( ")
              print(GREEN_FG, i)
              stack = CrayonStack()
              print(stack, " - ",j," ) ")
            elseif j in F´ 
              stack = CrayonStack()
              print(stack,"( ", i, " - ")
              print(GREEN_FG, j)
              print(stack, " ) ")
            else
              print("( $i - $j ) ", i, j)
            end
            #i = j
            #visited[j] = true
            #continue
          end
          #####################################
        end
      end
    end
    println("")
  end
  for i in 1:length(adj_list)
    #print(i, ": ", adj_list[i], "\n")
  end
end

function print_routes(solution)
  #x = 1:10; y = rand(10); # These are the plotting data
  #plot(x, y)
  #savefig("myplot.png")
  """
  n = 8
  g = wheel_digraph(n)
  edgelabel_dict = Dict()
  edgelabel_mat = Array{String}(undef, n, n)
  for i in 1:n
      for j in 1:n
          edgelabel_mat[i, j] = edgelabel_dict[(i, j)] = string("edge ", i, " to ", j)
      end
  end
  edgelabel_vec = edgelabel_mat[:]

  graphplot(g, names=1:n, edgelabel=edgelabel_dict, curves=false, nodeshape=:rect)  # Or edgelabel=edgelabel_mat, or edgelabel=edgelabel_vec.
  """
end

# checks the feasiblity of a solution
ed(i, j) = i < j ? (i, j) : (j, i)

function checksolution(data::DataGVRP, solution)
end

# read solution from file (CVRPLIB format)
function readsolution(app::Dict{String,Any})
end

# write solution in a file
function writesolution(solpath, solution)
end

# write solution as TikZ figure (.tex)
function drawsolution(tikzpath, data, solution)
end
