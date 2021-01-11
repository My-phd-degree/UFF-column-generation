mutable struct Solution
  cost::Union{Int,Float64}
  routes::Array{Array{Int}}
end

# build Solution from the variables x
function getsolution(data::DataGVRP, optimizer::VrpOptimizer, x, objval, app::Dict{String,Any})
  E, dim = edges(data), dimension(data)
  n = nb_vertices(data)
  V = [i for i in 1:n]
  T = data.T
  K = data.M
  adj_list = [[] for i in 1:dim]
  
  println("MATRIX\n\n")
  
  print("\t")
  for j in V
    print(j, "\t")
  end
  println("\n")
  for i in V
    print(i, "\t")
    for j in V
      if i != j && !((i, j) in data.E′) 
        print(data.G′.cost[ed(i, j)], "\t")
      else
        print("*.*", "\t")
      end
    end
    println("")
  end

  for k in K
    print("Route $k S: ", k)
    for i in V
      for j in V
        if i != j && !((i, j) in data.E′) 
          val = get_value(optimizer, x[i, j, k])
          if val > 0.5
            push!(adj_list[i], j)
            push!(adj_list[j], i)
            if val > 1.5
              push!(adj_list[i], j)
              push!(adj_list[j], i)
            end
            print("( $i - $j ) ", i, j)
          end
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
