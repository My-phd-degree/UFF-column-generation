function get_invalid_edges_1(data::DataGVRP) 
  edges = []
  for i in 1:length(data.G′.V′)
    #    for j in (i + 1):length(data.G′.V′)
    for j in 1:length(data.G′.V′)
      if i != j && (f(data, ed(i, j)) > data.β || t(data, ed(i, j)) + (data.G′.V′[i].service_time + data.G′.V′[j].service_time) > data.T)
        push!(edges, ed(i, j))
      end
    end
  end
  return edges
end

function get_invalid_edges_2(data::DataGVRP) 
  gvrp_afs_tree = data.gvrp_afs_tree
  edges = []
  fuel = f
  for i in data.C 
    if fuel(data, ed(i, data.depot_id)) > data.β/2.0
      valid = false
      for f in gvrp_afs_tree.F0
        if fuel(data, ed(f, i)) + fuel(data, ed(i, data.depot_id)) <= data.β && 
          gvrp_afs_tree.times[f] + t(data, ed(f, i)) + t(data, ed(i, data.depot_id)) <= data.T
          valid = true
          break
        end
      end
      if !valid
        push!(edges, ed(i, data.depot_id))
      end
    end
  end
  return edges
end

function get_invalid_edges_3(data::DataGVRP) 
  gvrp_afs_tree = data.gvrp_afs_tree
  edges = []
  fuel = f
  for i in data.C 
    for j in data.C 
      if i < j
        invalidEdge = true
        #check if edge (i, j) can be feasible
        for f in gvrp_afs_tree.F0 
          if f == data.depot_id || gvrp_afs_tree.pred[f] != f
            for r in gvrp_afs_tree.F0 
              #if the afss f and r are connected 
              if (r == data.depot_id || gvrp_afs_tree.pred[r] != r) &&
                fuel(data, ed(f, i)) + fuel(data, ed(i, j)) + fuel(data, ed(j, r)) <= data.β && 
                gvrp_afs_tree.times[f] + t(data, ed(f, i)) + t(data, ed(i, j)) + t(data, ed(j, r)) + gvrp_afs_tree.times[r] <= data.T 
                invalidEdge = false
                break
              end
            end
            if !invalidEdge
              break
            end
          end
        end
        if invalidEdge
          push!(edges, ed(i, j))
        end
      end
    end
  end
  return edges
end

function get_invalid_edges_4(data::DataGVRP) 
  gvrp_afs_tree = data.gvrp_afs_tree
  edges = []
  fuel = f
  for i in data.C 
    for f in gvrp_afs_tree.F0 
      invalidEdge = true
      if f == data.depot_id || gvrp_afs_tree.pred[f] != f 
        for r in gvrp_afs_tree.F0 
          if (r == data.depot_id || gvrp_afs_tree.pred[r] != r) &&
            fuel(data, ed(f, i)) + fuel(data, ed(i, r)) <= data.β && 
            gvrp_afs_tree.times[f] + t(data, ed(f, i)) + t(data, ed(i, r)) + gvrp_afs_tree.times[r] <= data.T 
              invalidEdge = false
              break
          end
        end
      end
      if invalidEdge
        push!(edges, ed(i, f))
      end
    end
  end
  return edges
end
