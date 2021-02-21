function calculateGvrpLBByImprovedMST(S₀::Array{Int64}, η::Dict{Int64,Float64}, π::Dict{Int64,Float64}, gvrpReducedGraph::Dict{Tuple{Int64,Int64}, Float64}) 

end

int utils::calculateGVRP_BPP_NRoutesLB(const Gvrp_instance& gvrp_instance, const vector<const Vertex *>& vertices, const vector<pair<double, double>>& closestsTimes, unsigned int execution_time_limit) {
  const int svertices = vertices.size();
  vector<double> items (svertices); 
  for (int i = 0; i < svertices; ++i) 
    items[i] = vertices[i]->serviceTime + (closestsTimes[i].first + closestsTimes[i].second)/2.0;
  BPP_instance bpp_instance (items, gvrp_instance.timeLimit);
  BPP_model bpp_model (bpp_instance, execution_time_limit);
  return bpp_model.run().first.size();
}
