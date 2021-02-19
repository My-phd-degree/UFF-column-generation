pair<vector<vector<double>>, vector<vector<double>>> utils::calculateGVRPReducedGraphs (const Gvrp_instance& gvrp_instance, const Gvrp_afs_tree& gvrp_afs_tree) {
  int depotId = gvrp_instance.depot.id;
  const int sc0 = gvrp_instance.customers.size() + 1.0,
              stotal = gvrp_instance.distances.size();
  vector<const Vertex *> c0 (sc0);
  vector<vector<double>> closestTime (stotal, vector<double>(stotal, -1.0)); 
  vector<vector<double>> closestDistance (stotal, vector<double>(stotal, -1.0)); 
  //build c0
  c0[0] = &gvrp_instance.depot;
  int i = 1;
  for (const Vertex& customer : gvrp_instance.customers) {
    c0[i] = &customer;
    ++i;
  }
  //for each pair
  for (int i = 0; i < sc0; ++i) {
    const Vertex * customerI =  c0[i];
    int customerIid = customerI->id;
    for (int j = i + 1; j < sc0; ++j) {
      const Vertex * customerJ =  c0[j];
      int customerJid = customerJ->id;
      if (customerIid != customerJid && 
          gvrp_instance.time(depotId, customerIid) + customerI->serviceTime + gvrp_instance.time(customerIid, customerJid) + customerJ->serviceTime + gvrp_instance.time(customerIid, depotId) <= gvrp_instance.timeLimit) {
        //origin afs
        for (int f = 0; f < gvrp_afs_tree.f0.size(); f++) {
          int fId = gvrp_afs_tree.f0[f]->id;
          double spentTime = gvrp_afs_tree.times[f] + gvrp_instance.time(fId, customerIid) + customerI->serviceTime;
          double spentFuel = gvrp_instance.fuel(fId, customerIid);
          if (spentFuel <= gvrp_instance.vehicleFuelCapacity && 
              spentTime + gvrp_instance.time(customerIid, customerJid) + customerJ->serviceTime + gvrp_instance.time(customerJid, depotId) <= gvrp_instance.timeLimit) {
            //destiny afs
            for (int r = 0; r < gvrp_afs_tree.f0.size(); r++) {
              int rId = gvrp_afs_tree.f0[r]->id;
              if (spentTime + gvrp_instance.time(customerIid, customerJid) + customerJ->serviceTime + gvrp_instance.time(customerJid, rId) + gvrp_afs_tree.times[r] <= gvrp_instance.timeLimit) {
                //direct travel
                if (spentFuel + gvrp_instance.fuel(customerIid, customerJid) + gvrp_instance.fuel(customerJid, rId) <= gvrp_instance.vehicleFuelCapacity) {
                  double time = gvrp_instance.time(customerIid, customerJid);
                  double cost = gvrp_instance.distances[customerIid][customerJid];
                  if (closestTime[customerIid][customerJid] < 0.0 || time < closestTime[customerIid][customerJid]) 
                    closestTime[customerIid][customerJid] = closestTime[customerJid][customerIid] = time;
                  if (closestDistance[customerIid][customerJid] < 0.0 || cost < closestDistance[customerIid][customerJid]) 
                    closestDistance[customerIid][customerJid] = closestDistance[customerJid][customerIid] = cost;
                } else {
                  //intermediate travel
                  for (int f_ = 0; f_ < gvrp_afs_tree.f0.size(); f_++) {
                    int f_Id = gvrp_afs_tree.f0[f_]->id;
                    //fuel and time feasible
                    if (spentFuel + gvrp_instance.fuel(customerIid, f_Id) <= gvrp_instance.vehicleFuelCapacity && 
                        spentTime + gvrp_instance.time(customerIid, f_Id) + gvrp_afs_tree.f0[f_]->serviceTime + gvrp_instance.time(f_Id, customerJid) + customerJ->serviceTime + gvrp_instance.time(customerJid, rId) + gvrp_afs_tree.times[r] <= gvrp_instance.timeLimit) {
                      //intermediate afs destiny
                      for (int r_ = 0; r_ < gvrp_afs_tree.f0.size(); r_++) {
                        int r_Id = gvrp_afs_tree.f0[r_]->id;
                        //time feasible
                        if (gvrp_instance.fuel(r_Id, customerJid) + gvrp_instance.fuel(customerJid, rId) <= gvrp_instance.vehicleFuelCapacity && 
                            spentTime + gvrp_instance.time(customerIid, f_Id) + gvrp_afs_tree.pairTimes[f_][r_] + gvrp_instance.time(r_Id, customerJid) + customerJ->serviceTime + gvrp_instance.time(customerJid, rId) + gvrp_afs_tree.times[r] <= gvrp_instance.timeLimit) {
                          double time = gvrp_instance.time(customerIid, f_Id) + gvrp_afs_tree.pairTimes[f_][r_] + gvrp_instance.time(r_Id, customerJid);
                          double cost = gvrp_instance.distances[customerIid][f_Id] + gvrp_afs_tree.pairCosts[f_][r_] + gvrp_instance.distances[r_Id][customerJid];
                          if (closestTime[customerIid][customerJid] < 0.0 || 
                              time < closestTime[customerIid][customerJid]) 
                            closestTime[customerIid][customerJid] = closestTime[customerJid][customerIid] = time;
                          if (closestDistance[customerIid][customerJid] < 0.0 || 
                              cost < closestDistance[customerIid][customerJid]) 
                            closestDistance[customerIid][customerJid] = closestDistance[customerJid][customerIid] = cost;
                        }
                      }
                    } 
                  }
                }
              } 
            }
          }
        }
      }
    }
  } 
  //set unfinded distances
  for (int i = 0; i < stotal; ++i) 
    closestDistance[i][i] = closestTime[i][i] = 0.0;
  /*
  for (int i = 1; i < sc0; ++i)
    for (int j = 1; j < sc0; ++j) {
      int iId = c0[i]->id,
          jId = c0[j]->id;
      if (closestDistance[iId][jId] < 0.0) {
        closestDistance[iId][jId] = closestDistance[iId][depotId] + closestDistance[depotId][jId];
        closestTime[iId][jId] = closestTime[iId][depotId] + closestTime[depotId][jId];
      }
    }
    */
  return {closestDistance, closestTime};
}
int utils::calculateGVRP_BPP_NRoutesLB(const Gvrp_instance& gvrp_instance, const vector<const Vertex *>& vertices, const vector<pair<double, double>>& closestsTimes, unsigned int execution_time_limit) {
  const int svertices = vertices.size();
  vector<double> items (svertices); 
  for (int i = 0; i < svertices; ++i) 
    items[i] = vertices[i]->serviceTime + (closestsTimes[i].first + closestsTimes[i].second)/2.0;
  BPP_instance bpp_instance (items, gvrp_instance.timeLimit);
  BPP_model bpp_model (bpp_instance, execution_time_limit);
  return bpp_model.run().first.size();
}
double utils::calculateGvrpLBByImprovedMST (const vector<const Vertex *>& vertices, const vector<pair<double, double>> closests, const vector<vector<double>>& gvrpReducedGraph) {
  const int svertices = vertices.size();
  const Vertex * nodeI, * nodeJ, * nodeK;
  double bestLB = 0.0, cost;
  //boruvka
  vector<vector<bool>> adjMatrix (svertices, vector<bool> (svertices, false));
  DSU dsu(svertices);
  int nTrees = svertices;
  double MSTCost = 0.0;
  int setI, setJ, setK;
//  cout<<"building mst"<<endl;
  vector<pair<int, int>> bestEdge (svertices, {-1, -1});
  vector<double> bestEdgeCost (svertices, DBL_MAX);
  while (nTrees > 1) {
//    cout<<"\tnTrees: "<<nTrees<<endl;
    //calculate best cuts
    for (int i = 0; i < svertices; ++i) {
      nodeI = vertices[i];
      setI = dsu.findSet(i);
      for (int j = 0; j < svertices; ++j) {
        nodeJ = vertices[j];
        cost = gvrpReducedGraph[nodeI->id][nodeJ->id];
        if (i != j && setI != dsu.findSet(j) && cost < bestEdgeCost[setI]) {
          bestEdgeCost[setI] = cost;
          bestEdge[setI] = {i, j};
        }
      }
    }
    //insert best edges
//    cout<<"Edges: "<<endl;
    for (int i = 0; i < svertices; ++i) {
      setI = dsu.findSet(i);
      //if edge exist
      if (bestEdge[setI].first != -1) {
        int j = bestEdge[setI].second;
        setJ = dsu.findSet(j);
        if (setI != setJ) {
          adjMatrix[i][j] = true;
          dsu.join(i, j);
          --nTrees;
          MSTCost += bestEdgeCost[setI];
 //         cout<<"("<<i<<", "<<j<<") ";
        }
        bestEdgeCost[setI] = bestEdgeCost[setJ] = DBL_MAX;
        bestEdge[setI] = bestEdge[setJ] = {-1, -1};
      }
    }
//    cout<<endl;
//    cout<<"Sets:"<<endl;
    set<int> sets;
    for (int i = 0; i < svertices; ++i) 
      sets.insert(dsu.findSet(i));
    /*
    for (int components : sets) {
      cout<<"\t"<<components<<": ";
      for (int i = 0; i < svertices; ++i) 
        if (dsu.findSet(i) == components)
          cout<<i<<", ";
      cout<<endl;
    }
    */
  }
  //greedy boruvka
//  cout<<"greedy borukva"<<endl;
  for (int i = 0; i < svertices; ++i) {
//    cout<<"\tIn vertex "<<i<<endl;
    nodeI = vertices[i];
    dsu.clean();
//    cout<<"\tremoving "<<i<<" from the tree"<<endl;
    //remove i from the graph   
    for (int j = 0; j < svertices; ++j) {
      nodeJ = vertices[j];
      if (adjMatrix[i][j] || adjMatrix[j][i]) {
        MSTCost -= gvrpReducedGraph[nodeI->id][nodeJ->id];
        adjMatrix[i][j] = adjMatrix[j][i] = false; 
      }
    }
//    cout<<"\tpopulating DSU"<<endl;
    //populate dsu
    for (int j = 0; j < svertices; ++j) 
      for (int k = 0; k < svertices; ++k) 
        if (adjMatrix[j][k])
            dsu.join(j, k);
//    cout<<"\tSets:"<<endl;
    set<int> sets;
    for (int i = 0; i < svertices; ++i) 
      sets.insert(dsu.findSet(i));
    /*
    for (int components : sets) {
      cout<<"\t\t"<<components<<": ";
      for (int i = 0; i < svertices; ++i) 
        if (dsu.findSet(i) == components)
          cout<<i<<", ";
      cout<<endl;
    }
    */
    //get number of components
//    cout<<"\tgetting number of components"<<endl;
    nTrees = sets.size();
//    cout<<"\tborukva again"<<endl;
    //borukva
    while (nTrees > 2) {
//      cout<<"\tnTrees: "<<nTrees<<endl;
      //calculate best cuts
      for (int j = 0; j < svertices; ++j) {
        if (j != i) {
          nodeJ = vertices[j];
          setJ = dsu.findSet(j);
          for (int k = 0; k < svertices; ++k) {
            if (k != i) {
              nodeK = vertices[k];
              cost = gvrpReducedGraph[nodeJ->id][nodeK->id];
              //cout<<"\t\t"<<j<<" "<<k<<" "<<nodeJ->id<<" "<<nodeK->id<<" "<<cost<<endl;
              if (j != k && setJ != dsu.findSet(k) && cost < bestEdgeCost[setJ]) {
                bestEdgeCost[setJ] = cost;
                bestEdge[setJ] = {j, k};
              }
            }
          }
        }
      }
      //insert best edges
//      cout<<"\tInserted edges: "<<endl<<"\t\t";
      for (int j = 0; j < svertices; ++j) {
        setJ = dsu.findSet(j);
        //if edge exist
        if (bestEdge[setJ].first != -1) {
          int k = bestEdge[setJ].second;
          setK = dsu.findSet(k);
          if (setJ != setK) {
            adjMatrix[j][k] = true;
            dsu.join(setJ, setK);
            --nTrees;
            MSTCost += bestEdgeCost[setJ];
//            cout<<"("<<k<<", "<<j<<"), ";
          }
          bestEdgeCost[setJ] = bestEdgeCost[setK] = DBL_MAX;
          bestEdge[setJ] = bestEdge[setK] = {-1, -1};
        }
      }
//      cout<<endl;
//      cout<<"\tSets:"<<endl;
      sets.clear();
      for (int i = 0; i < svertices; ++i) 
        sets.insert(dsu.findSet(i));
      /*
      for (int components : sets) {
        cout<<"\t\t"<<components<<": ";
        for (int i = 0; i < svertices; ++i) 
          if (dsu.findSet(i) == components)
            cout<<i<<", ";
        cout<<endl;
      }
      */
      bestLB = max(closests[i].first + closests[i].second + MSTCost, bestLB);
    }
  }
  return bestLB;
}
