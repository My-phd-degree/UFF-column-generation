Solution format:

#_of_routes

For each route: 
#_of_customers #_of_stations 
node sequence


Notes:

- the first node of each route is node 0 (start depot), the last node of each route is node n+1 (end depot) 

- #_of_customers includes the start depot but not the end depot, i.e., the node sequence contains a total of (#_of_customers + #_of_stations + 1) nodes

- nodes with negative id in the customer sequence are refueling stations, i.e., node -i corresponds to the i-th refueling station