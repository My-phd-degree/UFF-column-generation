# UFF-column-generation

## To make changes in the docker bapdock image
docker run -it -v /home/MYUSERFOLDERNAME/vrp_solver/GVRP:/GVRP --name GVRP bapdock
docker commit DOCKER-HASH-CODE bapdock

## To define the args array

args = ["GVRP/data/EMH/AS1_20c3sU10.txt", "--preprocessings", "GVRP/data/EMH/preprocessings/AS1_20c3sU10.txt_preprocessed_edges.txt", "--instance_type", "EMH"]
