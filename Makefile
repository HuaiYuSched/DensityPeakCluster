CC = g++
CXXFLAGS = -g

all: clean cluster

cluster: cluster.cpp
	$(CC) cluster.cpp -o cluster

srun: all
	./cluster
run: cluster
	./cluster
	gnuplot plotclust.p

plot:
	gnuplot plotclust.p

clean:
	rm -f *.o cluster Clusterfile Decision_graph

mrproper:
	rm -f *.o cluster Clusterfile Decision_graph *.jpeg
