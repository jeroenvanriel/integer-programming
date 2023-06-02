import sys
import gurobipy as gp
import pathlib
import numpy as np

###############################################################################################
##### HEADER ####
###############################################################################################

# Group Members:
#
# Student 1: Jeroen van Riel
# Student 2: Koen Moors

######################################
# AUXILIARY FUNCTIONS
######################################

def read_instance(filename):
    '''
    reads a graph

    filename - path to file containing the graph
    '''

    f = open(filename, 'r')

    edges = []
    for line in f:
        splitline = line.split()
        (u,v) = (int(splitline[0]), int(splitline[1]))
        edges.append((u-1,v-1))

    f.close()

    return edges

class Graph:

    def __init__(self, edges, x):
        self.n = max(max(e) for e in edges)
        self.m = len(edges)
        self.V = np.arange(2*self.n)
        self.E = np.zeros(shape=(2*self.n, 2*self.n))
        for (e0, e1) in edges:
            w = 1 - x[e0] - x[e1]
            self.E[e0,e1+self.n] = w
            self.E[e0+self.n,e1] = w
            self.E[e1,e0+self.n] = w
            self.E[e1+self.n,e0] = w
 
    def printSolution(self, dist):
        print("Vertex \t Distance from Source")
        for node in self.V:
            print(node, "\t\t", dist[node])
 
    # A utility function to find the vertex with
    # minimum distance value, from the set of vertices
    # not yet included in shortest path tree
    def minDistance(self, dist, sptSet):
 
        # Initialize minimum distance for next node
        min = 1e7
 
        # Search not nearest vertex not in the
        # shortest path tree
        for v in self.V:
            if dist[v] < min and sptSet[v] == False:
                min = dist[v]
                min_index = v
 
        return min_index
 
    # Function that implements Dijkstra's single source
    # shortest path algorithm for a graph represented
    # using adjacency matrix representation
    def dijkstra(self, src):
 
        dist = [1e7] * self.n * 2
        dist[src] = 0
        sptSet = [False] * self.n * 2
 
        for cout in self.V:
 
            # Pick the minimum distance vertex from
            # the set of vertices not yet processed.
            # u is always equal to src in first iteration
            u = self.minDistance(dist, sptSet)
 
            # Put the minimum distance vertex in the
            # shortest path tree
            sptSet[u] = True
 
            # Update dist value of the adjacent vertices
            # of the picked vertex only if the current
            # distance is greater than new distance and
            # the vertex in not in the shortest path tree
            for v in self.V:
                if (self.graph[u][v] > 0 and
                   sptSet[v] == False and
                   dist[v] > dist[u] + self.graph[u][v]):
                    dist[v] = dist[u] + self.graph[u][v]

        self.printSolution(dist)
        return dist

    def parallel_distances(self):
        dist = np.zeros(self.n)
        for v in self.V:
            dist[v] = self.dijkstra(v)[v+self.n]
        return dist
    
    def find_cycle(self):
        p_dist = self.parallel_distances()
        return np.where(p_dist == np.min(p_dist))[0]



##############################
# EXERCISE 5
##############################

def exercise5(edges):

    m = gp.Model("stableSet")

    m.Params.Cuts = 0
    m.Params.Presolve = 0
    m.Params.PreCrush = 1

    ###########################
    Graph(edges, x=1)
    ###########################

##############################
# MAIN FUNCTION
##############################

def main():

    path = None
    if len(sys.argv) > 1:
        path = pathlib.Path(sys.argv[1])
    else:
        path = pathlib.Path("graph1.txt")
        print(f"WARNING: Instance not specified. Defaulting to {path}.")
    assert path.exists(), f"Path {path} does not exist. ({path.resolve()})"
    print(f"Reading instance {path}")

    edges = read_instance(path)

    exercise5(edges)

if __name__ == "__main__":
    main()
