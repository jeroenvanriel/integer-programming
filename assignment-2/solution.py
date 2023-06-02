import sys
import gurobipy as gp
import pathlib
import numpy as np
import networkx as nx

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

def edit_path(path, n):
    """Given a path, it finds the shortest cycle in it. Moreover, reduces all elements modulo n.

    Examples:
    edit_path([0,1,2,3,0], 10)
    >>> [0,1,2,3]
    edit_path([0,1,2,3,1,0], 10)
    >>> [1,2,3]
    edit_path([0,11,2,13,1,10], 10)
    >>> [1,2,3]
    """
    path = [x % n for x in path]
    while path[0] == path[-1]:
        x = path[0]
        path = path[1:-1]
    return [x] + path


class MyGraph:

    def __init__(self, edges):
        self.n = max(max(e) for e in edges) + 1
        self.m = len(edges)
        self.edges = edges
        self.V = np.arange(2*self.n)
        self.E = np.zeros(shape=(2*self.n, 2*self.n))

    def update_weights(self, x):
        for i in range(2*self.n):
            for j in range(2*self.n):
                self.E[i,j] = np.infty
        for (e0, e1) in self.edges:
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
        for v in range(self.n*2):
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
 
        for cout in range(self.n*2):
 
            # Pick the minimum distance vertex from
            # the set of vertices not yet processed.
            # u is always equal to src in first iteration
            try:
                u = self.minDistance(dist, sptSet)
            except UnboundLocalError:
                continue
            # Put the minimum distance vertex in the
            # shortest path tree
            sptSet[u] = True
 
            # Update dist value of the adjacent vertices
            # of the picked vertex only if the current
            # distance is greater than new distance and
            # the vertex in not in the shortest path tree
            for v in range(self.n*2):
                if (self.E[u][v] > 0 and
                   sptSet[v] == False and
                   dist[v] > dist[u] + self.E[u][v]):
                    dist[v] = dist[u] + self.E[u][v]

        # self.printSolution(dist)
        return dist

    def parallel_distances(self):
        dist = np.zeros(self.n)
        for v in range(self.n):
            dist[v] = self.dijkstra(v)[v+self.n]
        return dist
    
    def find_cycle(self):
        p_dist = self.parallel_distances()
        m_dist = np.where(p_dist == np.min(p_dist))[0]
        assert len(m_dist) % 2 == 1, "Found cycle is even"
        return m_dist



##############################
# EXERCISE 5
##############################

def exercise5(edges):

    m = gp.Model("stableSet")

    m.Params.Cuts = 0
    m.Params.Presolve = 0
    m.Params.PreCrush = 1

    ###########################
    n = max(max(e) for e in edges) + 1

    x = m.addVars(range(n), vtype=gp.GRB.BINARY, name="x")
    for (i,j) in edges:
        m.addConstr(x[i] + x[j] <= 1)

    m.setObjective(
        gp.quicksum(1 * x[i] for i in range(n)),
        sense=gp.GRB.MAXIMIZE
    )

    G = MyGraph(edges)
    
    def callback(model, where):
        if where == gp.GRB.Callback.POLLING:
            # Polling callback: This is called every so-many milliseconds.
            return

        if where == gp.GRB.Callback.MIPNODE:
            # MIPNODE callback: A branch-and-bound tree node is solved to optimality
            # Get status
            status = model.cbGet(gp.GRB.Callback.MIPNODE_STATUS)
            if status != gp.GRB.Status.OPTIMAL:
                return

            # Get (fractional) solution
            x_frac = model.cbGetNodeRel(x)

            print()
            print(f"# Solved branch-and-bound tree node to optimality. Search violated cover inequality:")

            # Find a violated cover inequality
            G: nx.Graph = nx.Graph()
            for (i,j) in edges:
                w = 1 - x_frac[i] - x_frac[j]
                w = min(w , 1)
                w = max(w, 0)
                G.add_edge(u_of_edge=i, v_of_edge=j+n, weight=w)
                G.add_edge(u_of_edge=i+n, v_of_edge=j, weight=w)

            min_l = np.infty
            for i in range(n):
                l = nx.dijkstra_path_length(G, source=i, target=i+n, weight="weight")
                if l < min_l:
                    min_l = l
                    min_i = i

            path = nx.dijkstra_path(G, source=min_i, target=min_i+n, weight="weight")
            path = edit_path(path, n)   # first and last element are the same

            cut = gp.quicksum(x[i] for i in path) <= (len(path) - 1) / 2
            if cut is None:
                print("Search is unsuccessful")
                return
            else:
                print(cut)
                model.cbCut(cut)
            return

        if where == gp.GRB.Callback.MIPSOL:
            x_sol = model.cbGetSolution(x)
            # get the solution value; see https://www.gurobi.com/documentation/9.5/refman/cb_codes.html#sec:CallbackCodes
            x_sol_val = model.cbGet(gp.GRB.Callback.MIPSOL_OBJ)


        pass

    m.optimize(callback=callback)
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
