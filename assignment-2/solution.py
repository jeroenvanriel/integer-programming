import sys
import gurobipy as gp
import pathlib

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

    def __init__(self, vertices, edges, x):
        pass

##############################
# EXERCISE 5
##############################

def exercise5(edges):

    m = gp.Model("stableSet")

    m.Params.Cuts = 0
    m.Params.Presolve = 0
    m.Params.PreCrush = 1

    ###########################
    print(edges)
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
