import sys
import gurobipy as gp
import pathlib

import matplotlib.pyplot as plt
import numpy as np

###############################################################################################
##### HEADER ####
###############################################################################################

# Group Members:
#
# Student 1: Jeroen van Riel
# Student 2: Koen Moors
# ...

######################################
# AUXILIARY FUNCTIONS
######################################

def read_instance(filename):
    '''
    reads an instance of the polyhedral classifier problem

    filename - path to file containing instance
    '''

    f = open(filename, 'r')

    K = -1
    X = []
    Y = []

    for line in f:
        if line.startswith('k'):
            K = int(line.split()[1])
        elif line.startswith('x'):
            point = [int(i) for i in line.split()[1:]]
            X.append(point)
        elif line.startswith('y'):
            point = [int(i) for i in line.split()[1:]]
            Y.append(point)

    f.close()

    return K, X, Y

######################################
# EXERCISE 1
######################################

def exercise1(K, X, Y, eps):
    # this should be large w.r.t. to the input instance
    M = 10

    # dimension of the points
    N = len(X[0])

    plot_data(X, Y)

    m = gp.Model()

    # We use i to index constraints (lines), j to index points from Y
    # and 0 <= n < N to index coordinates.
    
    # constraint active / inactive decision variables
    u = [m.addVar(obj=1, vtype=gp.GRB.BINARY, name=f"u_{i}") for i in range(K)]
    # constraint orientations
    a = { (i, n):
            m.addVar(-float('inf'), float('inf'), obj=0, vtype=gp.GRB.CONTINUOUS, name=f"a_{i}_{n}")
            for n in range(N)
            for i in range(K) }
    # constraint offsets
    b = [m.addVar(obj=0, vtype=gp.GRB.CONTINUOUS, name=f"b_{i}") for i in range(K)]
    # constraint violation decision variables
    s = { (i, j):
            m.addVar(obj=0, vtype=gp.GRB.BINARY, name=f"s_{i}_{j}")
            for i in range(K)
            for j in range(len(Y)) }

    # constraint satisfaction for points in X
    for j in range(len(X)):
        for i in range(K):
            m.addConstr(gp.quicksum(a[i,n] * X[j][n] for n in range(N)) <= b[i])

    # constraint satisfaction for points in Y
    # provided that s[i,j] = 0 (no violation)
    for j in range(len(Y)):
        for i in range(K):
            m.addConstr(M * (s[i,j] + 1 - u[i]) + gp.quicksum(a[i,n] * Y[j][n] for n in range(N)) >= b[i] + eps)

    # each point in Y satisfies at least one constraint
    for j in range(len(Y)):
        m.addConstr(gp.quicksum(s[i,j] for i in range(K)) <= gp.quicksum(u[i] for i in range(K)) - 1)

    ### Bound a_i's away from zero vector ###
    # N.B.: these are not necessary for correctness, but they seem to speed up the
    # optimization for all the provided instances. Especially for the 'data_cg_cross4.cg'
    # instance, it reduces computation time from 96.59s to 10.35s.
    # # sign variables for a_i
    # p = { (i, n): m.addVar(obj=0, vtype=gp.GRB.BINARY, name=f"p_{i}_{n}") for i in range(K) for n in range(N) }
    # # sign constraints for a_i
    # for i in range(K):
    #     for n in range(N):
    #         m.addConstr(a[i,n] + p[i,n] * M >= 1)
    #         m.addConstr(a[i,n] - (1 - p[i,n]) * M <= -1)

    m.ModelSense = gp.GRB.MINIMIZE
    m.update()
    m.optimize()

    # Sample points from lines for the first two variables,
    # such that we can completely inspect problems with n=2.
    # Note that this may be also be used for n>2 problems,
    # by visualizing "cross sections" of the space.
    step = 0.5
    xs = np.arange(-10, 10, step)
    for i in range(K):
        if u[i].X < 0.5:
            # skip inactive (u_i == 0) lines
            continue
        if abs(a[i, 1].X) < eps / 2:
            # horizontal line
            ys = [0 for x in xs]
        else:
            # plot the first two coordinates
            ys = [(b[i].X - a[i,0].X * x) / a[i,1].X for x in xs]
        plt.plot(xs, ys, linestyle=('dotted' if u[i].X < 0.5 else 'solid'))
        plt.xlim([-10, 10])
        plt.ylim([-10, 10])
    plt.show()


def plot_data(X, Y):
    xs = np.array(X)
    ys = np.array(Y)
    plt.scatter(xs[:,0], xs[:,1], marker='x')
    plt.scatter(ys[:,0], ys[:,1], marker='o')

######################################
# EXERCISE 3
######################################

def exercise3(
        # provide input here
        ):
    pass


##############################
# EXERCISE 4
##############################

def exercise4(
        # provide input here
        ):
    pass


##############################
# MAIN FUNCTION
##############################

def main():

    path = None
    if len(sys.argv) > 1:
        path = pathlib.Path(sys.argv[1])
    else:
        path = pathlib.Path("data_cg_cube.cg")
        print(f"WARNING: Instance not specified. Defaulting to {path}.")
    assert path.exists(), f"Path {path} does not exist. ({path.resolve()})"
    
    print(f"Reading instance {path}")

    K, X, Y = read_instance(path)
    eps = 0.001

    ###############################################
    # only modify this function within this block
    ###############################################

    # add your input to the following functions

    exercise1(K, X, Y, eps)

    exercise3()

    exercise4()

    ###############################################
    # end of modifiable block
    ###############################################

if __name__ == "__main__":
    main()
