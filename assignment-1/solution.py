import sys
import gurobipy as gp
import pathlib

import matplotlib.pyplot as plt
import numpy as np

import time

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

def keyboard_terminate(model, where):
    """Terminate m.optimize() using `ctrl + c`."""
    try:
        pass
    except KeyboardInterrupt:
        model.terminate()

def plot_lines(u, a, b):
    """Sample points from lines for the first two variables,
    such that we can completely inspect problems with n=2.
    Note that this may be also be used for n>2 problems,
    by visualizing "cross sections" of the space."""
    step = 0.5
    xs = np.arange(-10, 10, step)
    for i in range(len(u)):
        if u[i].X < 0.5:
            # skip inactive (u_i == 0) lines
            continue
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

    m.ModelSense = gp.GRB.MINIMIZE
    m.update()
    m.optimize()

    plot_lines(u, a, b)

######################################
# EXERCISE 3
######################################

def singleton_cols(primal, Y):
    """Take all singletons {y} for all y in Y as columns."""
    cols = np.eye(len(Y), dtype=np.int16).tolist()
    vars = primal.addVars(len(Y), obj=1, vtype='C', name='vars')
    cons = primal.addConstrs((vars[j] >= 1 for j in range(len(Y))), name='c1')

    return cols, cons

def partition_cols(primal, pricing, Y):
    """Greedily construct columns by solving the problem from exercise 2."""
    m = pricing.copy()

    for j in range(len(Y)):
        s = m.getVarByName(f"s_{Y[j]}")
        s.Obj = 1
    
    cols = []
    removed = 0
    # while there are constraints for Y
    while removed < len(Y):
        m.update()
        m.optimize()

        col = [0 for j in range(len(Y))]
        for j in range(len(Y)):
            v = m.getVarByName(f"s_{Y[j]}")
            if v is not None:
                col[j] = v.X

        cols.append(col)
        # remove the y's from the problem
        for j in range(len(Y)):
            if col[j]:
                c = m.getConstrByName(str(j))
                m.remove(c)
                s = m.getVarByName(f"s_{Y[j]}")
                m.remove(s)
                removed += 1

    vars = primal.addVars(len(cols), obj=1, vtype='C', name='vars')
    cons = primal.addConstrs((gp.quicksum([vars[i] * cols[i][j] for i in range(len(cols))]) >= 1
                       for j in range(len(Y))), name='c1')

    return cols, cons

def exercise3(X, Y, eps):
    # dimension of the points
    N = len(X[0])

    # tolerance for pricing problem
    tolerance = 0.0001

    # construct primal problem with feasible basis
    # we could start by including every I_y := \{ y \}, for all y \in Y
    primal = gp.Model()
    primal.params.OutputFlag = 0  # for readability, output is turned off
    primal.ModelSense = gp.GRB.MINIMIZE

    # construct the pricing problem
    pricing = gp.Model()
    pricing.params.OutputFlag = 0 # for readability, output is turned off
    pricing.ModelSense = gp.GRB.MAXIMIZE

    # for every y, we have a variable indicating whether it is in the set
    s = { j : pricing.addVar(obj=0, vtype=gp.GRB.BINARY, name=f"s_{Y[j]}") for j in range(len(Y)) }

    # the feasible region encodes all the linear separable sets
    # so we need orientation and offset
    a = { n: pricing.addVar(float("-inf"), float("inf"), 
                            obj=0, name=f"a_{n}") for n in range(N) }
    b = pricing.addVar(obj=0, name="b")

    M = 10 # "large" number

    # separation constraints
    for x in X:
        pricing.addConstr(gp.quicksum([a[n] * x[n] for n in range(N)]) <= b)
    for j in range(len(Y)):
        pricing.addConstr(M * (1 - s[j]) + gp.quicksum([a[n] * Y[j][n] for n in range(N)]) >= b + eps, name=str(j))

    pricing.update()

    # initial variables and constraints
    cols, cons = partition_cols(primal, pricing, Y)
    # cols, cons = singleton_cols(primal, Y)

    pricing.Params.BestObjStop = 1 + tolerance

    while True:
        primal.update()
        primal.optimize(keyboard_terminate)

        # obtain the dual solution to update the pricing
        # problems objective (constr.pi is shadow price)
        for j in range(len(Y)):
            s[j].Obj = cons[j].pi

        # solve the pricing problem
        pricing.update()
        pricing.optimize(keyboard_terminate)

        # reduced costs non-negative? => optimal
        if pricing.ObjVal <= 1 + tolerance:
            print("optimal solution found!")
            break

        for k in range(pricing.SolCount):
            pricing.Params.SolutionNumber = k
            if pricing.PoolObjVal >= 1:
                # column s is incidence vector of y \in Y in pattern I
                col = [int(s[j].Xn) for j in range(len(Y))]
                # now add all y that violate the boundary ax<=b
                for j in range(len(Y)):
                    if sum(a[n].Xn * Y[j][n] for n in range(N)) >= b.Xn +eps:
                        col[j] = 1
                cols.append(col)
                # add new variable to the primal problem
                primal.addVar(obj=1, vtype='C', column=gp.Column(col, cons.values()))

    print(f"total number of columns: {len(cols)}")
    print(f"objective value: {primal.ObjVal}")
    return cols

##############################
# EXERCISE 4
##############################

def exercise4(cols):
    R = len(cols[0])  # len(Y)
    L = len(cols)

    m = gp.Model()
    m.ModelSense = gp.GRB.MINIMIZE

    vars = m.addVars(L, obj=1, vtype=gp.GRB.BINARY, name='vars')
    cons = m.addConstrs((gp.quicksum(cols[i][j]*vars[i] for i in range(L)) >= 1 for j in range(R)))

    m.update()
    m.optimize()
    print(f"objective value: {m.ObjVal}")

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
    t = time.time()
    exercise1(K, X, Y, eps)
    print(f"Time it took to complete exercise 1: {time.time()-t} s", end= '\n \n')
    t = time.time()
    cols = exercise3(X, Y, eps)
    print(f"Time it took to complete exercise 3: {time.time()-t} s", end= '\n \n')
    t = time.time()
    exercise4(cols)
    print(f"Time it took to complete exercise 4: {time.time()-t} s", end= '\n \n')
    
    ###############################################
    # end of modifiable block
    ###############################################

if __name__ == "__main__":
    main()
