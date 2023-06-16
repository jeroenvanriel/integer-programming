import sys
import pathlib
import gurobipy as gp
import numpy as np

###############################################################################################
##### HEADER ####
###############################################################################################

# Group Members:
#
# Student 1: Jeroen van Riel
# Student 2: Koen Moors

###############################################################################################
### MAIN PART ###
###############################################################################################

#######################
# AUXILIARY FUNCTIONS #
#######################

def readFLInstance(
    path                # specification of the file to be read
):
    file = open(path, 'r')
    lines = file.readlines()
    file.close()

    # read number of facilities and customers
    line = lines[0].split()

    nfacilities = int(line[0])
    ncustomers = int(line[1])
    facilities = [i for i in range(nfacilities)]
    customers = [j for j in range(ncustomers)]

    # read capacities and fixed costs
    capacity = []
    fixcost = []
    for i in range(1, nfacilities + 1):
        line = lines[i].split()
        capacity.append(float(line[0]))
        fixcost.append(float(line[1]))

    # read customers' data
    demand = []
    cost = []

    for i in customers:
        idx = nfacilities + 1 + i
        line = lines[idx].split()
        demand.append(float(line[0]))
        customercost = [float(k) for k in line[1:(nfacilities + 1)]]
        cost.append(customercost)

    return nfacilities, ncustomers, facilities, customers, capacity, fixcost, cost, demand

def readClassifierInstance(filename):
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

##############
# EXERCISE 1 #
##############

def exercise1(n, m, f, c):
    # note: c is a list of lists

    # this may be used for manual verification
    # m = 3
    # n = 2
    # f = f[0:2]
    # c = c[0:m*n]

    customers = list(range(m))
    facilities = list(range(n))

    print(f"{n} facilities")
    print(f"{m} customers")

    master = gp.Model()
    master.ModelSense = gp.GRB.MINIMIZE

    # first- and second-stage variables
    x = master.addMVar(n, vtype='B', obj=f)
    # second-stage objective variable
    z = master.addVar(obj=1)

    # ensure second-stage problem is feasible
    master.addConstr(gp.quicksum(x) >= 1)

    # TODO (maybe more efficient): use sparse matrices
    B = np.zeros((m + m*n, n))
    D = np.zeros((m + m*n, m*n + m + m*n))
    d = np.zeros((m + m*n))

    cy = np.zeros(m*n + m + m*n)

    for i in customers:
        for j in facilities:
            # "all demand satisfied" constraint
            D[i, i*n + j] = 1
            # slack variables
            D[i, m*n + i] = -1
            # right-hand side
            d[i] = 1

            # "only open facilities can serve" constraint
            D[m + i*n + j, i*n + j] = 1
            # slack variables
            D[m + i*n + j, m*n + m + i*n + j] = 1
            # first-stage variables
            B[m + i*n + j, j] = -1 
        
            # delivery cost
            cy[i*n + j] = c[i][j]

    # construct subproblems
    subfea = gp.Model()
    subfea.ModelSense = gp.GRB.MAXIMIZE
    w = subfea.addMVar(m + m*n, lb=float('-inf'))
    subfea.addConstr(w.transpose() @ D <= np.zeros_like(cy).transpose())

    subopt = gp.Model()
    subopt.ModelSense = gp.GRB.MAXIMIZE
    p = subopt.addMVar(m + m*n, lb=float('-inf'))
    subopt.addConstr(p.transpose() @ D <= cy.transpose())
    
    opt = False
    while not opt:
        opt = True

        master.update()
        master.optimize()

        print(f"--- solution = {x.X}, with objective = {master.ObjVal}")

        # get best solution
        x_opt = x.X
        z_opt = z.X

        subfea.setObjective(w.transpose() @ (d - B @ x_opt))
        # check whether feasibility cut is violated
        subfea.update()
        subfea.optimize()
        if subfea.ObjVal > 0:
            # violated, add constraint
            wj = w.X
            master.addConstr(-(B.transpose() @ wj) @ x <= -np.inner(wj, d))
            opt = False
            print("--- VIOLATED feasibility")
        else:
            # p.Obj = d - B @ x_opt
            subopt.setObjective(p.transpose() @ (d - B @ x_opt))
            # check whether optimality cut is violated
            subopt.update()
            subopt.optimize()
            if subopt.ObjVal > z_opt:
                # violated, add constraint
                pi = p.X
                master.addConstr(-(pi.transpose() @ B) @ x - z <= - np.inner(pi, d))
                opt = False
                print(f"--- VIOLATED optimality")
    

    # solve the second-stage
    second = gp.Model()
    second.ModelSense = gp.GRB.MINIMIZE
    # ignore the slack variables
    y = second.addMVar(m*n, obj=cy[0:m*n])

    # fullfillment constraints
    for i in customers:
        second.addConstr(gp.quicksum(y[i*n : (i+1)*n]) >= 1)
        # open constraints
        for j in facilities:
            second.addConstr(y[i*n + j] <= x_opt[j])
    
    second.update()
    second.optimize()

    opening_cost = x_opt @ f
    print(f"objective (verification) = {second.ObjVal + opening_cost}")
    print(f"solution y = {y.X}")

##############
# EXERCISE 3 #
##############

def exercise3(K, X, Y):
    n = len(X[0])
    print(f'K = {K}')
    print(f'{len(X)} x-points')
    print(f'{len(Y)} y-points')
    print(f'dimension of points = {n}')

    eps = 0.0001

    master = gp.Model()
    master.ModelSense = gp.GRB.MINIMIZE

    # first-stage variables
    u = master.addMVar(K, vtype='B', obj=1)
    s = master.addMVar((len(Y), K), vtype='B')

    # every y-point is cut off    
    master.addConstr(s @ np.ones(K) >= 1)

    # every y-point is only cut off by active inequalities
    for y in range(len(Y)):
        for i in range(K):
            master.addConstr(s[y, i] <= u[i])

    # construct subproblems
    B = n * max([abs(x[i]) for i in range(n) for x in X])
    M = B + eps + n * max([abs(y[i]) for i in range(n) for y in Y])

    subproblem = gp.Model()
    subproblem.ModelSense = gp.GRB.MINIMIZE
    subproblem.Params.InfUnbdInfo = 1
    constraints = []

    a = subproblem.addMVar(n, lb=-1, ub=1)
    b = subproblem.addVar(lb=-B, ub=B)

    for x in X:
        subproblem.addConstr(a @ np.array(x) <= b)
    for iy, y in enumerate(Y):
        constraints.append(
            subproblem.addConstr(a @ np.array(y) - b >= 0) # RHS is placeholder
        )

    done = False
    while not done:
        master.update()
        master.optimize()

        # get primal solution
        s_opt = s.X
        # print(f'--- SOLUTION: found solution s = {s_opt}')

        # check the subproblems
        done = True
        for i in range(K):
            # update the subproblem
            for iy in range(len(Y)):
                constraints[iy].RHS = eps - M * (1 - s_opt[iy][i])

            # solve it
            subproblem.update()
            subproblem.optimize()

            # add cut if not feasible
            if subproblem.status == gp.GRB.INFEASIBLE:
                done = False

                dual = subproblem.FarkasDual
                # print('--- DUAL')
                # print(dual)

                # We only need the dual coefficients corresponding to the last inequalities,
                # i.e., the ones corresponding to the y-points (hence 'len(X) + iy')
                # N.B. Because these inequalities are >=, the RHS needs to be negated.
                master.addConstr(gp.quicksum(
                    dual[len(X) + iy] * (-eps + M * (1 - s[iy][i])) for iy in range(len(Y))
                ) <= 0)

                break # to immediately continue solving the master problem once we found a cut

    print(f'--- DONE: solution objective {master.ObjVal}')

#################
# MAIN FUNCTION #
#################

def main():

    pathFL = None
    pathClassifier = None
    if len(sys.argv) > 2:
        pathFL = pathlib.Path(sys.argv[1])
        pathClassifier = pathlib.Path(sys.argv[2])
    else:
        pathFL = pathlib.Path("benders_test.dat")
        pathClassifier = pathlib.Path("data_cg_cross4.cg")
        print(f"WARNING: Instance not specified. Defaulting to {pathFL} and {pathClassifier}.")
    assert pathFL.exists(), f"Path {pathFL} does not exist. ({pathFL.resolve()})"
    assert pathClassifier.exists(), f"Path {pathClassifier} does not exist. ({pathClassifier.resolve()})"
    print(f"Reading instance {pathFL} and {pathClassifier}")

    nfactories, ncustomers, factories, customers, capacity, fixcost, cost, demand = readFLInstance(pathFL)
    K, X, Y = readClassifierInstance(pathClassifier)

    ###############################################
    # only modify this function within this block
    ###############################################

    # add your input to the following functions
    
    exercise1(nfactories, ncustomers, fixcost, cost)

    exercise3(K, X, Y)

    ###############################################
    # end of modifiable block
    ###############################################

if __name__ == "__main__":
    main()
