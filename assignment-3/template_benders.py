import sys
import pathlib
import gurobipy as gp

###############################################################################################
##### HEADER ####
###############################################################################################

# Group Members:
#
# Student 1:
# Student 2:
# ...

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

def exercise1():
    pass

##############
# EXERCISE 3 #
##############

def exercise3():
    pass

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
    
    exercise1()

    exercise3()

    ###############################################
    # end of modifiable block
    ###############################################

if __name__ == "__main__":
    main()
