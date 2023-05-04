import gurobipy as gp
import numpy as np

m = gp.Model()

budget = 500
# objective vector
C = np.array([170, 180, 120, 110, 120, 100])
# price vector
P = np.array([190, 240, 160, 120, 150, 140])

xs = m.addMVar(len(P), vtype=gp.GRB.BINARY, obj=C)

m.addConstr(C @ xs <= budget)
# 3 requires 2
m.addConstr(xs[2 - 1] >= xs[3 - 1])
# 5 and 6 not together
m.addConstr(1 - xs[6 - 1] >= xs[5 - 1])
m.addConstr(1 - xs[5 - 1] >= xs[6 - 1])

m.ModelSense = gp.GRB.MAXIMIZE
m.update()
m.optimize()

print(m.ObjVal)
print(xs)

