import pyscipopt as scip
# Create a new SCIP instance
solver = scip.Model("MyLP")
v = {}
s = {}
p = {}
e = {}
y = {}
x = {}
z = {}
n = {}
alpha = {}
for i in I:
  v[i] = solver.addVar(vtype="B", name="v(%s)" % i)
  for j in J:
    for mu in range(maxnint[j]):
      y[i,j,mu] = solver.addVar(vtype="B",      name="y(%s,%s,%s)" % (i,j,mu))
      s[i,j,mu] = solver.addVar(vtype="I",lb=1, name="s(%s,%s,%s)" % (i,j,mu))
      p[i,j,mu] = solver.addVar(vtype="I",lb=0, name="p(%s,%s,%s)" % (i,j,mu))
      e[i,j,mu] = solver.addVar(vtype="I",lb=1, name="e(%s,%s,%s)" % (i,j,mu))
      for t in H:
        x[i,j,mu,t] = solver.addVar(vtype="B", name="x(%s,%s,%s,%s)" % (i,j,mu,t))
    for t in H:
      z[i,j,t]     = solver.addVar(vtype="I",lb=0, name="z(%s,%s,%s)" % (i,j,t))
      n[i,j,t]     = solver.addVar(vtype="I",lb=0, name="n(%s,%s,%s)" % (i,j,t))
      alpha[i,j,t] = solver.addVar(vtype="B",  name="alpha(%s,%s,%s)" % (i,j,t))
      
# ------------ OBJECTIVE FUNCTION

solver.setObjective(scip.quicksum(e[(i,j,maxnint[j]-1)]- ds[i]*v[i] for i in I), "minimize")

# ------------ CONSTRAINTS

for i in I:
  solver.addCons(v[(i)] * Mt      >= e[(i,j,maxnint[j]-1)]-ds[i])
  solver.addCons(-(1-v[(i)]) * Mt <= e[(i,j,maxnint[j]-1)]-ds[i])

for i in I:
  for j in J:
    solver.addCons(scip.quicksum(n[(i,j,t)] for t in H) == proc[i,j])
    for t in H:
      solver.addCons(z[(i,j,t)]   == scip.quicksum(x[(i,j,mu,t)] for mu in range(maxnint[j])))
      solver.addCons(n[(i,j,t)]   <= Mr*z[(i,j,t)])
      solver.addCons(n[(i,j,t)]   >= 5-(1-z[(i,j,t)])*Mr)
    for mu in range(maxnint[j]):
      solver.addCons(p[(i, j, mu)]>= minSplit*y[(i, j, mu)])
      solver.addCons(p[(i, j, mu)]<= Mt*y[(i, j, mu)])
      solver.addCons(e[(i, j, mu)]== s[(i, j, mu)]+p[(i, j, mu)])
      solver.addCons(p[(i,j,mu)]  == scip.quicksum(x[(i,j,mu,t)] for t in H))
      for t in H:
        solver.addCons(s[(i, j, mu)]-t <= Mt*(1-x[(i,j,mu,t)]))
        solver.addCons(t-e[(i,j,mu)] <= Mt*(1-x[(i,j,mu,t)]))
    if maxnint[j] >1:
      for mu in range(maxnint[j]-1):
        solver.addCons(e[(i, j, mu)]<= s[(i, j, mu+1)] )
        solver.addCons(y[(i, j, mu)]>= y[(i, j, mu+1)] )
        solver.addCons(e[(i, j, mu)]>= s[(i, j, mu+1)]-Mt*y[(i, j, mu+1)] )


for t in H:
  for j in J:
    solver.addCons(scip.quicksum(n[(i,j,t)] for i in I)<=ra[j])


for i in I:
  for j in range(5):
    solver.addCons(e[(i,j,maxnint[j]-1)]+1<=e[(i,j+1,maxnint[j+1]-1)])
    for t in range(1,Mt-1):
      solver.addCons(scip.quicksum(n[(i,j,h)] for h in range(1,t+1))  >= rho[j]*proc[i,j]-alpha[i,j,t]*proc[i,j])
      solver.addCons(scip.quicksum(n[(i,j+1,h)] for h in range(1,t))<= (1-alpha[i,j,t])*proc[i,j+1])
      
for i in I:
  for j in range(2,6):
    for j_ in range(j-1):
      for t in H:
        solver.addCons(s[(i,j,0)]>=t*z[(i,j_,t)])

# ---------- PRINT CONSTRAINTS
#solver.printStatistics()
solver.setIntParam("parallel/maxnthreads", 2)  
# ----------
solver.optimize()

import pandas as pd
dates = pd.DataFrame(list(zip(
  [(i+1)  for i in I for j in J for mu in range(maxnint[j]) ],
  [(j+1) for i in I for j in J for mu in range(maxnint[j]) ],
  [(mu+1) for i in I for j in J for mu in range(maxnint[j]) ],
  [solver.getVal(s[i,j,mu]) for i in I for j in J for mu in range(maxnint[j]) ],
  [solver.getVal(p[i,j,mu]) for i in I for j in J for mu in range(maxnint[j]) ],
  [solver.getVal(e[i,j,mu]) for i in I for j in J for mu in range(maxnint[j]) ],
  [solver.getVal(y[i,j,mu]) for i in I for j in J for mu in range(maxnint[j]) ])),
                                              columns =['I','J','Mu','S','P','E','Y']
)
energy = pd.DataFrame(list(zip(
  [(i+1)  for i in I for j in J for t in H],
  [(j+1) for i in I for j in J   for t in H ],
  [ t for i in I for j in J   for t in H],
  [solver.getVal(n[i,j,t]) for i in I for j in J for t in H ],
  [solver.getVal(z[i,j,t]) for i in I for j in J for t in H ],
  [solver.getVal(alpha[i,j,t]) for i in I for j in J for t in H ]
  )),
                                              columns =['I','J','week','N','Z',"alpha"]
)

X = pd.DataFrame(list(zip(
  [(i+1)  for i in I for j in J for mu in range(maxnint[j])  for t in H],
  [(j+1) for i in I for j in J for mu in range(maxnint[j])  for t in H],
  [(mu+1) for i in I for j in J for mu in range(maxnint[j])  for t in H],
  [ t  for i in I for j in J for mu in range(maxnint[j]) for t in H],
  [solver.getVal(x[i,j,mu,t]) for i in I for j in J for mu in range(maxnint[j]) for t in H],
  [solver.getVal(z[i,j,t]) for i in I for j in J for mu in range(maxnint[j])  for t in H])),
                                              columns =['Project','job','Mu','week','X','Z']
)  




