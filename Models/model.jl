#LIST OF VARIABLES
@variable(m, s[i in I,j in J,mu in 1:maxnint[j]] >=0, Int)       # Starting time of each slice mu of job j, project i
@variable(m, p[i in I,j in J,mu in 1:maxnint[j]] >=0, Int)       # Processing time of each slice mu of job j, project i 
@variable(m, e[i in I,j in J,mu in 1:maxnint[j]] >=0, Int)       # End of sub-job mu of j of project i 
@variable(m, f[I] >= 0,                               Int)       # finish date of projects
@variable(m,        n[i in I, j in J, t in H] >=0)               # Number of hours allocated from resource k to process slice mu of job j, project i 
#@variable(m, c[i in I, j in J, t in H])                          # consommation resources

@variable(m,     z[i in I, j in J, t in H],                Int)  # Binary variable for processing job j of i in t
@variable(m, v[I] >= 0,                                    Bin)  # variable max(0,f-ds)
@variable(m, x[i in I, j in J,mu in 1:maxnint[j], t in H], Bin)  # Binary variable for processing subjob mu of (j , i) in t
@variable(m, y[i in I, j in J,mu in 1:maxnint[j]],         Bin)  # Binary auxilary variable 
@variable(m, alpha[i in I, j in J, t in H],                Bin)
#@variable(m,  beta[i in I, j in J, t in H],                Bin)
#@variable(m, gamma[i in I, j in J, t in H],                Bin)  
#OBJECTIVE FUNCTION
@objective(m, Min, sum( f[i]-ds[i]*v[i] for i in I) )           # minimize number of delayed items and time periods of delay 
                                    
#LIST OF CONSTRAINTS

# Project Delay
@constraint(m, numbdelay1[i in I],                                 Mt*v[i]     >=  f[i]-ds[i])  
@constraint(m, numbdelay2[i in I],                                 -Mt*(1-v[i]) <=  f[i]-ds[i])  
@constraint(m, finishProj[i in I, j in J, mu in 1:maxnint[j]  ],   f[i]>=e[i,j,mu])
# Spans of possible resource allocation
@constraint(m, slice_prec[i in I, j in J, mu in 1:maxnint[j]-1],   e[i,j,mu]<= s[i, j, mu+1])                                             #  slice precedence                                                 #  max(0,f-ds)
@constraint(m, slice_Inc1[i in I, j in J, mu in 1:maxnint[j]  ],   p[i, j, mu] >= minSplit*y[i, j, mu])                                   #  slices incremental order 1
@constraint(m, slice_Inc2[i in I, j in J, mu in 1:maxnint[j]-1],   y[i, j, mu]>=y[i, j, mu+1])                                            #  slices incremental order 2
@constraint(m, slice_Inc3[i in I, j in J, mu in 1:maxnint[j]  ],   p[i, j, mu]<=Mt*y[i, j, mu])                                           #  slices incremental order 3
@constraint(m, end__slice[i in I, j in J, mu in 1:maxnint[j]  ],   e[i, j, mu]==p[i, j, mu]+s[i, j, mu])                                  #  define end of slices
@constraint(m, dummyTails[i in I, j in J, mu in 2:maxnint[j]  ],   e[i,j,mu-1]>=s[i,j,mu]-Mt*y[i,j,mu])
# Defining mask of resources allocation
@constraint(m, procwind1[i in I,j in J,t in H,mu in 1:maxnint[j]], s[i,j,mu]-t <= Mt*(1-x[i,j,mu,t]))                                     #  no resources to allocate before s 
@constraint(m, procwind2[i in I,j in J,t in H,mu in 1:maxnint[j]], t-e[i,j,mu]+ 1 <= Mt*(1-x[i,j,mu,t]))                                 #  no resources to allocate after e
@constraint(m, satisfyslice[i in I,j in J,mu in 1:maxnint[j]],     sum(x[i,j,mu,t] for t in H)==p[i,j,mu])                                #  call for resources to process slice mu of job j, i
@constraint(m, satisfy__job[i in I,j in J,t in H],                 sum(x[i,j,mu,t] for mu in 1:maxnint[j])==z[i,j,t])                  #  call for resources to process slice mu of job j, i
# Constraint resources  to mask and limits
@constraint(m, reslimits[t in H],                           sum(n[i,j,t] for i in I)<=ra[j])                       #  respect resources limits
@constraint(m, allocress[i in I,j in J],                           sum(n[i,j,t] for t in H) == proc[i][j])                   #  allocate enough resources 
@constraint(m, holdresources[i in I, j in J, t in H],              n[i,j,t]-Mr*z[i,j,t]<=0)                                             #  not allowing resources if there is no job being processed

# Overlapping conditions 
@constraint(m, minWork2_Pass[i in I,j in 1:5, t in H],             sum(n[i,j,  h] for h in 1:t)>=rho[j]*proc[i][j]-alpha[i,j,t]*proc[i][j])
@constraint(m, BlockPassBmin[i in I,j in 1:5, t in H],             sum(n[i,j+1,h] for h in 1:t)<=       proc[i][j+1]*(1-alpha[i,j,t]))
#Precedence conditions
@constraint(m, precedenc[i in I, j in 1:5],                        e[i,j,maxnint[j]]<=e[i,j+1,maxnint[j+1]])
#@constraint(m, [i in I,j in J,   t in H],                          sum(n[i,j  ,h] for h in 1:t)>=proc[i][j]-beta[i,j,t]*proc[i][j])
#@constraint(m, [i in I,j in 1:5, t in H],                          sum(n[i,j+1,h] for h in 1:t)<=proc[i][j+1])