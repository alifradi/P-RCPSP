

df_gantt = DataFrame(s       = [value.(s[i,j,mu]) for i in I for j in J[i] for mu in 1:maxnint[i][j]],
p       = [value.(p[i,j,mu]) for i in I for j in J[i] for mu in 1:maxnint[i][j]],
e       = [value.(e[i,j,mu]) for i in I for j in J[i] for mu in 1:maxnint[i][j]],
pjt     = [i for i in I for j in J[i] for mu in 1:maxnint[i][j]],
act     = [j for i in I for j in J[i] for mu in 1:maxnint[i][j]],
mu      = [mu for i in I for j in J[i] for mu in 1:maxnint[i][j]],
endact  = [value.(e[i,j,maxnint[i][j]]) for i in I for j in J[i] for mu in 1:maxnint[i][j]]) ;




df_ressource = DataFrame(R1=[sum(value.(n[i,j,1,t]) for i in I for j in J[i]) for t in H],
    R2=[sum(value.(n[i,j,2,t]) for i in I for j in J[i]) for t in H],
    R3=[sum(value.(n[i,j,3,t]) for i in I for j in J[i]) for t in H],
    R4=[sum(value.(n[i,j,4,t]) for i in I for j in J[i]) for t in H],
    R5=[sum(value.(n[i,j,5,t]) for i in I for j in J[i]) for t in H],
    R6=[sum(value.(n[i,j,6,t]) for i in I for j in J[i]) for t in H],
    R7=[sum(value.(n[i,j,7,t]) for i in I for j in J[i]) for t in H]) ;





df_ressource_job = DataFrame(R1=[sum(value.(n[i,j,1,t])) for i in I for j in J[i] for t in H],
       R2=[sum(value.(n[i,j,2,t])) for i in I for j in J[i] for t in H],
       R3=[sum(value.(n[i,j,3,t])) for i in I for j in J[i] for t in H],
       R4=[sum(value.(n[i,j,4,t])) for i in I for j in J[i] for t in H],
       R5=[sum(value.(n[i,j,5,t])) for i in I for j in J[i] for t in H],
       R6=[sum(value.(n[i,j,6,t])) for i in I for j in J[i] for t in H],
       R7=[sum(value.(n[i,j,7,t])) for i in I for j in J[i] for t in H],
       task=[(i,j) for i in I for j in J[i] for t in H],
       time=[t for i in I for j in J[i] for t in H]) ;

