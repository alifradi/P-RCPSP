include "globals.mzn";
int: I;
int: J;
int: Mr;
int: minSplit;
% Planning horizon
int: t_max ;
array [1..J-1] of float: rho;
float: nu;
array [1..I, 1..J] of int: proc;
array [1..J] of int: ra;
array [1..J] of int: maxnint;
array [1..I] of int: ds;






array [1..I] of var bool:                        V;
array [1..I] of var int:                         F;
array [1..I, 1..J, 1..6] of var int:              S;
array [1..I, 1..J, 1..6] of var int:              E;
array [1..I, 1..J, 1..6] of var int:              P;
array [1..I, 1..J, 1..6] of var bool:             Y;
array [1..I, 1..J, 1..t_max] of var bool:         Z;
array [1..I, 1..J, 1..6, 1..t_max] of var bool:   X;
array [1..I, 1..J,       1..t_max] of var int:  N;
array [1..I, 1..J, 1..t_max] of var bool: alpha;
%array [1..I, 1..J, 1..t_max] of var int:      Cons;
%array [1..I, 1..J, 1..t_max] of var bool: gamma;
%array [1..I, 1..J, 1..t_max] of var bool: beta;

% Domain of variables

constraint forall(i in 1..I, j in 1..J, mu in 1..maxnint[j])(S[i,j,mu]>=0);
constraint forall(i in 1..I, j in 1..J, mu in 1..maxnint[j])(P[i,j,mu]>=0);
constraint forall(i in 1..I, j in 1..J, mu in 1..maxnint[j])(E[i,j,mu]>=0);
constraint forall(i in 1..I, j in 1..J, mu in 1..maxnint[j])(E[i,j,mu]<=F[i]);
constraint forall(i in 1..I)(F[i]>=0);
constraint forall(i in 1..I, j in 1..J,k in 1..J, t in 1..t_max)(N[i,j,t]>=0);
var 1..1000: obj;





%constraint obj = sum(t in 1..t_max, j in 1..J)(ra[j]-sum(i in 1..I)(N[i,j,t]));
constraint obj = sum(i in 1..I)(F[i]-ds[i]*V[i]);
%solve maximize obj;
solve minimize obj;

% Unconstrained multi project scheduling allowing preemption

constraint forall(i in 1..I)( t_max*V[i]     >=F[i]-ds[i]);
constraint forall(i in 1..I)(-t_max*(1-V[i]) <=F[i]-ds[i]);
constraint forall(i in 1..I, j in 1..J, mu in (1..maxnint[j]-1))(E[i,j,mu]<=S[i,j,mu+1]);
constraint forall(i in 1..I, j in 1..J, mu in 1..maxnint[j])(P[i,j,mu]>=minSplit*Y[i,j,mu]);
constraint forall(i in 1..I, j in 1..J, mu in 1..maxnint[j]-1)(Y[i,j,mu]>=Y[i,j,mu+1]);
constraint forall(i in 1..I, j in 1..J, mu in 1..maxnint[j])(P[i,j,mu]<=Mr*Y[i,j,mu]);
constraint forall(i in 1..I, j in 1..J, mu in 1..maxnint[j])(E[i,j,mu]==P[i,j,mu]+S[i,j,mu]);
constraint forall(i in 1..I, j in 1..J, mu in 2..maxnint[j])(E[i,j,mu-1]>=S[i,j,mu]-t_max*Y[i,j,mu]);

% Resources constraints

constraint forall(i in 1..I, j in 1..J, t in 1..t_max, mu in 1..maxnint[j])(S[i, j, mu] - t <= t_max*(1-X[i, j, mu, t]));
constraint forall(i in 1..I, j in 1..J, t in 1..t_max, mu in 1..maxnint[j])(t - E[i, j, mu] + 1 <= t_max*(1-X[i, j, mu, t]));
constraint forall(i in 1..I, j in 1..J, mu in 1..maxnint[j])(sum(t in 1..t_max)(X[i, j, mu, t])==P[i, j, mu]); 
constraint forall(i in 1..I, j in 1..J, t in 1..t_max)(sum(mu in 1..maxnint[j])(X[i, j, mu, t])==Z[i, j, t]); 
constraint forall(t in 1..t_max, j in 1..J)(sum(i in 1..I)(N[i,j,t])<=ra[j]);
constraint forall(i in 1..I, j in 1..J)(sum(t in 1..t_max)(N[i,j,t]) == proc[i,j]);
constraint forall(i in 1..I, j in 1..J,  t in 1..t_max)(N[i,j,t]-Mr*Z[i, j, t]<=0);
constraint forall(i in 1..I, j in 1..J,  t in 1..t_max)(N[i,j,t]>=5-(1-Z[i, j, t])*Mr);


%% Resources dates
% overlapping conditions

constraint forall(i in 1..I, j in 1..J-1,t in 1..t_max)(sum( h in 1..t)(N[i,j,h]) >= rho[j]*proc[i,j]-alpha[i,j,t]*proc[i,j]);
constraint forall(i in 1..I, j in 1..J-1,t in 1..t_max-1)(sum( h in 1..t+1)(N[i,j+1,h]) <= proc[i,j+1]*(1-alpha[i,j,t]));

constraint forall(i in 1..I, j in 1..J-1)(E[i,j,maxnint[j]]<=E[i,j+1,maxnint[j+1]]);
%constraint forall(i in 1..I, j in 1..J,t in 1..t_max)(Cons[i,j,t] >= proc[i,j]-beta[i,j,t]*proc[i,j]);
%constraint forall(i in 1..I, j in 1..J-1,t in 1..t_max-1)(Cons[i,j+1,t+1] + nu <= proc[i,j+1]+Mr*(1-beta[i,j,t]));
constraint forall(i in 1..I,t in 1..t_max,j in 3..6, j_ in 1..j-2)(S[i,j,1]>=t*Z[i,j_,t]);



output
[
"Solution \n"
]
++
[
  "p: "++show(i)++" j: "++show(j)++" [S-(P)-E]: [" ++ show(S[i,j,mu])  ++ "-(" ++ show(P[i,j,mu]) ++ ")-" ++ show(E[i,j,mu]) ++ "] " ++ 
 " F["++show(i)++"]= "++show(F[i])++ " Y["++show(i)++", "++show(j)++", "++show(mu)++"]= "++show(bool2int(Y[i,j,mu]))
 ++" \n"
  | i in 1..I, j in 1..J, mu in 1..maxnint[j] where(fix(P[i,j,mu]>0))
]

++
[
"Allocated work hours/activity/time \n"
]++
[  
"N["++show(i)++","++show(j)++", "++show(t)++"]= " ++show(N[i,j,t]) ++ "\n "
| i in 1..I, j in 1..J,t in 1..t_max where(fix(N[i,j,t]>0))
]
%++["Consumed resources \n"]++
%["pjt: "++show(i)++" job: "++show(j)++" d:"++show(t)++"   "++ show(Cons[i,j,t])++" threshold is: "++show(rho[j]*proc[i,j])++ " sucessor : "++show(Cons[i,j+1,t])++ " alpha= " ++ show(bool2int(alpha[i,j,t]))++ "\n" | i in 1..I, j in 1..J-1, t in 1..t_max]
++
["\n Resources allocation: R [R1, R2, R3, R4, R5, R6] "]++


[ "\n week: "++ show(t)++ 
" R: ["++show(sum(i in 1..I)(N[i,1,t]))++ ", "++
show(sum(i in 1..I)(N[i,1,t]))++ ", "++
show(sum(i in 1..I)(N[i,2,t]))++ ", "++
show(sum(i in 1..I)(N[i,3,t]))++ ", "++
show(sum(i in 1..I)(N[i,4,t]))++ ", "++
show(sum(i in 1..I)(N[i,5,t]))++ ", "++
show(sum(i in 1..I)(N[i,6,t]))++ 
"]"++ "\n "
| t in 1..t_max 
]++[
"Objective value: \n"++show(obj)
]
;
