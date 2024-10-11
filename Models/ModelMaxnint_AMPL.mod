param I;
param J;
param t_max;
param Mr;
param maxnint{1..J};
param ds{1..I};
param ra{1..J};
param rho{1..J-1};
param minSplit;
param proc{1..I,1..J};
param nu;

var S{1..I, 1..J, 1..3} integer;
var E{1..I, 1..J, 1..3} integer;
var P{1..I, 1..J, 1..3} integer;
var Y{1..I, 1..J, 1..3} binary;
var Z{1..I, 1..J, 1..t_max} binary;
var X{1..I, 1..J, 1..4, 1..t_max} binary;
var N{1..I, 1..J, 1..t_max} integer;
var alpha{1..I, 1..J, 1..t_max} binary;

minimize TotalCost:
    sum{i in 1..I}(E[i,J,maxnint[J]]-ds[i]);

s.t. S_con{i in 1..I, j in 1..J, mu in 1..maxnint[j]}: S[i,j,mu] >= 0;
s.t. E_con{i in 1..I, j in 1..J, mu in 1..maxnint[j]}: E[i,j,mu] >= 0;
s.t. P_con{i in 1..I, j in 1..J, mu in 1..maxnint[j]}: P[i,j,mu] >= 0;
s.t. C1 {i in 1..I, j in 1..J, t in 1..t_max}: N[i,j,t] >= 0;
#s.t. C2 {i in 1..I} t_max*V[i]     >= F[i]-ds[i];
#s.t. C3 {i in 1..I} -t_max*(1-V[i]) <= F[i]-ds[i];
#s.t. C4 {i in 1..I, j in 1..J, mu in 1..maxnint[j]} E[i,j,mu] <= F[i];
s.t. C5 {i in 1..I, j in 1..J, mu in (1..maxnint[j]-1)}: E[i,j,mu] <= S[i,j,mu+1];
s.t. C6 {i in 1..I, j in 1..J, mu in 1..maxnint[j]}: P[i,j,mu] >= minSplit*Y[i,j,mu];
s.t. C7 {i in 1..I, j in 1..J, mu in 1..maxnint[j]-1}: Y[i,j,mu] >= Y[i,j,mu+1];
s.t. C8 {i in 1..I, j in 1..J, mu in 1..maxnint[j]}: P[i,j,mu] <= Mr*Y[i,j,mu];
s.t. C9 {i in 1..I, j in 1..J, mu in 1..maxnint[j]}: E[i,j,mu] = P[i,j,mu]+S[i,j,mu];
s.t. C10 {i in 1..I, j in 1..J, mu in 2..maxnint[j]}: E[i,j,mu-1] >= S[i,j,mu]-t_max*Y[i,j,mu];
s.t. C11 {i in 1..I, j in 1..J, t in 1..t_max, mu in 1..maxnint[j]}:S[i,j,mu] - t <= t_max*(1-X[i,j,mu,t]);
s.t. C12 {i in 1..I, j in 1..J, t in 1..t_max, mu in 1..maxnint[j]}:t - E[i,j,mu] + 1 <= t_max*(1-X[i,j,mu,t]);
s.t. C13 {i in 1..I, j in 1..J, mu in 1..maxnint[j]}:sum{t in 1..t_max}(X[i,j,mu,t]) == P[i,j,mu];
s.t. C14 {i in 1..I, j in 1..J, t in 1..t_max}:sum{mu in 1..maxnint[j]}(X[i,j,mu,t]) == Z[i,j,t];
s.t. C15 {t in 1..t_max, j in 1..J}:sum{i in 1..I}(N[i,j,t]) <= ra[j];
s.t. C16 {i in 1..I, j in 1..J}:sum{t in 1..t_max}(N[i,j,t]) == proc[i,j];
s.t. C17 {i in 1..I, j in 1..J, t in 1..t_max}:N[i,j,t] - Mr*Z[i,j,t] <= 0;
s.t. C18 {i in 1..I, j in 1..J, t in 1..t_max}:N[i,j,t] >= nu - (1 - sum{mu in 1..maxnint[j]}(X[i,j,mu,t]))*Mr;

s.t. C19 {i in 1..I, j in 1..2, t in 1..t_max}: sum{h in 1..t-1}(N[i,j,h]) >= rho[j]*proc[i,j]-alpha[i,j,t]*proc[i,j];
s.t. C20 {i in 1..I, j in 1..2, t in 1..t_max-1}: sum{h in 1..t+1}(N[i,j+1,h]) <= proc[i,j+1]*(1-alpha[i,j,t]);
s.t. C21 {i in 1..I, j in 3..J-1}: E[i,j,maxnint[j]] <= S[i,j+1,1];
s.t. C22 {i in 1..I, j in 1..J-1}: E[i,j,maxnint[j]] <= E[i,j+1,maxnint[j+1]];