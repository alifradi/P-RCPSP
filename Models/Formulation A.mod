/*********************************************
 * OPL 22.1.0.0 Model
 * Author: alifr
 * Creation Date: 14 déc. 2023 at 12:00:53
 *********************************************/

 
  //DATA INPUT
int I = ...;
int J = ...;
int t_max =...;
float Mr=...;
int maxnint[1..J]=...;
string design[1..I]=...;
//int dl[1..I]=...;
int ds[1..I]=...;
float ra[1..J]=...;
float rho[1..J-1]=...;
int minSplit = ...;
float proc[1..I][1..J]=...;
//int q[1..J][1..J] = ...;
float nu =...;
//DECISION VARIABLES
//dvar boolean V[1..I];
//dvar int F[1..I];
dvar int T[1..I];
dvar int S[1..I][1..J][1..10];
dvar int E[1..I][1..J][1..10];
dvar int P[1..I][1..J][1..10];
dvar boolean Y[1..I][1..J][1..10];
dvar float Z[1..I][1..J][1..t_max];
dvar boolean X[1..I][1..J][1..10][1..t_max];
dvar float N[1..I][1..J][1..t_max];
dvar boolean alpha[1..I][1..J][1..t_max];
//dvar boolean beta[1..I][1..J][1..t_max];
//dvar boolean gamma[1..I][1..J][1..t_max];
//MAIN PROGRAM

//OBJECTIVE FUNCTION
//minimize sum(i in 1..I)(2*E[i,J,maxnint[J]] -ds[i]*V[i]);
//minimize sum(i in 1..I)(E[i,J,maxnint[J]] +100*V[i]);
minimize sum(i in 1..I)(E[i,J,maxnint[J]] +T[i]);
// CONSTRAINTS
subject to{
forall(i in 1..I, j in 1..J, mu in 1..maxnint[j]) S[i,j,mu]>=0;
forall(i in 1..I, j in 1..J, mu in 1..maxnint[j]) P[i,j,mu]>=0;
forall(i in 1..I, j in 1..J, mu in 1..maxnint[j]) E[i,j,mu]>=0;
forall(i in 1..I, j in 1..J,k in 1..J, t in 1..t_max) N[i,j,t]>=0;
//forall(i in 1..I)( t_max*V[i]     >=E[i,J,maxnint[J]]-ds[i]);
//forall(i in 1..I)(-t_max*(1-V[i]) <=E[i,J,maxnint[J]]-ds[i]);
forall(i in 1..I)(T[i]      >= E[i,J,maxnint[J]]-ds[i]);
forall(i in 1..I)(T[i]      >= 0);
forall(i in 1..I, j in 1..J, mu in (1..maxnint[j]-1)) E[i,j,mu]<=S[i,j,mu+1];
forall(i in 1..I, j in 1..J, mu in 1..maxnint[j]) P[i,j,mu]>=minSplit*Y[i,j,mu];
forall(i in 1..I, j in 1..J, mu in 1..maxnint[j]-1) Y[i,j,mu]>=Y[i,j,mu+1];
forall(i in 1..I, j in 1..J, mu in 1..maxnint[j]) P[i,j,mu]<=Mr*Y[i,j,mu];
forall(i in 1..I, j in 1..J, mu in 1..maxnint[j]) E[i,j,mu]==P[i,j,mu]+S[i,j,mu];
forall(i in 1..I, j in 1..J, mu in 2..maxnint[j]) E[i,j,mu-1]>=S[i,j,mu]-t_max*Y[i,j,mu];
forall(i in 1..I, j in 1..J, t in 1..t_max, mu in 1..maxnint[j]) S[i, j, mu] - t <= t_max*(1-X[i, j, mu, t]);
forall(i in 1..I, j in 1..J, t in 1..t_max, mu in 1..maxnint[j]) t - E[i, j, mu] + 1 <= t_max*(1-X[i, j, mu, t]);
forall(i in 1..I, j in 1..J, mu in 1..maxnint[j])sum(t in 1..t_max)(X[i, j, mu, t]) ==P[i, j, mu]; 
forall(i in 1..I, j in 1..J, t in 1..t_max) sum(mu in 1..maxnint[j])(X[i, j, mu, t]) ==Z[i, j, t]; 
forall(t in 1..t_max, j in 1..J) sum(i in 1..I)(N[i,j,t]) <=ra[j];
forall(i in 1..I, j in 1..J) sum(t in 1..t_max)(N[i,j,t]) >= proc[i,j];
forall(i in 1..I, j in 1..J) sum(t in 1..t_max)(N[i,j,t]) <= proc[i,j]+1;
forall(i in 1..I, j in 1..J, t in 1..t_max)N[i,j,t]-proc[i,j]*Z[i, j, t]<=0;
forall(i in 1..I, j in 1..J,  t in 1..t_max)(N[i,j,t]>=0.99-(1-Z[i, j, t])*2*proc[i,j]);
forall(i in 1..I, j in {1,2,5},t in 1..t_max) sum(h in 1..t-1)(N[i,j,h]) >= rho[j]*proc[i,j]-alpha[i,j,t]*proc[i,j];
forall(i in 1..I, j in {1,2,5},t in 1..t_max) sum(h in 1..t)(N[i,j+1,h]) <= proc[i,j+1]*(1-alpha[i,j,t]);
forall(i in 1..I, j in 1..J-1) E[i,j,maxnint[j]]+1<=E[i,j+1,maxnint[j+1]];
forall(i in 1..I, j in {3,4,5}) E[i,j,maxnint[j]]<=S[i,j+1,1];
//forall(i in 1..I, j in 3..J-1) E[i,j,maxnint[j]]<=S[i,j+1,1];
//forall(i in 1..I,t in 1..t_max,j in 3..6, j_ in 1..j-2)(S[i,j,1]>=t*Z[i,j_,t]);
//forall(i in 1..I, j in {1,2,3},t in 1..t_max)(N[i,j,t]<=80);

}



execute {
  var f = new IloOplOutputFile("A_EpT_lenSplit"+minSplit+"_"+I+"p_"+t_max+"_w"+nu+"_minE"+maxnint[1]+maxnint[2]+maxnint[6]+"splits.csv");
  f.writeln("design,  activityId, mu, week, S, P, N, dueDate, capacity_for_week");
  for(var p=1; p<=I; p++){
    for(var j=1; j<=J; j++){
      for(var mu=1; mu<=maxnint[j]; mu++){
        for(var t=1; t<=t_max; t++){
        f.writeln(design[p]+", "+j+", "+mu+", "+t+", "+S[p][j][mu]+", "+P[p][j][mu]+", "+N[p][j][t]+", "+ds[p]+", "+ ra[j]);
      }        
      }
    }
  }
  f.close();
}
 