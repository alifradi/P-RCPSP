/*********************************************
 * OPL 12.6.0.0 Model
 * Author: Ali
 * Creation Date: Feb 6, 2023 at 5:25:07 AM
 *********************************************/
 //DATA INPUT
int I = ...;
int J = ...;
int T =...;
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
dvar int     e[1..I][1..J];
dvar int     p[1..I][1..J];
dvar int     s_[1..I][1..J];
dvar boolean d[1..I];
dvar float   x[1..I][1..J][1..T];
dvar float   c[1..I][1..J][1..T];
dvar boolean s[1..I][1..J][1..T];
dvar boolean f[1..I][1..J][1..T];
dvar float   n[1..I][1..J][1..T];
dvar boolean Ind[1..I][1..J][1..T];
//MAIN PROGRAM

//OBJECTIVE FUNCTION
minimize sum(i in 1..I)(2*e[i,J]-ds[i]*d[i]);
// CONSTRAINTS
subject to{
forall(j in 1..J,t in 1..T)                  (sum(i in 1..I)(x[i,j,t]*proc[i,j])<=ra[j]); // respect budget

forall(i in 1..I, j in 1..J,t in 1..T)       (x[i,j,t]<=1);          // define ratio UB 
forall(i in 1..I, j in 1..J,t in 1..T)       (x[i,j,t]>=0);          // define ratio LB
forall(i in 1..I, j in 1..J)                 (c[i,j,T]==1);          // define ratio sum to 1

forall(i in 1..I, j in 1..J)(sum(t in 1..T)  (s[i,j,t])==1);         // start  j only once
forall(i in 1..I, j in 1..J)(sum(t in 1..T)  (f[i,j,t])==1);         // finish j only once

forall(i in 1..I, j in 1..J,t in 1..T)       (x[i,j,t]<=sum(h in 1..t)(s[i,j,h])); // forcing indication of starting the activity at the first ratio >0
forall(i in 1..I, j in 1..J,t in 1..T)       (f[i,j,t]<=c[i,j,t]);                 // preventing end before satisfiying  progress ratio


forall(i in 1..I, j in 1..J,t in 1..T)       (e[i,j]==sum(t in 1..T)(t*f[i,j,t])+1); // calculate possible end of activity j
forall(i in 1..I, j in 1..J-1,t in 1..T)     (s[i,j+1,t]<=1+c[i,j,t]-rho[j]);// prevent j+1 to begin while rho[j]% of j is incomplete
forall(i in 1..I, j in 1..J,t in 1..T)       (sum(h in 1..t)(f[i,j,h]) >= 0.1 - sum(h in t+1..T)(x[i,j,h])); // end of j if no strict positif ratio to allocate 
forall(i in 1..I, j in 1..J,t in 1..T)       (sum(h in 1..t)(x[i,j,h]) >= 0.1 * sum(h in 1..t)(s[i,j,h]));

forall(i in 1..I, j in 1..J-1,t in 1..T)     (c[i,j+1,t]<=c[i,j,t]);         // respect order finishing tasks: task j finishs before task j+1 ends

//forall(i in 1..I)(d[i]>=e[i,J]-ds[i]);                                       // d_i = max(0,e_iJ-ds_i)
//forall(i in 1..I)(d[i]>=0);                                                  // d_i = max(0,e_iJ-ds_i)

forall(i in 1..I)( T*d[i]     >=e[i,J]-ds[i]);
forall(i in 1..I)(-T*(1-d[i]) <=e[i,J]-ds[i]);
  
  
//forall(i in 1..I, j in 1..J)                 (sum(t in 1..T)(Ind[i,j,t])<=p[i,j]/2);  // biggest bound of suspension number
forall(i in 1..I, j in 1..J, t in 1..T)      (Ind[i,j,t]>=x[i,j,t]);                  // activate indicator on allocated effort
forall(i in 1..I, j in 1..J, t in 1..T)      (Ind[i,j,t]<=100*x[i,j,t]);              // turn off mask if no allocated work
forall(i in 1..I, j in 1..J, t in 1..T)      (sum(h in 1..t)(Ind[i,j,h])<=T*sum(h in 1..t)(s[i,j,h]));// 
forall(i in 1..I, j in 1..J, t in 1..T)      (sum(h in t..T)(Ind[i,j,h])<=T*sum(h in t..T)(f[i,j,h]));// 
forall(i in 1..I, j in 1..J)                 (p[i,j]-sum(t in 1..T)(Ind[i,j,t])<=maxnint[j]);// 

// REPORTING
forall(i in 1..I,j in 1..J,t in 1..T)        (n[i,j,t] == x[i,j,t]*proc[i,j]);         // report used energy in each period per time unit
forall(i in 1..I, j in 1..J)                 (s_[i,j]  == sum(t in 1..T)(t*s[i,j,t])); // report starting time 
forall(i in 1..I, j in 1..J,t in 1..T)       (c[i,j,t] == sum(h in 1..t)(x[i,j,h]));   // compute progress in ij until time t 
forall(i in 1..I, j in 1..J,t in 1..T)       (p[i,j]   == e[i,j]-s_[i,j]);             // compute time used to finish tasks

//x[1,1,1]==0;
//x[1,1,2]==0;
//x[1,1,4]==0;
//x[1,1,3]==0.1;
//c[1,1,6]==1;
}


execute {
  var f = new IloOplOutputFile(I+"p_"+T+"w"+"_minE"+maxnint[1]+maxnint[2]+maxnint[3]+maxnint[4]+maxnint[5]+maxnint[6]+"splits.csv");
  f.writeln("design,  activityId,  week, S, P, N, dueDate, capacity_for_week");
  for(var i=1; i<=I; i++){
    for(var j=1; j<=J; j++){
        for(var t=1; t<=T; t++){
        f.writeln(design[i]+", "+j+", "+t+", "+s_[i][j]+", "+p[i][j]+", "+n[i][j][t]+", "+ds[i]+", "+ ra[j]);
      }       
           
    }
  }
  f.close();
}


