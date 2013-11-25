#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "IL.h"
#include "HP.h"
#include "convert_Vienna.h"
#include "backtrack.h"
#include "RNAlocopt.h"

PRIVATE int backtrack_init();
PRIVATE int backtrack_ML(int i, int j, int p, int q);
PRIVATE int backtrack_M1(int i, int j, int p, int q);
PRIVATE int backtrack_M1E(int i, int j, int p, int q);
PRIVATE int backtrack_MLC(int i, int j);
PRIVATE int backtrack_EL(int j, int q);
PRIVATE int backtrack_M1q(int i, int j, int q);
PRIVATE int backtrack_M1qE(int i, int j, int q);
PRIVATE int backtrack_M1p(int i, int j, int q);
PRIVATE int backtrack_M1pE(int i, int j, int q);
PRIVATE int backtrack_M1long(int i, int j);
PRIVATE int backtrack_M1longE(int i, int j);
PRIVATE int backtrack_Zstar(int i, int j);
PRIVATE int backtrack_ZstarE(int i, int j);
PRIVATE int backtrack_IL(int i, int j, int k);

PRIVATE char structure[1000];
PRIVATE double energy;
PRIVATE int debug_flag;

PUBLIC char *perform_backtrack(int debug_flag_in){
  int i;

  debug_flag=debug_flag_in;
  strcpy(structure, "x");
  strcat(structure, sequence);
  for (i=1;i<=n;i++){
    structure[i]='.';
  }
  energy=0.0;
  backtrack_init();
  //  printf("structure is now %s\n", structure);
  //  printf("energy is %f.\n", energy);
  return structure;
}

PRIVATE int backtrack_init(){
  double sum, part, ran;
  int q;

  sum=0.;
  for (q=0;q<=Tails+1;q++){
    sum+=Z_EL_1[n][q];
  }
  
  ran=(double) rand()/RAND_MAX*sum;
  part=0.;
  for (q=0;q<=Tails+1;q++){
    part+=Z_EL_1[n][q];
    if (ran<=part){
      backtrack_EL(n,q);
      return 0;
    }
  }
  return 1;
}

PRIVATE backtrack_ML(int i, int j, int p, int q){
  int k, a, minvar;
  double ran, part;
  
  if (debug_flag){
    printf("ML, i,j,p,q;%d,%d,%d,%d.\n",i,j,p,q);
  }
  ran=(double) rand()/RAND_MAX*Z_ML[i][j][p][q];
  part=0.;
  for (k=i+4; k<=j-4; k++){
    minvar=MIN(Tails+1, k-i+1);
    for (a=0;a<=minvar; a++){
      part+=Z_ML[i][k][p][a]*Z_M1_val(k+1-a,j,a,q);
      if (ran<part){
	backtrack_ML(i,k,p,a);backtrack_M1(k+1-a,j,a,q);return 0;
      } 
      part+=Z_M1_val(i,k,p,a)*Z_M1_val(k+1-a,j,a,q);
      if (ran<part){
	backtrack_M1(i,k,p,a);
	backtrack_M1(k+1-a,j,a,q);
	return 0;
      } 
      if (k+1-Tails>i){
	part+=Z_ML[i][k][p][Tails+1]*Z_M1_val(k+1-Tails,j,Tails,q);
	if (ran<part){
	  backtrack_ML(i,k,p,Tails+1);
	  backtrack_M1(k+1-Tails,j,Tails,q);return 0;
	}
	part+=Z_M1_val(i,k,p,Tails+1)*Z_M1_val(k+1-Tails,j,Tails,q);
	if (ran<part){
	  backtrack_M1(i,k,p,Tails+1);
	  backtrack_M1(k+1-Tails,j,Tails,q);return 0;
	}
      }
    }
  }
  return 1;
}

PRIVATE backtrack_M1(int i, int j, int p, int q){
 
  if (debug_flag){
    printf("M1, i,j,p,q;%d,%d,%d,%d.\n",i,j,p,q);
  }
  if (p==Tails+1)
    if (q==Tails+1)
      backtrack_M1long(i,j);
    else
      backtrack_M1q(i,j,q);
  else if (q==Tails+1)
    backtrack_M1p(i,j,p);
  else{
    energy+=ML_base*(p+q);
    backtrack_Zstar(i+p,j-q);
  }
  if (debug_flag){
    printf("got to end of M1");
  }

  return 0;
}

PRIVATE backtrack_M1E(int i, int j, int p, int q){
 
  if (p==Tails+1)
    if (q==Tails+1)
      backtrack_M1longE(i,j);
    else
      backtrack_M1qE(i,j,q);
  else if (q==Tails+1)
    backtrack_M1pE(i,j,p);
  else
    backtrack_ZstarE(i+p,j-q);
  return 0;
}


PRIVATE backtrack_MLC(int i, int j){
  int p,q,prange,qrange;
  double ran,part;
  
  if (debug_flag){
    printf("MLC, i, j, are %d, %d.\n", i,j);
  }

  part=0.;
  ran=(double) rand()/RAND_MAX*Z_MLC[i][j];
  prange=MIN(Tails+1,j-i-5);
  for (p=0;p<=prange;p++){
    qrange=MIN(Tails+1,j-i-5-p);
    for (q=0;q<=qrange;q++){
      if (delta_min(i,j,p,q)){
	part+=Z_ML[i+1][j-1][p][q]*exp(-(ML_close+ML_Energy(j,i,S0))/kT);
	if (ran<part){
	  energy+=ML_close+ML_Energy(j,i,S0);
	  backtrack_ML(i+1,j-1,p,q);
	  return 0;
	}
      }
    }
  }
  return 1;
}

PRIVATE backtrack_EL(int j, int q){
  int minvar, a, k;
  double ran, part;

  if (debug_flag){
    printf("EL, j, q, are %d, %d.\n", j,q);
  }
  ran=(double) rand()/RAND_MAX*Z_EL_1[j][q];
  part=0.;
  for (k=0;k<=j-4;k++){
    minvar=MIN(Tails, k);
    for (a=0;a<=minvar;a++){
      part+=Z_EL_1[k][a]*Z_M1_E_val(k+1-a,j,a,q);
      if (ran<part){
	backtrack_EL(k,a);
	backtrack_M1E(k+1-a,j,a,q);
	return 0;
      }
    }
    if (k+1-Tails>=1){
      part+=Z_EL_1[k][Tails+1]*Z_M1_E_val(k+1-Tails,j,Tails,q);
      if (ran<part){
	backtrack_EL(k,Tails+1);
	backtrack_M1E(k+1-Tails,j,Tails,q);
	return 0;
      }
    }
  }
  return 0;//corresponds to the empty structure, up to j.
}

PRIVATE backtrack_M1q(int i, int j, int q){
  double ran;

  ran=(double) rand()/RAND_MAX*Z_M1q[i][j][q];
  if (ran<exp(-ML_base/kT)*Z_M1q[i+1][j][q]){
    energy+=ML_base;
    backtrack_M1q(i+1,j,q);
  }
  else
    backtrack_Zstar(i+Tails+1,j-q);
  return 0;
}

PRIVATE backtrack_M1qE(int i, int j, int q){
  double ran;

  ran=(double) rand()/RAND_MAX*Z_M1q_E[i][j][q];
  if (ran<Z_M1q_E[i+1][j][q])
    backtrack_M1qE(i+1,j,q);
  else
    backtrack_ZstarE(i+Tails+1,j-q);
  return 0;
}

PRIVATE backtrack_M1p(int i, int j, int p){
  double ran;

  ran=(double) rand()/RAND_MAX*Z_M1p[i][j][p];
  if (ran<exp(-ML_base/kT)*Z_M1p[i][j-1][p]){
    energy+=ML_base;
    backtrack_M1p(i,j-1,p);
  }
  else
    backtrack_Zstar(i+p,j-Tails-1);
  return 0;
}

PRIVATE backtrack_M1pE(int i, int j, int p){
  double ran;

  ran=(double) rand()/RAND_MAX*Z_M1p_E[i][j][p];
  if (ran<Z_M1p_E[i][j-1][p])
    backtrack_M1pE(i,j-1,p);
  else
    backtrack_ZstarE(i+p,j-Tails-1);
  return 0;
}

PRIVATE backtrack_M1long(int i, int j){
  double ran;

  ran=(double) rand()/RAND_MAX*Z_M1long[i][j];
  if (ran<exp(-ML_base/kT)*Z_M1long[i][j]){
    energy+=ML_base;
    backtrack_M1long(i,j-1);
  }
  else
    backtrack_M1q(i,j-1,Tails);
  return 0;
}

PRIVATE backtrack_M1longE(int i, int j){
  double ran;

  ran=(double) rand()/RAND_MAX*Z_M1long_E[i][j];
  if (ran<Z_M1long_E[i][j])
    backtrack_M1longE(i,j-1);
  else
    backtrack_M1qE(i,j-1,Tails);
  return 0;
}

PRIVATE backtrack_Zstar(int i, int j){
  int start, end, k, ip, jp;
  double part, ran;

  if (debug_flag){
    printf("Zstar, i, j, are %d, %d.\n", i,j);
  }
  part=0.;
  ran=(double) rand()/RAND_MAX*Zstar[i][j];
  start=IL_start(i,j);
  end=IL_end(i,j);

  for (k=start;k<=end;k++){
    ip=IL[k].ip;jp=IL[k].jp;
    if (ML_Energy(i,j,S0)+IL[k].energy<ML_Energy(ip,jp,S0)){
      part+=IL[k].part*exp(-ML_Energy(i,j,S0)/kT);
      if (ran<part){
	energy+=ML_Energy(i,j,S0);
	structure[i]='('; structure[j]=')';
	structure[ip]='('; structure[jp]=')';
	if (debug_flag){
	  printf("Zstar, i, j, ip, jp are %d, %d, %d, %d.\n", i, j, ip, jp);
	  int jdeb;
	  printf("<<%s>>\n", structure);
	  int lb=0, rb=0;
	  for (jdeb=0;jdeb<=n;jdeb++){
	    if (structure[jdeb]=='(') lb++;
	    if (structure[jdeb]==')') rb++;
	  }
	  printf("%d, %d.\n", lb, rb);
	}
	backtrack_IL(i,j,k);
	return 0;
      }
    }
  }
  return 1;
}

PRIVATE backtrack_ZstarE(int i, int j){
  int start, end, k, ip, jp;
  double part, ran;

  if (debug_flag){
    printf("ZstarE, i, j, are %d, %d.\n", i,j);
  }
  part=0.;
  ran=(double) rand()/RAND_MAX*ZEstar[i][j];
  start=IL_start(i,j);
  end=IL_end(i,j);

  for (k=start;k<=end;k++){
    ip=IL[k].ip;jp=IL[k].jp;
    if (EL_Energy(i,j,S0)+IL[k].energy<EL_Energy(ip,jp,S0)){
      part+=IL[k].part*exp(-EL_Energy(i,j,S0)/kT);
      if (ran<part){
	energy+=EL_Energy(i,j,S0);
	if (debug_flag){
	  printf("ZstarE, i, j, ip, jp are %d, %d, %d, %d.\n", i, j, ip, jp);
	  int jdeb;
	  printf("<<%s>>\n", structure);
	  int lb=0, rb=0;
	  for (jdeb=0;jdeb<=n;jdeb++){
	    if (structure[jdeb]=='(') lb++;
	    if (structure[jdeb]==')') rb++;
	  }
	  printf("%d, %d.\n", lb, rb);
	}
	structure[i]='('; structure[j]=')';
	structure[ip]='('; structure[jp]=')';
	backtrack_IL(i,j,k);
	return 0;
      }
    }
  }
  return 1;
}

PRIVATE backtrack_IL(int i, int j, int k){
  int ip, jp, l, ipp, jpp;
  double ran, part, E0;

  if (debug_flag){
    printf("IL, i, j, are %d, %d.\n", i,j);
  }
  ran=(double) rand()/RAND_MAX*IL[k].part;
  ip=IL[k].ip;jp=IL[k].jp;E0=IL[k].energy;
  part=0.;
  //3 possibilities, hairpin, ML, or another IL.
  if (HP_E(ip,jp)+E0<=HP_E(i,j)+EPSILON){
    if (HP_is_min(ip, jp)){
      if (debug_flag){
	printf("IL, hairpin, i, j, ip, jp are %d, %d %d %d.\n", i,j,ip, jp);
      }
      part+=exp(-(HP_E(ip,jp)+E0)/kT);
      if (ran<part){
	energy+=HP_E(ip,jp)+E0;
	return 0;
      }
    }
  }
  for (l=IL_start(ip,jp);l<=IL_end(ip,jp);l++)
    if (E0+IL[l].energy<=IL_Energy(i,j,IL[l].ip,IL[l].jp, S0)+EPSILON){
      part+=IL[l].part*exp(-E0/kT);
      if (ran<part){
	energy+=E0;
	ipp=IL[l].ip;jpp=IL[l].jp;
	if (debug_flag){
	  printf("before IL, ipp, jpp, %d, %d %d %d.\n", ipp, jpp);
	  printf("X<%s>X\n", structure);
	}
	structure[ipp]='(';structure[jpp]=')';
	if (debug_flag){
	  printf("IL, ipp, jpp, %d, %d.\n", ipp, jpp);
	  int jdeb;
	  printf("X<%s>X\n", structure);
	  int lb=0, rb=0;
	  for (jdeb=0;jdeb<=n;jdeb++){
	    if (structure[jdeb]=='(') lb++;
	    if (structure[jdeb]==')') rb++;
	  }
	  printf("%d, %d.\n", lb, rb);
	}
	backtrack_IL(ip, jp, l);
	return 0;
      }
    }
  if (delta_MLCmin(i,j,ip,jp)){
    part+=exp(-E0/kT)*Z_MLC[ip][jp];
    if (ran<part){
      energy+=E0;
      if (debug_flag){
	printf("IL, hairpin, i, j, ip, jp are %d, %d, %d, %d.\n", i,j,ip, jp);
      }
      backtrack_MLC(ip,jp);return 0;
    }
  }
  if (debug_flag){
    printf("got to end of IL.\n");
  }
  return 1;
}
