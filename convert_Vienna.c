/* Program to convert from my standard notation for energy of IL's, HP's and
ML's from existing Vienna package subroutines. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fold.h"
#include "energy_const.h"
#include "fold_vars.h"
#include "pair_mat.h"
#include "convert_Vienna.h"
#include "params.h"

#define PRIVATE static

void Initialize_Params(){
  P = scale_parameters();//from params.c, gets our parameters, for a given temp
}

float HP_Energy(int i, int j, short *S0, char* sequence){
  int type, energy_int;

  type=pair[S0[i]][S0[j]];
  energy_int=HairpinE(j-i-1, type, S0[i+1], S0[j-1], sequence+i-1);
  //  printf("i,j,HP_energy are %d,%d,%d.\n", i, j, energy_int);
  return ((float) energy_int)/100.;
}
float IL_Energy(int i, int j, int ip, int jp, short *S0){
  int type, type2, energy_int, k, l;

  type=pair[S0[i]][S0[j]];
  type2=pair[S0[jp]][S0[ip]];
  energy_int=LoopEnergy(ip-i-1, j-jp-1, type, type2,
			S0[i+1], S0[j-1], S0[ip-1], S0[jp+1]);
  //  printf ("energy_int is now %d\n", energy_int);
  //getchar();
  return ((float) energy_int)/100.;
}  
float EL_Energy(int i, int j, short *S0){
  int type, energy_int;

  type=pair[S0[i]][S0[j]];
  energy_int=P->MLintern[type]-P->MLintern[1]; /* 0 or AU penalty */
  return (float) energy_int/100.;
  //return 3.5;
}
float ML_Energy(int i, int j, short *S0){
  int type, energy_int;
  int k;
  
  type=pair[S0[i]][S0[j]];
  energy_int=P->MLintern[type]; 
  //  printf("i,j,EL_energy are %d,%d,%d.\n",i,j,energy_int);
  return (float) energy_int/100.;
}
			
  
