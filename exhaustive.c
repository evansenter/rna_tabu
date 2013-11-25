/* Basic program to get partition of exhaustive hairpins calculated 8-07 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fold.h"
#include "energy_const.h"
#include "fold_vars.h"
#include "pair_mat.h"

#define PRIVATE static
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

PRIVATE int num_bp, **pair_list, *pairTable, n, count_NJ=0;
PRIVATE char sequence[200], structure[200], *pt_structure;
PRIVATE short *seq_num;
PRIVATE float kT;
PRIVATE double partition;
PRIVATE char dummy[200];
PRIVATE int debug=0;

double chk_if_min();
void get_next_structure();
int chk_add_bp(int i_in, int j_in);
int chk_next_same_num_bp();
PRIVATE short *encode_seq(const char *seq);

int main(int arg, char *argv[]){
  float energy, part;
  int i, j, count, count_all;

  printf("energy_set is %d.\n", energy_set);
  make_pair_matrix();//needed for pair matching

  //  temperature=150.;
  kT = (temperature+K0)*GASCONST/1000.0;
  printf("Enter nucleotide sequence\n");
  scanf("%s", sequence);
  n=strlen(sequence);
  for (i=0; i<n; i++)
    sequence[i] = toupper((int)sequence[i]);
  seq_num = encode_seq(sequence);
  for (i=0; i<=n; i++){
    printf("%d", seq_num[i]);
  }
  printf("length of sequence is %i.\n", strlen(sequence));
  for (i=0; i<n; i++){
    structure[i]='.';
  }
  structure[n]='\0';
  printf("structure is:\n %s\n", structure);
  printf("length of structure is %i.\n", strlen(structure));
  pt_structure=structure-1; //Index from 1.
  pair_list=(int **) malloc((n/2+1) * sizeof (int *));
  //  printf("got here");
  for (i=0; i<(n/2+1); i++){
    pair_list[i]=(int *) malloc(2*sizeof(int));
  }
  pairTable=(int *) calloc(n+1, sizeof(int));
  dangles=0;
  num_bp=0;
  partition=0.0;
  count=0;count_all=0;
  do {
    get_next_structure();count_all++;
    //    printf("next structure is \n%s.\n%s\n", structure, sequence);
    //    getchar();
    if (part=chk_if_min()) {
      count++;
      printf("%s - %f\n", structure, part);
    }
  } while (num_bp>0);
  printf("partition is %f.\n", partition);
  printf("number of minimal is %d.\n", count);
  printf("number of structures is %d.\n", count_all);
  printf("number of NJ minimal is %d.\n", count_NJ);
  return 0;
}

void get_next_structure(){
  if (chk_add_bp(0,0)){
    return;
  } else
    while (num_bp>0 && !chk_inc_bp())
      ;
  //the num_bp>0 causes the last structure enumerated to be the empty one.
  return;
}

double chk_if_min(){
  float E, E0;
  int i, j, k, NJ_flag;

  
  E0=energy_of_struct(sequence, structure);
  //Add arbitrary base pair.
  NJ_flag=1;
  for (i=1;i<=n-4;i++){
    if (pairTable[i]) continue;
    j=i+1;
    while (j<=n){
      if (pt_structure[j]=='('){
	j=pairTable[j]+1;	
      }  else if (pt_structure[j]==')'){
	break;
      } else if (j<i+4){
	j++;
      } else if (pair[seq_num[i]][seq_num[j]]){
	NJ_flag=0;
	pt_structure[i]='(';pt_structure[j]=')';
	E=energy_of_struct(sequence, structure);
	pt_structure[i]='.';pt_structure[j]='.';
	if (E<E0) {
	  return 0;
	}
	j++;	
      } else {
	j++;
      }
    }
  }

  if (NJ_flag){
    count_NJ++;
  }
  //Remove arbitrary base pair.
  for (k=1; k<=num_bp; k++){
    i=pair_list[k][0];j=pair_list[k][1];
    pt_structure[i]='.';pt_structure[j]='.';
    E=energy_of_struct(sequence, structure);
    pt_structure[i]='(';pt_structure[j]=')';
    if (E<E0) {
      return 0;
    }
  }

  //  printf("locally min structure is %s\n", structure);
  //  printf("adding %f to partition.\n ", exp((double)-E0/kT));
  partition+=exp( (double) -E0/kT);
  //  printf ("energy, partition are %f, %f\n",E0, (float) partition);
  return exp( (double) -E0/kT);
}

int chk_add_bp(int i_in, int j_in){
  int i0, j0, i, j;
  i=i_in;j=j_in;

  if (i==0){
    if (num_bp>0){
      i0=pair_list[num_bp][0];
      j0=pair_list[num_bp][1];
    } else {
      i0=0;
      j0=n+1;
    }
    i=i0+1;j=i0+2;
  }
  while (i<n){
    if (pairTable[j] || pairTable[i] || j>n){
      i++;j=i+1;
    } else if (j<i+4)
      j++;
    else if (pair[seq_num[i]][seq_num[j]]){
      num_bp++;
      pair_list[num_bp][0]=i;pair_list[num_bp][1]=j;
      pt_structure[i]='(';pt_structure[j]=')';
      pairTable[i]=j;pairTable[j]=i;
      return 1;
    } else j++;
  }
  return 0;
}

int chk_inc_bp(){
  int i0, j0,i,j, border, border2;
  
  if (num_bp==0){
    printf ("problem, looks like no structures.");
    return 1;
  }
  i0=pair_list[num_bp][0];
  j0=pair_list[num_bp][1];
  pt_structure[i0]='.';pt_structure[j0]='.';
  pairTable[i0]=0;pairTable[j0]=0;
  num_bp--;
  //remove base pair and then see if we can add another one..
  return (chk_add_bp(i0,j0+1));
}


//Below modified from fold.c
PRIVATE short *encode_seq(const char *seq) {
  
  unsigned int i,l;
  short *seq_num_out;
  
  l = strlen(seq);
  seq_num_out = (short *) space(sizeof(short)*(l+2));
  seq_num_out[0]=l;
  for (i=1; i<=l; i++) { /* make numerical encoding of seq */
    seq_num_out[i]= (short) encode_char(toupper(seq[i-1]));
  }
  return seq_num_out;
}
      
