/* Program to use recursion to get partition 8-07 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "fold.h"
#include "energy_const.h"
#include "fold_vars.h"
#include "pair_mat.h"
#include "convert_Vienna.h"
#include "IL.h"
#include "params.h"
#include "HP.h"
#include "RNAlocopt.h"
#include "backtrack.h"

PRIVATE void initialize(int argc, char *argv[]);
PRIVATE short *encode_seq(const char *seq);

int main(int argc, char *argv[]){
  int j_minus_i, i, j, q;
  double total_part, total_count;
  int fin_index;

  char struc[100][500], struc1[maxSeqLen], struc2[maxSeqLen];
  char smpStruc[1000];
  int sample_size;
  
  int dstart, dend, di, lb, rb;
  double dtest, dtest2;

  initialize(argc, argv);
  /* This simply moves most initialization away from the front
     of the code, so that the main subroutine is more compact
     This initializes several globals (that I should be able to 
     switch to variables just visible by module RNAlocopt.c),
     n, the sequence length
     kT,
     sequence, from argv[1]
     S0, the sequence in a format that Vienna package, more specifically
         fold.c, likes.
     Tails, we just copy this from tailses, defined at the top of RNAlocopt.h,
         (This is done so we can assign arrays without malloc, ironically
	 because this worked a lot better with the debugger DDD.)  Tails
	 gives the length of tails extending into a multiloop to store. 
	 To get exactly the right answer, use 10 (in definition of tailses
	 of RNAlocopt.h) but 3 gives very close to the right answer and is much
	 faster    
  */

  if (argc>2){
    printf("argv[2] is %s\n", argv[2]);
    sample_size=atoi(argv[2]);
    printf("sample size is %d.\n", sample_size);
  } else
    sample_size=10;
  /* This is placed here, instead of in initialize, because sample_size 
     is defined local to this subroutine.
   */
  dangles=0;
  /* Declared in fold_vars.c.  
     As of 2-03-09, I'm almost sure dangles arent' implemented, this flag sets
     dangles=0 in fold.c, for consistency.  Set dangles=1 to use dangling
     energies, say when getting a return value from energy_of_struct().
  */
  for (i=0;i<=0;i++){
    printf("i is %d.\n", i);
  }
  printf("dangles is %d.\n", dangles);
  for (j_minus_i=4; j_minus_i<n;j_minus_i++)
    for (i=1; i<=(n-j_minus_i); i++){
      j=j_minus_i+i;
      //      printf("in main loop, i, j are %d, %d.\n", i, j);
      if (pair[S0[i]][S0[j]]){
	Induct_IL(i,j);
	Induct_Zstar(i,j);
	Induct_MLC(i,j);
      }
      Induct_M1(i,j);
      Induct_ML_EL(i,j);
      //getchar();
    }
  total_part=0.;
  total_count=0.;
  for (q=0;q<=Tails+1;q++){
    total_part+=Z_EL_1[n][q];
    total_count+=count_EL_1[n][q];
  }
  printf("total_part is %e.\n", total_part);
  printf("total_count is %e.\n", total_count);
  printf("length of sequence is %i.\n", strlen(sequence));
  fin_index=IL_index();
  printf("index maxSeqLen is %d.\n", fin_index);
  printf("IL portion of mem usage.....index/n^2 is %f.\n", 
	 ((double) fin_index)/ ((double) (n*n)) );
  for (i=0;i<sample_size;i++){
    float etemp;
    //    if (i==28)
    //      strcpy(smpStruc,perform_backtrack(1));
    //    else
      strcpy(smpStruc,perform_backtrack(0));
    //    if (i==28){
    int jdeb;
    //printf("<<%s\n>>", smpStruc);
    int lb=0, rb=0;
    for (jdeb=0;jdeb<=n;jdeb++){
      if (smpStruc[jdeb]=='(') lb++;
      if (smpStruc[jdeb]==')') rb++;
    }
    //printf("%d, %d.\n", lb, rb);
    //exit(0);
    //}
    etemp=energy_of_struct(sequence,smpStruc+1);
    printf("%s %8.3f\n", smpStruc+1, etemp);
  }

  return 0;
}     

int Get_Edit_Distance_from_faa(char* faaStr1, char* faaStr2){
  struct bpStruct bpStr1, bpStr2;
  
  basepairs_from_faa(&bpStr1, faaStr1);
  basepairs_from_faa(&bpStr2, faaStr2);
  return Get_Edit_Distance_from_bp(&bpStr1, &bpStr2);
} 

void basepairs_from_faa(struct bpStruct *bpStr, char* faaStr){
  int i, stackL[maxSeqLen/2], stackPos;

  stackPos=0;
  bpStr->size=0;
  for (i=1;i<strlen(faaStr);i++){
    if (faaStr[i]=='('){
      stackPos++;
      stackL[stackPos]=i;
    } else if (faaStr[i]==')') {
      bpStr->left[bpStr->size]=stackL[stackPos];
      bpStr->right[bpStr->size]=i;
      bpStr->size++;
      stackPos--;
    }
  }
  return;
}
int Get_Edit_Distance_from_bp(struct bpStruct *bpStr1, struct bpStruct *bpStr2){
  int i, j, match, mismatch, bp1R, bp2R;
  
  /*  printf("bp1 list is\n");
  for (i=0;i<bpStr1->size;i++){
    printf("(%d,%d)", bpStr1->left[i], bpStr1->right[i]);
  }
  printf("\n");
  printf("bp2 list is\n");
  for (i=0;i<bpStr2->size;i++){
    printf("(%d,%d)", bpStr2->left[i], bpStr2->right[i]);
  }
  printf("\n");
  */
  i=j=0;
  match=0;
  mismatch=0;
  while(i<bpStr1->size || j<bpStr2->size){
    bp1R=bpStr1->right[i];
    bp2R=bpStr2->right[j];
    if (i==bpStr1->size){
      j++; mismatch++;
    } else if (j==bpStr2->size){
      i++; mismatch++;
    } else if (bp1R==bp2R){
      if (bpStr1->left[i]==bpStr2->left[j])
	match++;
      else
	mismatch=mismatch+2;
      i++;j++;
    } else if (bp1R<bp2R){
      mismatch++; i++;
    } else {
      mismatch++; j++;
    }
    //printf("i, j, match, mismatch, %d, %d, %d, %d.\n", i, j, match, mismatch);
  }

  //printf("mismatch is %d.\n", mismatch);

  return mismatch;
}


void Induct_Zstar(int i, int j){
  int ip, jp, k, l, start, end, num_open_bp;
  int k_last, l_last, l_max;
  double EL_ij, ML_ij, E0, E_diff;
  
  int debug;
  
  start=IL_start(i,j);
  if (start==0) return;
  end=IL_end(i,j);
  ML_ij=ML_Energy(i,j,S0);
  EL_ij=EL_Energy(i,j,S0);
  //Get Zstar, ZEstar.
  for (k=start;k<=end;k++){
    E0=IL[k].energy;ip=IL[k].ip;jp=IL[k].jp;
    if (ML_ij+E0<ML_Energy(ip,jp,S0)+EPSILON){
      Zstar[i][j]+=IL[k].part*exp(-ML_ij/kT);
      Zstar[j][i]+=IL[k].count;
    }
    if (EL_ij+E0<EL_Energy(ip,jp,S0)+EPSILON){
      ZEstar[i][j]+=IL[k].part*exp(-EL_ij/kT);
      ZEstar[j][i]+=IL[k].count;
    }
  }
  return;
}

void Induct_M1(int i, int j){
  int p, q;

  if (j-i<1){
    printf("j<=i in Induct_M1.\n");
    exit(EXIT_FAILURE);
  }
  //Z_M1p has short tail on left, Z_M1q has short tail on right, Z_M1long has long tails
  for (p=0;p<=Tails;p++){
    Z_M1p[i][j][p]=exp(-ML_base/kT)*(Z_M1p[i][j-1][p]+Z_M1_val(i,j-1,p,Tails));
    Z_M1p[j][i][p]=Z_M1p[j-1][i][p]+Z_M1_val(j-1,i,p,Tails);
    Z_M1p_E[i][j][p]=Z_M1p_E[i][j-1][p]+Z_M1_E_val(i,j-1,p,Tails);
    Z_M1p_E[j][i][p]=Z_M1p_E[j-1][i][p]+Z_M1_E_val(j-1,i,p,Tails);
  }
  for (q=0;q<=Tails;q++){
    Z_M1q[i][j][q]=exp(-ML_base/kT)*(Z_M1q[i+1][j][q]+Z_M1_val(i+1,j,Tails,q));
    Z_M1q[j][i][q]=Z_M1q[j][i+1][q]+Z_M1_val(j,i+1,Tails,q);
    Z_M1q_E[i][j][q]=Z_M1q_E[i+1][j][q]+Z_M1_E_val(i+1,j,Tails,q);
    Z_M1q_E[j][i][q]=Z_M1q_E[j][i+1][q]+Z_M1_E_val(j,i+1,Tails,q);
  }
  Z_M1long[i][j]=exp(-ML_base/kT)*(Z_M1long[i][j-1]+Z_M1q[i][j-1][Tails]);
  Z_M1long[j][i]=Z_M1long[j-1][i]+Z_M1q[j-1][i][Tails];
  Z_M1long_E[i][j]=Z_M1long_E[i][j-1]+Z_M1q_E[i][j-1][Tails];
  Z_M1long_E[j][i]=Z_M1long_E[j-1][i]+Z_M1q_E[j-1][i][Tails];
  return;
}


void Induct_ML_EL(int i, int j){
  int p, q, k, a, minvar;
  double temp;
  
  for (p=0;p<=Tails+1;p++)
    for (q=0;q<=Tails+1;q++){
      for (k=i+4; k<=j-4;k++){
	minvar=MIN(Tails, k-i+1);
	for (a=0; a<=minvar; a++){
	  Z_ML[i][j][p][q]+=(Z_ML[i][k][p][a]+Z_M1_val(i,k,p,a))*Z_M1_val(k+1-a,j,a,q);
	  Z_ML[j][i][p][q]+=(Z_ML[k][i][p][a]+Z_M1_val(k,i,p,a))*Z_M1_val(j,k+1-a,a,q);
	}
	if (k+1-Tails>i){
	  Z_ML[i][j][p][q]+=(Z_ML[i][k][p][Tails+1]+Z_M1_val(i,k,p,Tails+1))*
	    Z_M1_val(k+1-Tails,j,Tails,q);
	  Z_ML[j][i][p][q]+=(Z_ML[k][i][p][Tails+1]+Z_M1_val(k,i,p,Tails+1))*
	    Z_M1_val(j,k+1-Tails,Tails,q);
	}
      }
    }
  if (i==1){
    for (q=0; q<=Tails+1; q++){
      for (k=0;k<=j-4;k++){
	minvar=MIN(Tails, k-i+1);
	for (a=0;a<=minvar;a++){
	  Z_EL_1[j][q]+=Z_EL_1[k][a]*Z_M1_E_val(k+1-a,j,a,q);
	  count_EL_1[j][q]+=count_EL_1[k][a]*Z_M1_E_val(j,k+1-a,a,q);
	}
	if (k+1-Tails>=1){
	  Z_EL_1[j][q]+=Z_EL_1[k][Tails+1]*Z_M1_E_val(k+1-Tails,j,Tails,q);
	  count_EL_1[j][q]+=count_EL_1[k][Tails+1]*Z_M1_E_val(j,k+1-Tails,Tails,q);
	}
      }
    }
  }
  return;
}

PUBLIC double Z_M1_val(i,j,p,q){
  int ip, jp;
  if (i<j){
    ip=i;jp=j;
  } else {
    ip=j;jp=i;
  }
  if ( (jp-q)-(ip+p)<4) return 0.;
  if (p==Tails+1){
    if (q==Tails+1)
      return Z_M1long[i][j];
    else
      return Z_M1q[i][j][q];
  } else if (q==Tails+1)
    return Z_M1p[i][j][p];
  // If not above, must check if minimal.
  if (!min_ML(ip,jp,p,q)) return 0.0;
  if (i<j)
    return exp(-ML_base*(p+q)/kT)*Zstar[i+p][j-q];
  else
    return Zstar[i-q][j+p];
}
 
PRIVATE int min_ML(i,j,p,q){
  int k, l, ip, jp;
  double ML_ipjp, diff;

  ip=i+p;jp=j-q;
  if (jp<=ip) return 0;//confusing check
  ML_ipjp=ML_Energy(ip,jp,S0);
  
  for (k=ip-1;k>=i;k--)
    for (l=jp+1;l<=j;l++){
      if (pair[S0[k]][S0[l]]){
	diff=IL_Energy(k,l,ip,jp,S0)+ML_Energy(k,l,S0)-ML_ipjp;
	if (diff<0.0-EPSILON) return 0;
      }
    }
  return 1;
}
PUBLIC double Z_M1_E_val(i,j,p,q){
  int ip, jp;
  if (i<j){
    ip=i;jp=j;
  } else {
    ip=j;jp=i;
  }
  if ( (jp-q)-(ip+p)<4) return 0.;
  if (p==Tails+1){
    if (q==Tails+1)
      return Z_M1long_E[i][j];
    else
      return Z_M1q_E[i][j][q];
  } else if (q==Tails+1)
    return Z_M1p_E[i][j][p];
  // If not above, must check if minimal.
  if (!min_EL(ip,jp,p,q)) return 0.0;
  if (i<j)
    return ZEstar[i+p][j-q];
  else
    return ZEstar[i-q][j+p];
}

PRIVATE int min_EL(i,j,p,q){
  int ip, jp, k, l;
  double EL_ipjp, diff;
  
  ip=i+p;jp=j-q;

  EL_ipjp=EL_Energy(ip,jp,S0);
  
  for (k=ip-1;k>=i;k--)
    for (l=jp+1;l<=j;l++){
      if (pair[S0[k]][S0[l]]){
	diff=IL_Energy(k,l,ip,jp,S0)+EL_Energy(k,l,S0)-EL_ipjp;
	if (diff<0.0-EPSILON) return 0;
      }
    }
  return 1;
}
 
PRIVATE void Induct_MLC(i,j){
  int p, q, prange, qrange;
  
  prange=MIN(Tails+1, (j-1)-(i+1)-3);
  for (p=0;p<=prange;p++){
    qrange=MIN(Tails+1, (j-1)-(i+1)-3-p);
    for (q=0;q<=qrange;q++){
      Z_MLC[i][j]+=delta_min(i,j,p,q)*Z_ML[i+1][j-1][p][q];
      Z_MLC[j][i]+=delta_min(i,j,p,q)*Z_ML[j-1][i+1][p][q];
    }
  }
  Z_MLC[i][j]*=exp(-(ML_close+ML_Energy(j,i,S0))/kT);
}

PUBLIC int delta_min(i,j,p,q){
  int k, l;
  double EMLC_ij;
  
  EMLC_ij=ML_Energy(j,i,S0);
  for (k=i+1;k<=i+p;k++)
    for (l=j-1;l>=j-q;l--)
      if (pair[S0[k]][S0[l]])
	if (ML_Energy(l,k,S0)+IL_Energy(i,j,k,l,S0)<EMLC_ij-EPSILON)
	  return 0;
  return 1;
}

PRIVATE void Induct_IL(int i, int j){
  double partition, exp_ext_energy, E0, count;
  int ip, jp;

  for (ip=i+1;ip<=MIN(i+30, j);ip++)
    for (jp=j-1; jp>=MAX( ip+4, j-(30-(ip-i)) ); jp-- ){
      if (!pair[S0[ip]][S0[jp]]) 
	continue;
      E0=IL_Energy(i,j,ip,jp, S0);
      if (Not_minimal_IL(i,j,ip,jp,E0)) 
	continue;
      partition=Find_part_IL(i, j, ip, jp, &count, E0);//count modified too.
      if (partition!=0.)
	IL_add(partition, E0, i,j,ip,jp,count);
    }
  return;
}

PRIVATE double Find_part_IL(int i, int j, int ip, int jp,
			   double *count_pt, double E0){
  double part, sum, count;
  int start, end, k, l;
  int k_last, l_last, l_max, num_open_bp;
  double E_diff, ML_jpip;
  part=0.;
  count=0;
  if (HP_E(ip,jp)+E0<=HP_E(i,j)+EPSILON){
    if (HP_is_min(ip, jp)){
      part+=exp(-(HP_E(ip,jp)+E0)/kT);
      count=count+1.;
    }
  } 
  sum=0.;
  start=IL_start(ip, jp);
  if (start){
    end=IL_end(ip,jp);
    for (k=IL_start(ip,jp);k<=end;k++)
       if (E0+IL[k].energy<=IL_Energy(i,j,IL[k].ip,IL[k].jp, S0)+EPSILON){
	sum+=IL[k].part;
	count+=IL[k].count;
      }
  }
  if (sum)
    part+=sum*exp(-E0/kT);

  //Now for ML part.

  if (delta_MLCmin(i,j,ip,jp)){
    part+=exp(-E0/kT)*Z_MLC[ip][jp];
    count+=Z_MLC[jp][ip];
  }
  *count_pt=count;
  return part;
}

PUBLIC int delta_MLCmin(i,j,ip,jp){

  if (ML_Energy(jp,ip,S0)+IL_Energy(i,j,ip,jp,S0)>ML_Energy(j,i,S0)+EPSILON)
	return 0;
  return 1;
}
  

PRIVATE int Not_minimal_IL(int i, int j, int ip, int jp, double E0){
  int k, l;

  for (k=i+1;k<ip;k++)
    for (l=j-1;l>jp;l--)
      if (pair[S0[k]][S0[l]])
	if (  (IL_Energy(i,j,k,l,S0)+IL_Energy(k,l,ip,jp,S0) )+EPSILON<E0)
	  return 1;
  return 0;
}

PRIVATE void initialize(int argc, char *argv[]) {
  int i, k, l, m;

  Tails=tailses;
  kT = (temperature+K0)*GASCONST/1000.0;
  /* Removing interactive option as of 2-03-09
     if (argc<2){
     printf("Enter nucleotide sequence\n");
     scanf("%s", sequence);
     } else */
  if (argc==2 || argc==3) {
    strcpy(sequence, argv[1]);
  } else {
    printf("Usage:\n" );
    printf(" %s [RNA_sequence] [sample_size]\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  n = (int) strlen(sequence);
  if (n>maxSeqLen){
    printf("Limited to sequences of length %d.  Sorry too long.", maxSeqLen);
    exit(EXIT_FAILURE);
  }
  for (i=0; i<n; i++)
    sequence[i] = toupper((int)sequence[i]);
  printf("sequence is %s.\n", sequence);
  printf("maxSeqLen, n, %d, %d.\n", maxSeqLen, n);
  //strcpy(sequence, "ACGUACGUACGU");
  //strcpy(sequence, "CCCAAAAAAAAAAAGGG");
  if (n!=maxSeqLen) printf("problem\n");

  S0 = encode_seq(sequence);
  /* S0 is a variable in vienna package format.  It's a global variable 
     declared in RNAlocopt.h, but it's in the format needed for several routines
     in fold.c.  encode_seq is actually a routine that is
     copied from fold.c (it's private to fold.c, so a copy was made here, so
     that no changes need be made to fold.c).  
  */  
  make_pair_matrix();//needed for pair matching
  /* That is, this is used in fold.c, and here, to initialize some matrices
     that are used to evaluate pair matching for members of the arraay S0 
  */
  for (l=0; l<=Tails; l++){
    Z_EL_1[l][l]=1.;
    count_EL_1[l][l]=1.;
  }
  for (l=Tails+1; l<=n; l++){
    Z_EL_1[l][Tails+1]=1.;
    count_EL_1[l][Tails+1]=1.;
  }
  /*  Above are initial conditions on the partition function.
   */
  for (i=0; i<=n; i++){
    printf("%d", S0[i]);
  }
  printf("\n");
  printf("initializing params...\n");
  Initialize_Params();
  /* This subroutine is just a 1-liner in convert_Vienna.c.  convert_Vienna.c
     is code written by me, Andy Lorenz, to change the format of Vienna into
     something I liked.  Looking at this code, it becomes clear how to use
     the original Vienna subroutines as well.  Initialize_Params imports some 
     parameters from params.c to be used in the module convert_Vienna.c.
  */
  ML_base=(double) P->MLbase/100.;
  ML_close=(double) P->MLclosing/100.;//only locally used params.
  IL_initialize(n);
  /* Internal loops are stored in a table, here we malloc the necessary space,
     this space grows dynamically as the program runs, see IL.c
  */
  find_irred_HP(S0,sequence,kT,n);
  /* This is another initial condition.  We determine exhaustively all irreducible
     hairpins.  See HP.c, there's very little code there.
  */
  return;
}

//Below modified from fold.c
PRIVATE short *encode_seq(const char *seq) {
  
  unsigned int k,l;
  short *S0_out;
  
  l = strlen(seq);
  S0_out = (short *) space(sizeof(short)*(l+2));
  S0_out[0]=l;
  for (k=1; k<=l; k++) { /* make numerical encoding of seq */
    S0_out[k]= (short) encode_char(toupper(seq[k-1]));
  }
  return S0_out;
}
