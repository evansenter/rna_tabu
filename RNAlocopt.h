#define EPSILON 0.00001
#define maxSeqLen 500
#define tailses 3
#define PRIVATE static
#define PUBLIC
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

PUBLIC int n, Tails;
PUBLIC double kT;
PUBLIC char sequence[1000];
PUBLIC short *S0;

PUBLIC double ML_base, ML_close;
PUBLIC double Z_ML[maxSeqLen+1][maxSeqLen+1][tailses+2][tailses+2], Z_MLC[maxSeqLen+1][maxSeqLen+1];
PUBLIC double Z_M1p[maxSeqLen+1][maxSeqLen+1][tailses+2], Z_M1q[maxSeqLen+1][maxSeqLen+1][tailses+2], Z_M1long[maxSeqLen+1][maxSeqLen+1];
PUBLIC double Z_EL_1[maxSeqLen+1][tailses+2], count_EL_1[maxSeqLen+1][tailses+2];
PUBLIC double Z_M1p_E[maxSeqLen+1][maxSeqLen+1][tailses+2], Z_M1q_E[maxSeqLen+1][maxSeqLen+1][tailses+2], \
  Z_M1long_E[maxSeqLen+1][maxSeqLen+1];
PUBLIC double Zstar[maxSeqLen+1][maxSeqLen+1], ZEstar[maxSeqLen+1][maxSeqLen+1];

PRIVATE void Induct_Zstar(int i,int j);
PRIVATE void Induct_M1(int i,int j);
PRIVATE void Induct_ML_EL(int i,int j);
PUBLIC double Z_M1_val(int i, int j, int p, int q);
PRIVATE int min_ML(int i,int j,int p,int q);
PUBLIC double Z_M1_E_val(int i,int j,int p,int q);
PRIVATE int min_EL(int i,int j,int p,int q);
PRIVATE void Induct_MLC(int i, int j);
PUBLIC int delta_MLCmin(int i, int j, int ip, int jp);
PUBLIC int delta_min(int i, int j, int p, int q);
PRIVATE void Induct_IL(int i,int j);
PRIVATE double Find_part_IL(int i, int j, int ip, int jp,double *count_pt,double E0);
PRIVATE int Not_minimal_IL(int i, int j, int ip, int jp, double E0);


struct bpStruct {
  int left[maxSeqLen/2], right[maxSeqLen/2];
  int size;
};

PRIVATE int Get_Edit_Distance_from_faa(char* faaStr1, char* faaStr2);
PRIVATE void basepairs_from_faa(struct bpStruct *bpStr, char* faaStr);
PRIVATE int Get_Edit_Distance_from_bp(struct bpStruct *bpStr1, struct bpStruct *bpStr2);

