struct eReturn{
  int start, end, max;
};
typedef struct eReturn eReturn;

/* function from fold.c */
extern float  fold(const char *sequence, char *structure);
/* calculate mfe-structure of sequence */
extern float  energy_of_struct(const char *string, const char *structure);
/* calculate energy of string on structure */
extern void   free_arrays(void);           /* free arrays for mfe folding */
extern void   initialize_fold(int length); /* allocate arrays for folding */
extern void   update_fold_params(void);    /* recalculate parameters */
extern char  *backtrack_fold_from_pair(char *sequence, int i, int j);
extern float energy_of_circ_struct(const char *string, const char *structure);
extern int stack_energy_ext(int i, const char *string);
extern short encode_seq_ext(const char *sequence);
extern int ML_Energy_fold(int i, int is_extloop);
extern eReturn getBarrierEnergy(char *sequence, char *start_str, char *end_str, short route[1000][500], int *rl, int init_w, int th);
extern int local_energy(int i, const char *string);
extern int get_best_step(short* start, short* add, short* subtract, 
			 char *sequence,
			 int* return_index, int* addflag, int* emin, int **lock, int counter, float w, int nc, int bc, int emax, int e2);
int can_add(short* start, short* add, const int i);
int find_energy_add(short* start, short* add, const int i, char *sequence);
int find_energy_add2(short* start, const int i, const int j, char *sequence);
int find_energy_sub(short* start, const int i, const int j, char *sequence);
int disp_ptable(short*, char* name);
int add_ptable_bp(short*, const int, const int);
int remove_ptable_bp(short*, const int, const int);
int get_energy(short*, char*);


