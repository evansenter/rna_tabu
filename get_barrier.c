
#include "get_barrier.h"
#include "foldMod.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "utils.h"

#define MAXLINE 3000

int disp_ptable2(short intable[500]){
  
  int x;
  short local[500];
  
  memcpy(local,intable,500*sizeof(short));
  
//  for (x=0;x<=local[0];x++)
  //  printf(" %d ",local[x]);

  for (x=1;x<=local[0];x++){
    if (local[x]==0)
      printf(".");
    else if (local[x]>x){
      printf("(");
    //  intable[local[x]]=1;
    }else
      printf(")");
  }
  printf("\n");

}

int main(int argc, char *argv[]){
  char sequence[500], start_str[500], end_str[500];
  int energy, e1, e2, barrier, tmp;
  eReturn eret;
  //for solution storage
  short route[1000][500];
  short best_route[1000][500];
  int best_energy=1000000;
  int k,i,t,j,l,a;
  int n;
  int route_len;
  int best_route_len;
  int best_k;
  int reps,st;
  int k1,k2;
  int th;
  FILE *fp;
  char lines[4][MAXLINE], fname[100], comment[100];
  int fromFAA=0;

  if (argc==8){
    strcpy(sequence, argv[1]);
    strcpy(start_str, argv[2]);
    n = strlen(sequence);
    energy_of_struct(sequence, start_str);//May initialize some things.
    strcpy(end_str, argv[3]);
    reps = atoi(argv[4]);
    k1= atoi(argv[5]);
    k2 = atoi(argv[6]);
    th = (int)(atof(argv[7])*100);
  } else if (argc==7 && strcmp(argv[1], "-f")==0){
    fromFAA=1;
    strcpy(fname, argv[2]);
    fp = fopen(fname, "r");
    for (i=0;i<4;i++){
     if (fgets(lines[i],MAXLINE,fp)==NULL){
	printf("Cannot open fasta file.\n  Check format, need start and stop structure after sequence.\n", fname);
	exit(1);
      }
    }
    strcpy(comment, lines[0]);
    sscanf(lines[1], "%s", sequence);
    sscanf(lines[2], "%s", start_str);
    sscanf(lines[3], "%s", end_str);
    n = strlen(sequence);
    energy_of_struct(sequence, start_str);//May initialize some things.
    reps = atoi(argv[3]);
    k1= atoi(argv[4]);
    k2 = atoi(argv[5]);
    th = (int)(atof(argv[6])*100);
  } else {
    printf("Usage:\n");
    printf(" %s sequence start_structure end_structure iterations lb_init_weight ub_init_weight energy_bound\n", argv[0]);
    printf("    (use parenthesis around structures)\n");
    printf(" or \n");
    printf(" %s -f fasta_file_with_start_end_structure iterations lb_init_weight ub_init_weight energy_bound\n", argv[0]);
    exit(1);
  }
  for(t=0;t<reps;t++){
     best_energy=1000000; 
     for(k=k1;k<k2;k+=2){  
      srand( (unsigned int) time( NULL )+t+k);
      //srand( (unsigned int) reps);
      //printf("Try for k = %d\n",k);  
      
      eret=getBarrierEnergy(sequence, start_str, end_str, route, &route_len,k,th);
      barrier=eret.max;
      st= eret.start;
      //printf(" Barrier1 = %d (max = %d start = %d\n", barrier,eret.max,eret.start);
      //for(i=0;i<route_len;i++)
      // disp_ptable2(route[i]);
      //getchar();
      if(barrier<best_energy){
         best_energy=barrier;
         best_route_len = route_len;
         best_k=k;
         for(i=0;i<route_len;i++)
           memcpy(best_route[i],route[i],(n+1)*sizeof(short));
      }
      eret=getBarrierEnergy(sequence, end_str, start_str, route, &route_len,k,th);
      barrier=eret.max;
      //printf(" Barrier2 = %d (max = %d start = %d\n", barrier,eret.max,st);
      //printf(" Barrier2 = %d \n", barrier);
      //for(i=0;i<route_len;i++)
      // disp_ptable2(route[i]);
      //getchar();
      if(barrier<best_energy){
         best_energy=barrier;
         best_route_len = route_len;
         best_k = k;
         for(i=0;i<route_len;i++)
           memcpy(best_route[i],route[route_len-i-1],(n+1)*sizeof(short));
      }
    }
   
    if(best_energy>50000)
      printf("No path found\n");
    else{
      if (fromFAA){
	printf("%s", comment);
      }
      printf("%s\n",sequence);
      for(i=0;i<best_route_len;i++)
       disp_ptable2(best_route[i]);
      printf("barrier is %5.2f -- route length = %d -- best k = %d th = %d\n", (float) best_energy/100.,best_route_len,best_k,th);
    //getchar();
    }
  }
  return 0;
}

