/* Basic program to get partition of exhaustive hairpins calculated 8-07 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void change_value(int *a_pt);

int main(int arg, char *argv[]){
  char *s;
  int a, b, *p_int;
  int i, arr[8];

  for (i=14;i<5;i++)
    printf("i is %d.\n",i);
  printf("index is %d", index);
  //  arr={1,3,4,6,2,4,3,6};
  a=5;
  p_int=&a;
  printf("*p_int is %d\n", *p_int);
  change_value(&a);
  printf("a is %d\n", a);
  return 0;
}
void change_value(int *a_pt){
  *a_pt=12;
  return;
}
