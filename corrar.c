#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void cambia(double *a, double *b, double *c);

int main(){
  
  double *d = malloc(2*sizeof(double));
  double *e = malloc(2*sizeof(double));
  double *f = malloc(2*sizeof(double));

  d[0] = 1.0;
  d[1] = 1.58;

  e[0] = 1.0;
  e[1] = 1.58;

  f[0] = 0.0;
  f[1] = 0.0;
  

  cambia(d,e,f);
  printf("%f\n" , f[0]);
  printf("%f\n" , f[1]);

  return 0;
}


void cambia(double *a, double *b, double *c){

  int i;
  for(i=0;i<2;i++){
    c[i] = a[i]+b[i];
  }

}
