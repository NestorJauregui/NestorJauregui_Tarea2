 56#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "shock.h"


int main(){

 double dx = L/n_dis;
 double dt = 1E-6;


 double *rho = malloc((n_dis)*sizeof(double));
 double *P = malloc((n_dis)*sizeof(double));
 double *u = malloc((n_dis)*sizeof(double));

 double *U = malloc((n_dis*3.0)*sizeof(double));
 double *UB = malloc((n_dis*3.0)*sizeof(double));
 double *UA = malloc((n_dis*3.0)*sizeof(double));

 double *F = malloc((n_dis*3.0)*sizeof(double));
 double *FB = malloc((n_dis*3.0)*sizeof(double));
 double *FA = malloc((n_dis*3.0)*sizeof(double));


 init_array(U,UB,UA);

 calc_F(U,F);





return 0;
}



void init_array( double *U, double *UB,double *UA){
  int i;
  int j;
  
  for(i=0; i<n_dis; i++){
    if(i<=(n_dis-1)/2){
      U[matrix(0,i)] = rho1;
      U[matrix(1,i)] = rho1 * U1;
      U[matrix(2,i)] = (P1/(GAMMA+1))+rho1*U1*U1/2;
      UA[matrix(0,i)] = rho1;
      UA[matrix(1,i)] = rho1 * U1;
      UA[matrix(2,i)] = (P1/(GAMMA+1))+rho1*U1*U1/2;
      UB[matrix(0,i)] = rho1;
      UB[matrix(1,i)] = rho1 * U1;
      UB[matrix(2,i)] = (P1/(GAMMA+1))+rho1*U1*U1/2;
     }
  }
}

void calc_F( double *U, double *F){

  for(i=0; i<n_dis ; i++){
    
    F[matrix(0,i)] = U[matrix(1,i)];
    F[matrix(2,i)] = pow(U[matrix(1,i)],2)/U[matrix(0,i)] + (GAMMA-1) * (U[matrix(2,i)] -0.5* pow(U[matrix(1,i)],2)/U[matrix(0,i)]);
    F[matrix(3,i)] = (U[matrix(1,i)] / U[matrix(0,i)]) *( U[matrix(2,i)] + (GAMMA-1) * (U[matrix(2,i)] - 0.5 * pow(U[matrix(1,i)],2)/U[matrix(0,i)]) ) ;
    
  }

}

void lax(double *U, double *F,double *UA, double *FA,double *UB, double *FB){


  for(i=0;i<n_dis-1;i++){
    for(j=0;j<3;j++){
      for(){

      }
      UB[matrix(j,i)] = U[matrix(j,i)]-(dt/dx)(F[matrix(j,i+1)]-F[matrix(j,i)]);
      
      UA[matrix(j,i)] = 0.5*(U[matrix(j,i)] + UB[matrix(j,i)] - (dt/dx)());

    }
  }
}


int matrix(int fila, int columna){
  return (n_dis*fila)+columna;
}
