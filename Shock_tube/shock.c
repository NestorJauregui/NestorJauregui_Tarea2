#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "shock.h"


int main(){

  double dx = L/((float)(n_dis)-1);
  double dt = dx/3.0;
  int i;
  double rhoff;
  double Pff;
  double uff;
  int iter;

  double *rho = malloc(n_dis*sizeof(double));
  double *P = malloc(n_dis*sizeof(double));
  double *u = malloc(n_dis*sizeof(double));
  double *c = malloc(n_dis*sizeof(double));


  double *U = malloc((n_dis*3.0)*sizeof(double));
  double *UB = malloc((n_dis*3.0)*sizeof(double));
  double *UA = malloc((n_dis*3.0)*sizeof(double));

  double *F = malloc((n_dis*3.0)*sizeof(double));
  double *FB = malloc((n_dis*3.0)*sizeof(double));
  double *FA = malloc((n_dis*3.0)*sizeof(double));


 init_array(U,UB,UA);

 lax(U,F,UA,UB,FB,c,dt,dx);

 calc_var(U,rho,u,P);

 
 FILE *salida = fopen("datos", "w+");

 for(i=0; i<n_dis; i++)
       {
      rhoff = rho[i];
      uff = u[i];
      Pff = P[i];
      fprintf(salida, "%f %f %f \n", rhoff, uff, Pff);
    }
 fclose(salida);


   
 

return 0;

}

/*inicializo la matriz U y sus asociadas*/
void init_array( double *U, double *UB,double *UA){
  int i;
  int j;
  
  for(i=0; i<n_dis; i++){
    if(i<=(n_dis-1)/2){
      U[matrix(0,i)] = rho1;
      U[matrix(1,i)] = rho1 * U1;
      U[matrix(2,i)] = (P1/(GAMMA-1))+rho1*U1*U1/2;
      UA[matrix(0,i)] = rho1;
      UA[matrix(1,i)] = rho1 * U1;
      UA[matrix(2,i)] = (P1/(GAMMA-1))+rho1*U1*U1/2;
      UB[matrix(0,i)] = rho1;
      UB[matrix(1,i)] = rho1 * U1;
      UB[matrix(2,i)] = (P1/(GAMMA-1))+rho1*U1*U1/2;
     }
    if(i>(n_dis-1)/2){

      U[matrix(0,i)] = rho5;
      U[matrix(1,i)] = rho5 * U5;
      U[matrix(2,i)] = (P5/(GAMMA-1))+rho5*U5*U5/2;
      UA[matrix(0,i)] = rho5;
      UA[matrix(1,i)] = rho5 * U5;
      UA[matrix(2,i)] = (P5/(GAMMA-1))+rho5*U5*U5/2;
      UB[matrix(0,i)] = rho5;
      UB[matrix(1,i)] = rho5 * U5;
      UB[matrix(2,i)] = (P5/(GAMMA-1))+rho1*U5*U5/2;

    }


  }
}

/*Calcula la matriz F a partir de U unicamente*/
void calc_F( double *U, double *F){
  int i;
  for(i=0; i<n_dis ; i++){
    
    F[matrix(0,i)] = U[matrix(1,i)];
    F[matrix(1,i)] =( pow(U[matrix(1,i)],2.0)/U[matrix(0,i)] ) + (GAMMA-1) * (U[matrix(2,i)] - 0.5* pow(U[matrix(1,i)],2.0)/U[matrix(0,i)]);
    F[matrix(2,i)] = (U[matrix(1,i)] / U[matrix(0,i)]) *( U[matrix(2,i)] + (GAMMA-1) * (U[matrix(2,i)] - 0.5 * pow(U[matrix(1,i)],2)/U[matrix(0,i)]) ) ;
    
  }

}

/*Calcula u estrella en el algoritmo de Lax W. y me genera f estrella tambien*/
void calc_star(double *U,double *UB,double *F,double *FB, double dt, double dx)
{
  int i;
  int j;
  for(i=0;i<n_dis-1;i++){
    for(j=0;j<3;j++){
      UB[matrix(j,i)] = (U[matrix(j,i)]+U[matrix(j,i+1)] - (dt/dx) * ( F[matrix(j,i+1)]-F[matrix(j,i)]))/2.0;
    }
  }
  calc_F(UB,FB);  
}

/*calcula el avance temporal de U y actualiza U* a U*/

void calc_UA(double *U, double *UA, double *UB, double *FB, double dt, double dx )
{
    int i;
    int j;

    for(i=1;i<n_dis;i++)
    {
      for(j=0;j<3;j++)
	{
	  UA[matrix(j,i)] = U[matrix(j,i)] -(dt/dx)*(FB[matrix(j,i)]-FB[matrix(j ,i-1)]);
	}
    }

    for (i = 0; i < n_dis; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  U[matrix(j,i)] = UA[matrix(j,i)];
	}
    }
}

/*calcula la velocidad maxima en un tiempo fijo*/
double calc_umax(double *U,double *c)
{
  int i;
  double temp=0.0;
  double temp2=0.0;
  double umax=0.0;
  for ( i = 0; i < n_dis; i++)
    {
      temp = U[matrix(1,i)]/U[matrix(0,i)];
      temp2 = c[i];
      if (temp+temp2 > umax)
	{
	  umax = temp+temp2;
	}
    }
  return umax;
}

/*calculala velocidad del sonido en cada punto*/

void *s_speed(double *U, double *c){
  
  double rh;
  double pres;
  int i;
  for(i=0; i<n_dis; i++ ){
    
    rh = U[matrix(0,i)];
    pres = (GAMMA-1) * (U[matrix(2,i)] - 0.5* pow(U[matrix(1,i)],2.0)/U[matrix(0,i)]) ;

    c[i] = pow(GAMMA * pres/rh ,0.5);

  }

}


/*algoritmo Lax Wendroff*/

void lax(double *U, double *F,double *UA,double *UB, double *FB, double *c, double dt, double dx)
{

  double t = 0.0;
  while(t<time){ 

    calc_F(U,F);
    calc_star(U,UB,F,FB,dt,dx);
    calc_UA(U,UA,UB,FB,dt,dx);
    s_speed(U,c);
    if(calc_umax(U,c) > 0.0){
      dt = 0.5*dx/calc_umax(U,c);
    }
    printf(" %f  %f  %f   %f  %f  %f  %f  %f  %f  %f  %f \n",  F[matrix(0,245)], F[matrix(0,355)], F[matrix(0,465)], F[matrix(1,135)], F[matrix(1,245)], F[matrix(1,465)], F[matrix(2,135)], F[matrix(2,245)], F[matrix(2,355)], calc_umax(U,c), t);
    t = t + dt;
  }

}

/*calcula las variables finales*/

void calc_var(double *U, double *rho, double *u, double *P)
{
	int i;
	for (i = 0; i < n_dis; i++)
	{
	  rho[i] = U[matrix(0,i)];
	  u[i] = U[matrix(1,i)]/U[matrix(0,i)];
	  P[i] = (GAMMA-1) * (U[matrix(2,i)] -0.5* pow(U[matrix(1,i)],2)/U[matrix(0,i)]);

	}

}

/*para usar notacion de matriz en arreglos*/
int matrix(int fila, int columna){
  return (n_dis*fila)+columna;
}
