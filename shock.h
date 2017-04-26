
#define GAMMA  1.4

#define rho1 1.0
#define P1 1.0
#define U1 0.0

#define P5 0.1
#define rho5 0.125
#define U5 0.0
#define x0 0.5
#define L 1.0 
#define time 0.2379915616352639

#define n_dis 5001

void init_array( double *U, double *UB,double *UA);
void calc_F( double *U, double *F);
void calc_star(double *U,double *UB,double *F,double *FB, double dt, double dx);
void calc_UA(double *U, double *UA,double *UB,double *FB, double dt, double dx);
double calc_umax(double *U, double *c);
void lax(double *U, double *F,double *UA,double *UB, double *FB, double *c, double dt, double dx);
void calc_var(double *U, double *rho, double *u, double *P);
int matrix(int fila, int columna);


