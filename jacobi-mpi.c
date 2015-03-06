// HW0_Problem2_Fang for HPC Spring15

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
/* timing in util.h requires -lrt flag to compile */
#include "util.h"
#define ff() printf("here is fine!\n");

double norm_res(double *u, int n, double *A); //calculate norm of residual |Au-f|
void Jacobi(double *u, int n, double *A);
void Gauss_Seidel(double *u, int n, double *A);

static int n_iter = 0;

int main(int argc, char **argv){

    if (argc != 2) {
        fprintf(stderr, "Functions need one input as number of discretization!\n");
        abort();
    }

    double *u;

    int n;
    n = atoi(argv[1]);
    double h = 1. / n;

    /* A[0] is diagonal of matrix A 
        A[1] is off-diagonal of matrix A */
    double A[] = {2. / pow(h,2),    -1. / pow(h,2)}; 

    u = (double *) calloc(n, sizeof(double)); // initialize to zeroes
    //printf("The u_i are: "); for (int i = 0; i < n; ++i) { printf("%f\t", u[i]); } printf("\n");
    //
    double ini_res, res;
    ini_res = norm_res(u, n, A);  //initial residual
    res = ini_res;
    
  timestamp_type time1, time2;
  get_timestamp(&time1);

  int iter_term = 0;
  static int step = 0;
  if (n == 1e3) {iter_term = 1e7 ; step = 1e5;}
  if (n == 1e5) {iter_term = 1e5; step = 1e3;}

    while ((res > ini_res * 1e-6) && (n_iter < iter_term)){ // terminate after iter_term iterations or res/ini_res < 1e-6
    Jacobi(u, n, A);
    if (n_iter % step == 0){ 
        res = norm_res(u, n, A);
        printf("After %d iterations, the norm of residule is %e, decreased %e \n", n_iter, res, res/ini_res);
    }
    }

    //while ((res > ini_res * 1e-6) && (n_iter < iter_term)){ // terminate after iter_term iterations or res/ini_res < 1e-6
    //Gauss_Seidel(u, n, A);
    //if (n_iter % step == 0){ 
    //    res = norm_res(u, n, A);
    //    printf("After %d iterations, the norm of residule is %e, decreased %e \n", n_iter, res, res/ini_res);
    //}
    //}

  get_timestamp(&time2);
  double elapsed = timestamp_diff_in_seconds(time1,time2);

  printf("Time elapsed is %f seconds.\n", elapsed);

    free(u);
    return 0;
}

double norm_res(double *u, int n, double *A){
    double *res;
    res = (double *) malloc(sizeof(double) * n);

    res[0] = A[0]*u[0] + A[1]*u[1]  - 1;
    for (int i = 1; i < n-1; i ++){
    res[i] =  A[1]*u[i-1] + A[0]*u[i] + A[1]*u[i+1] - 1;
    }
    res[n-1] = A[1]*u[n-2] + A[0]*u[n-1] - 1;
    
    // calculate the norm
    double n_res = 0;
    for (int i = 0; i < n; ++i){
        n_res += pow(res[i],2);
    }
    n_res = sqrt(n_res);

    //printf("The res are: "); for (int i = 0; i < n; ++i) { printf("%f\t", res[i]); } printf("\n");

//    printf("After %d iterations, the norm of residule is %e, decreased %f%% \n", n_iter, n_res, n_res/ini_res);

    free(res);
    return n_res;
}

void Jacobi(double *u, int n, double *A){
    ++n_iter;
    double old_pre = 0, old_cur;
    old_cur = u[0];
    u[0] = (1 -  A[1]*u[1]) / A[0];
    old_pre = old_cur;
    for (int i = 1; i < n-1; ++i) {
        old_cur = u[i];
        u[i] = (1 - (A[1]*old_pre + A[1]*u[i+1]))/A[0];
        old_pre = old_cur;
    }
        u[n-1] = (1 - A[1]*old_pre)/A[0];
    //printf("The u_i are: "); for (int i = 0; i < n; ++i) { printf("%f\t", u[i]); } printf("\n");
}

void Gauss_Seidel(double *u, int n, double *A){
    ++n_iter;
    u[0] = (1 -  A[1]*u[1]) / A[0];
    for (int i = 1; i < n-1; ++i) {
        u[i] = (1 - (A[1]*u[i-1]+ A[1]*u[i+1]))/A[0];
    }
        u[n-1] = (1 - A[1]*u[n-2])/A[0];
    //printf("The u_i are: "); for (int i = 0; i < n; ++i) { printf("%f\t", u[i]); } printf("\n");
}
