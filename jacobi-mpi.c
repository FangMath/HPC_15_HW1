// HPC_15_HW1 jacobi-mpi.c
#include<stdio.h>
#include<stdlib.h>
#include <unistd.h>
#include<math.h>
#include "util.h"
#include <mpi.h>
#define ff() printf("here is fine!\n");

void Jacobi_chunk(double *u, int n, double *A);

int main(int argc, char **argv){

    if (argc != 3) {
        fprintf(stderr, "Functions need one input as number of discretization!\n");
        abort();
    }

    int N, N_iter, n_iter;;
    N = atoi(argv[1]);
    double h = 1. / N;
    N_iter = atoi(argv[2]);

    int rank, tag, p;
    MPI_Status status;

    char hostname[1024];
    gethostname(hostname, 1024);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p); 

    if (N%p != 0) {
        fprintf(stderr, "Grid point N must be a multiple of processor number p!\n");
        abort();
    }
    
    int Np = N/p;

  double message_out_f, message_out_b;
  double message_in_f, message_in_b;

  tag = 99;

    double *u;
    u = (double *) calloc(Np+2, sizeof(double)); // initialize to zeroes, global ghost nodes of U[-1]=0 and U[N+1]=0

    /* A[0] is diagonal of matrix A 
        A[1] is off-diagonal of matrix A */
    double A[] = {2. / pow(h,2),    -1. / pow(h,2)}; 
    //printf("The u_i are: "); for (int i = 0; i < n; ++i) { printf("%f\t", u[i]); } printf("\n");
    
  timestamp_type time1, time2; //time it
  get_timestamp(&time1);

    for (n_iter = 0; n_iter < N_iter; ++n_iter){ // terminate after N_iter iterations 
        // do computations
        Jacobi_chunk(u, Np+2, A);
        // communications
        if (rank == 0)
        {
            message_out_f = u[Np];
            MPI_Send(&message_out_f, 1, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD);
            MPI_Recv(&message_in_b, 1, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &status);
            u[Np+1] = message_in_b;
        }
        else if (rank == p-1)
        {
            message_out_b = u[1];
            MPI_Send(&message_out_b, 1, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD);
            MPI_Recv(&message_in_f, 1, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD, &status);
            u[0] = message_in_f;
        }
        else
        {
            message_out_f = u[Np];
            MPI_Send(&message_out_f, 1, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD);
            MPI_Send(&message_out_b, 1, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD);
            MPI_Recv(&message_in_b, 1, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &status);
            u[Np+1] = message_in_b;

            message_out_b = u[1];
            MPI_Recv(&message_in_f, 1, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD, &status);
            u[0] = message_in_f;
        }
    }

//    printf("The %dth communication, rank %d received from %d the message %d\n", nN, rank, origin, (int)*message_in);
  MPI_Finalize();

//  printf("Rank %d hosted on %s runs time %f seconds.\n", rank, hostname, elapsed);

  get_timestamp(&time2);
  double elapsed = timestamp_diff_in_seconds(time1,time2);

  printf("Time elapsed is %f seconds.\n", elapsed);

    free(u);
    return 0;
}

void Jacobi_chunk(double *u, int n, double *A){
    double old_pre = 0, old_cur;
    int i;
    old_pre = u[0];
    for (i = 1; i < n-1; ++i) {
        old_cur = u[i];
        u[i] = (1 - (A[1]*old_pre + A[1]*u[i+1]))/A[0];
        old_pre = old_cur;
    }
    //printf("The u_i are: "); for (int i = 0; i < n; ++i) { printf("%f\t", u[i]); } printf("\n");
}

