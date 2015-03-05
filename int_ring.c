/* Communication ring:
 * Exchange between messages between mpirank
 * 0 -> 1 -> 2 -> 3, ..., p -> 0
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "util.h"
#include <mpi.h>

int main( int argc, char *argv[])
{

    if (argc != 2) {
        fprintf(stderr, "Functions need one input as number of discretization!\n");
        abort();
    }

    int N, nN; // frequency of communications
    N = atoi(argv[1]);

  int rank, n, tag, origin, destination, num_proc;
  MPI_Status status;

  char hostname[1024];
  gethostname(hostname, 1024);

  timestamp_type time1, time2; //time it
  get_timestamp(&time1);

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_proc); 

  int msgN = 1;
  double *message_out = calloc(msgN, sizeof(double));
  double *message_in = calloc(msgN, sizeof(double));

  tag = 99;

    if(rank == 0) //first processor
    {
    destination = rank + 1;
    origin = num_proc-1;
    }
    else if (rank == num_proc-1) //last processor
    {
    destination = 0;
    origin = rank - 1;
    }
    else
    {
    destination = rank + 1;
    origin = rank - 1;
    }



    for (nN = 0; nN < N; ++nN)
    {
        if (rank == 0)
        {
            if (nN ==0) //head processor
            {
                for(n = 0; n < msgN; ++n) message_out[n] = 0;
            }
            else
            {
                for(n = 0; n < msgN; ++n) message_out[n] = message_in[n] + rank;
            }

            MPI_Send(message_out, msgN, MPI_DOUBLE, destination, tag, MPI_COMM_WORLD);
            MPI_Recv(message_in,  msgN, MPI_DOUBLE, origin, tag, MPI_COMM_WORLD, &status);
        }

        else //other processors
        {

            MPI_Recv(message_in,  msgN, MPI_DOUBLE, origin, tag, MPI_COMM_WORLD, &status);
            for(n = 0; n < msgN; ++n) message_out[n] = message_in[n] + rank;
            MPI_Send(message_out, msgN, MPI_DOUBLE, destination, tag, MPI_COMM_WORLD);
        }

//      printf("rank %d sent to %d the message %d\n", rank, destination, (int)*message_out);
//      printf("rank %d received from %d the message %d\n", rank, origin, (int)*message_in);
    }



  /*
   * Write output to a file
   */
  if (msgN == 1)
  {
//     printf("rank %d sent to %d the message %d\n", rank, destination, (int)*message_out);
      printf("rank %d received from %d the message %d\n", rank, origin, (int)*message_in);

      //printf("rank %d hosted on %s sent to %d the message %d\n", rank, hostname, destination, (int)*message_out);
      //printf("rank %d hosted on %s received from %d the message %d\n", rank, hostname, origin, (int)*message_in);
  }
  else
  {
    FILE* fd = NULL;
    char filename[256];
    snprintf(filename, 256, "output%02d.txt", rank);
    fd = fopen(filename,"w+");

    if(NULL == fd)
    {
      printf("Error opening file \n");
      return 1;
    }

    fprintf(fd, "rank %d received from %d the message:\n", rank, origin);
    for(n = 0; n < msgN; ++n)
      fprintf(fd, "  %d\n", (int)message_in[n]);

    fclose(fd);

  }

  MPI_Finalize();

  get_timestamp(&time2);
  double elapsed = timestamp_diff_in_seconds(time1,time2);
  printf("Time elapsed is %f seconds.\n", elapsed);

  free(message_out);
  free(message_in);
  return 0;
}
