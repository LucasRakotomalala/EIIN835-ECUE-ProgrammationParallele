#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <mpi.h>

int main() {
        srand(time(NULL));
        int rank, numprocs;
        MPI_Init(NULL, NULL);
		
        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		
        if (rank == 0) {
                int number_send = rand() % 100;
                for (size_t i = 1; i < numprocs; i++) {
                        MPI_Send(&number_send,
                                 1,
                                 MPI_INT,
                                 i,
                                 0,
                                 MPI_COMM_WORLD);
                        printf("number_send : %d to id : %d\n", number_send,i);
                }
                for (size_t i = 1; i < numprocs; i++) {
                        int number_recv_procs0 = malloc(1*sizeof(int));
                        MPI_Recv(&number_recv_procs0,
                                 1,
                                 MPI_INT,
                                 i,
                                 0,
                                 MPI_COMM_WORLD,
                                 MPI_STATUS_IGNORE);
                        printf("number_recv_procs0 : %d from %d\n", number_recv_procs0,i);
                }
        }
        else {
                int number_recv = malloc(1*sizeof(int));
                MPI_Recv(&number_recv,
                         1,
                         MPI_INT,
                         MPI_ANY_SOURCE,
                         0,
                         MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
                printf("rank : %d, number_recv : %d \n", rank,number_recv);
                int number_send_proc_1 = number_recv*2;
                MPI_Send(&number_send_proc_1,
                         1,
                         MPI_INT,
                         0,
                         0,
                         MPI_COMM_WORLD);
        }
        MPI_Finalize();
        return 0;
}