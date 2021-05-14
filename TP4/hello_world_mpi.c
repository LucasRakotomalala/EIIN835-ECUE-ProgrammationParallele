#include <stdio.h>
#include <mpi.h>

#include <time.h>
#include <stdlib.h>

void token_ring(int valeur, int emetteur, int rank, int numprocs);

int main(int argc, char *argv[])
{
    int rank, numprocs;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    /**
     Hello World
     */
    //printf("Hello World from %d\n", rank);
    
    /**
     Hello World 2
     */
    /*if (rank == 0)
        printf("Hello World (%d)\n", rank);
    else if (rank == numprocs -1)
        printf("Good bye (%d)\n", rank);
    else
        printf("... (%d)\n", rank);*/
    
    /**
     Échange de données
     */
    /*srand(time(NULL));
    
    int number = rand() % 10;
    if (rank == 0) {
        for (int i = 1; i < numprocs; i++) {
            MPI_Send(&number, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
        for (int i = 1; i < numprocs; i++) {
            MPI_Recv(&number, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("Process 0 received number %d from process %d\n", number, i);
        }

    } else {
        MPI_Recv(&number, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Process %d received number %d from process 0\n", rank, number);
        number *= 2;
        MPI_Send(&number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }*/
    
    /**
     Anneau de processus
     */
    /*int number = (rank + 1) % numprocs;
    MPI_Send(&rank, 1, MPI_INT, number, 0, MPI_COMM_WORLD);
    MPI_Recv(&number, 1, MPI_INT, (rank - 1) % numprocs, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("Mon rang est %d et j'ai reçu %d\n", rank, number);*/
    
    /**
     1. Si tous les process 'send' en même temps, il n'y personne qui 'receive' donc tout le process est bloqué
     2. Cela peut être résolu en ayant un processus non synchronisé avec les autres
     */
    
    srand(time(NULL));
    int number = rand() % 20;
    
    for (int i = 0; i < numprocs; i++) {
        token_ring(number, i, rank, numprocs);
    }

    MPI_Finalize();
        
    return 0;
}

void token_ring(int value, int emetteur, int rank, int numprocs) {
    if (emetteur == rank) {
        MPI_Send(&rank, 1, MPI_INT, (rank + 1) % numprocs, 0, MPI_COMM_WORLD);
        MPI_Recv(&value, 1, MPI_INT, (rank - 1) % numprocs, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Mon rang est %d et j'ai reçu %d\n", rank, value);
    } else {
        MPI_Recv(&value, 1, MPI_INT, (rank-  1) % numprocs, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Mon rang est %d et j'ai reçu %d\n", rank, value);
        MPI_Send(&rank, 1, MPI_INT, (rank + 1) % numprocs, 0, MPI_COMM_WORLD);
    }
    
    /*if (rank == emetteur) {
        int number = value;
        printf("Mon rang est %d et je tourne %d\n", rank, number);
        
        int number_received;
        MPI_Send(&number, 1, MPI_INT, (rank + 1) % numprocs, emetteur, MPI_COMM_WORLD);
        MPI_Recv(&number_received, 1, MPI_INT, (rank - 1) % numprocs, emetteur, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Mon rang est %d et j'ai reçu %d\n", rank, number);
    }
    else {
        int number;
        MPI_Recv(&number, 1, MPI_INT, (rank - 1) % numprocs, emetteur, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&number, 1, MPI_INT, (rank+1) % numprocs, emetteur, MPI_COMM_WORLD);
        printf("Mon rang est %d et j'envoie %d à %d\n", rank, number, (rank + 1) % numprocs);
    }*/
}
