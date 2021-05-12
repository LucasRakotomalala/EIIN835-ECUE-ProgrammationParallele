#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <omp.h>

#define INF 999

#define TAG_SIZES 123
#define TAG_SCATTER 12
#define TAG_GATHER 16
#define TAG_CIRCULATE 15

struct Matrix {
    int** data;
    int columns;
    int rows;
};

void printMatrix(struct Matrix *matrix) {
    for (int y = 0; y < matrix->rows; y++) {
        for (int x = 0; x < matrix->columns; x++) {
            if (matrix->data[y][x] >= INF)
                printf("i ");
            else
                printf("%d ", matrix->data[y][x]);
        }
        printf("\n");
    }
}

struct Matrix* allocateMatrix(int columns, int rows) {
    struct Matrix *tmp = malloc(sizeof(struct Matrix));
    
    tmp->columns = columns;
    tmp->rows = rows;
    tmp->data = (int**) malloc(sizeof(int*) * rows);
        
    #pragma omp parallel for
    for (int i = 0; i < rows; i++) {
        tmp->data[i] = (int*) malloc(sizeof(int) * columns);
    }
    
    return tmp;
}

void freeMatrix(struct Matrix *matrix) {
    #pragma omp parallel for
    for (int i = 0; i < matrix->rows; i++) {
        free(matrix->data[i]);
    }
    
    //free(matrix->data);
    free(matrix);
}

/**
 * Calcule la transposée de la matrice passée en paramètre
 */
struct Matrix* transpose(struct Matrix* matrix) {
    struct Matrix *transpose = allocateMatrix(matrix->rows, matrix->columns);

    #pragma omp parallel for
    for (int y = 0; y < matrix->rows; y++) {
        #pragma omp parallel for
        for (int x = 0; x < matrix->columns; x++) {
            transpose->data[x][y] = matrix->data[y][x];
        }
    }

    return transpose;
}

void product(struct Matrix* W_row, struct Matrix* W_column, struct Matrix* result, int startZ) {
    result->data[0][startZ] = 0;
    //#pragma omp parallel for // Ne marche pas en parallèle
    for (int x = 0; x < W_row->columns; x++) {
        result->data[0][startZ] += (W_row->data[0][x] * W_column->data[0][x]);
    }
}

void floyd(struct Matrix* W_row, struct Matrix* W_column, struct Matrix* result, int startZ) {
    result->data[0][startZ] = INF;
    //#pragma omp parallel for // Ne marche pas tout le temps
    for (int x = 0; x < W_row->columns; x++) {
        if (W_row->data[0][x] + W_column->data[0][x] < result->data[0][startZ])
            result->data[0][startZ] = W_row->data[0][x] + W_column->data[0][x];
    }
}

/**
 * Lecture du fichier pour construire la matrice
 */
struct Matrix* parseFileAndFillMatrix(char* filePath) {
    FILE* file;
    
    if ((file = fopen(filePath, "r")) == NULL) {
        printf("Erreur sur l'ouverture du fichier\n");
        exit(1);
    }
    
    int size = 1;
    char c;
    int nb;
        
    for (c = fgetc(file); c != '\n'; c = fgetc(file)) {
        if (c == ' ') {
            size++;
        }
    }
            
    struct Matrix* matrix = allocateMatrix(size, size); // Parce que nous sommes dans le cas des matrices carrées
    
    // Retour au début du fichier
    rewind(file);
    
    // Remplissage de la matrice
    fscanf(file, "%d", &nb);
      for (int y = 0; !feof(file); y++){
        for (int x = 0; x < size; x++) {
            matrix->data[y][x] = nb;
            fscanf(file, "%d", &nb);
        }
      }
    
    fclose(file);
    
    return matrix;
}

struct Matrix* transformToW(struct Matrix* A) {
    struct Matrix *W = allocateMatrix(A->columns, A->rows);

    #pragma omp parallel for
    for (int y = 0; y < A->rows; y++) {
        #pragma omp parallel for
        for (int x = 0; x < A->columns; x++) {
            if (x == y)
                W->data[x][y] = 0;
            else if (A->data[x][y] > 0)
                W->data[x][y] = A->data[x][y];
            else
                W->data[x][y] = INF;
        }
    }

    return W;
}

/**
 * Scatter une matrice et stocke la 1ère ligne dans le 2ème paramètre
 */
void scatter(struct Matrix* W, struct Matrix* matrix, int tab_size, int rank, int nbr_procs, int next) {
    matrix->data[0] = W->data[0];
    for (int i = nbr_procs - 1; i > 0; i--) {
        MPI_Send(W->data[i], tab_size, MPI_INT, next, TAG_SCATTER, MPI_COMM_WORLD);
    }
}

void receiveAndSend(struct Matrix* matrix, int tab_size, int rank, int nbr_procs, int next, int previous) {
    MPI_Status status;
    struct Matrix* tmp = allocateMatrix(tab_size, 1);
    
    for (int i = nbr_procs - rank - 1; i > 0; i--) {
        MPI_Recv(tmp->data[0], tab_size, MPI_INT, previous, TAG_SCATTER, MPI_COMM_WORLD, &status);
        MPI_Send(tmp->data[0], tab_size, MPI_INT, next, TAG_SCATTER, MPI_COMM_WORLD);
    }
    
    MPI_Recv(matrix->data[0], tab_size, MPI_INT, previous, TAG_SCATTER, MPI_COMM_WORLD, &status);
    
    freeMatrix(tmp);
}

void gather(struct Matrix* result, int tab_size, int rank, int nbr_procs, int next, int previous) {
    MPI_Status status;
    
    if (rank == 0) {
        for (int i = 1; i < nbr_procs; i++) {
            MPI_Recv(result->data[i], tab_size, MPI_INT, next, TAG_GATHER, MPI_COMM_WORLD, &status);
        }
    } else {
        //MPI_Send(result->data[rank], tab_size, MPI_INT, previous, TAG_GATHER, MPI_COMM_WORLD);
        MPI_Send(result->data[0], tab_size, MPI_INT, previous, TAG_GATHER, MPI_COMM_WORLD);
        for (int i = 1; i < (nbr_procs - rank); i++) {
            MPI_Recv(result->data[0], tab_size, MPI_INT, next, TAG_GATHER, MPI_COMM_WORLD, &status);
            MPI_Send(result->data[0], tab_size, MPI_INT, previous, TAG_GATHER, MPI_COMM_WORLD);
        }
    }
}

void circulate(struct Matrix* W_column, int tab_size, int next, int previous) {
    MPI_Status status;
    
    MPI_Send(W_column->data[0], tab_size, MPI_INT, next, TAG_CIRCULATE, MPI_COMM_WORLD);
    MPI_Recv(W_column->data[0], tab_size, MPI_INT, previous, TAG_CIRCULATE, MPI_COMM_WORLD, &status);
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf("Fichier manquant en paramètre\n");
        exit(1);
    }
    
    int rank;
    int nbr_procs;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nbr_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int previous = ((rank - 1 + nbr_procs) % nbr_procs);
    int next = ((rank + 1) % nbr_procs);
    
    int matrix_size;
    int tab_size;
    
    struct Matrix* W_row;
    struct Matrix* W_column;
    
    struct Matrix* result;
        
    if (rank == 0) {
        // Parse du fichier
        struct Matrix* A = parseFileAndFillMatrix(argv[1]);
        struct Matrix* W = transformToW(A);
        
        //printf("Matrice W:\n");
        //printMatrix(W);
        //printf("\n");
                
        tab_size = W->rows;
        matrix_size = tab_size * tab_size;
        
        W_row = allocateMatrix(tab_size, 1);
        W_column = allocateMatrix(tab_size, 1);
        
        result = allocateMatrix(tab_size, tab_size);

        // Envoie de la taille de la matrice à P1
        MPI_Send(&matrix_size, 1, MPI_INT, next, TAG_SIZES, MPI_COMM_WORLD);
        MPI_Send(&tab_size, 1, MPI_INT, next, TAG_SIZES, MPI_COMM_WORLD);
        
        // Scatter W en lignes et en colonnes
        scatter(W, W_row, tab_size, rank, nbr_procs, next);
        scatter(transpose(W), W_column, tab_size, rank, nbr_procs, next);

        for (int i = nbr_procs; i > 0; i--) {
            floyd(W_row, W_column, result, (rank + i) % nbr_procs);
            circulate(W_column, tab_size, next, previous);
        }
        
        gather(result, tab_size, rank, nbr_procs, next, previous);
        
        printMatrix(result);
        
        freeMatrix(A);
        freeMatrix(W_row);
        freeMatrix(W_column);
        freeMatrix(result);
    } else {
        // Réception de la taille de la matrice du processus précedents
        MPI_Recv(&matrix_size, 1, MPI_INT, previous, TAG_SIZES, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&tab_size, 1, MPI_INT, previous, TAG_SIZES, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        // Allocation de la matrice
        W_row = allocateMatrix(tab_size, 1);
        W_column = allocateMatrix(tab_size, 1);
        result = allocateMatrix(tab_size, 1);
        
        // Envoie de la taille de la matrice au processus suivant si le prochain n'est pas 0
        if (rank != nbr_procs - 1) {
            MPI_Send(&matrix_size, 1, MPI_INT, next, TAG_SIZES, MPI_COMM_WORLD);
            MPI_Send(&tab_size, 1, MPI_INT, next, TAG_SIZES, MPI_COMM_WORLD);
        }
        
        receiveAndSend(W_row, tab_size, rank, nbr_procs, next, previous);
        receiveAndSend(W_column, tab_size, rank, nbr_procs, next, previous);
        
        for (int i = 0; i < nbr_procs; i++) {
            floyd(W_row, W_column, result, (rank - i + nbr_procs) % nbr_procs);
            circulate(W_column, tab_size, next, previous);
        }
            
        gather(result, tab_size, rank, nbr_procs, next, previous);
        
        freeMatrix(W_row);
        freeMatrix(W_column);
        freeMatrix(result);
    }
    
    MPI_Finalize();
    
    return 0;
}
