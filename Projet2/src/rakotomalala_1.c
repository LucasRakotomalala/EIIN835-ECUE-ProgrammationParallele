#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <omp.h>

#define INF 999

#define TAG_SIZES 123
#define TAG_SCATTER 12
#define TAG_GATHER 16

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

struct Matrix* allocateMatrix(int size) {
    struct Matrix *tmp = malloc(sizeof(struct Matrix));
    
    tmp->columns = size;
    tmp->rows = size;
    tmp->data = (int**) malloc(sizeof(int*) * size);
        
    #pragma omp parallel for
    for (int i = 0; i < size; i++) {
        tmp->data[i] = (int*) malloc(sizeof(int) * size);
    }
    
    return tmp;
}

void freeMatrix(struct Matrix *matrix) {
    #pragma omp parallel for
    for (int i = 0; i < matrix->rows; i++) {
        free(matrix->data[i]);
    }
    
    free(matrix->data);
    free(matrix);
}

/**
 * Produit matriciel distribué, servant de base à floyd()
 */
void product(struct Matrix* A, struct Matrix* result, int startZ, int tab_size) {
    #pragma omp parallel for
    for (int y = 0; y < A->columns; y++) {
        result->data[startZ][y] = 0;
        for (int x = 0; x < A->rows; x++) {
            result->data[startZ][y] += (A->data[x][y] * A->data[startZ][x]);
        }
    }
}

/*void floyd(struct Matrix* W, struct Matrix* result, int startZ, int tab_size) {
    #pragma omp parallel for
    for (int i = 0; i < tab_size; i++)
        for (int j = 0; j < tab_size; j++)
            result->data[i][j] = INF;
 
    #pragma omp parallel for
    for (int y = 0; y < tab_size; y++) {
        for (int x = 0; x < tab_size; x++) {
            if (result->data[startZ][y] > (W->data[x][y] + W->data[startZ][x]))
                result->data[startZ][y] = W->data[x][y] + W->data[startZ][x];
        }
    }
}*/

void floyd(struct Matrix* W, struct Matrix* result, int startZ, int tab_size) {
    #pragma omp parallel for
    for (int y = 0; y < tab_size; y++) {
        result->data[startZ][y] = INF;
        #pragma omp parallel for
        for (int x = 0; x < tab_size; x++) {
            if (result->data[startZ][y] > (W->data[x][y] + W->data[startZ][x]))
                result->data[startZ][y] = (W->data[x][y] + W->data[startZ][x]);
        }
    }
}

/**
 * Séquentiel
 */
/*void floyd(struct Matrix* W, struct Matrix* result, int startZ, int tab_size) {
    
    for (int i = 0; i < tab_size; i++)
        for (int j = 0; j < tab_size; j++)
            result->data[i][j] = W->data[i][j];
    
    #pragma omp parallel for
    for (int k = 0; k < tab_size; k++) {
        for (int i = 0; i < tab_size; i++) {
            for (int j = 0; j < tab_size; j++) {
                if (W->data[i][k] + W->data[k][j] < result->data[i][j])
                    result->data[i][j] = W->data[i][k] + W->data[k][j];
            }
        }
    }
}*/

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
            
    struct Matrix* matrix = allocateMatrix(size); // Parce que nous sommes dans le cas des matrices carrées
    
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


/**
 * Calcule la transposée de la matrice passée en paramètre
 */
struct Matrix* transpose(struct Matrix* matrix) {
    struct Matrix *transpose = allocateMatrix(matrix->columns);

    #pragma omp parallel for
    for (int y = 0; y < matrix->rows; y++) {
        #pragma omp parallel for
        for (int x = 0; x < matrix->columns; x++) {
            transpose->data[x][y] = matrix->data[y][x];
        }
    }

    return transpose;
}

struct Matrix* transformToW(struct Matrix* A) {
    struct Matrix *W = allocateMatrix(A->columns);

    # pragma omp parallel for
    for (int y = 0; y < A->rows; y++) {
        # pragma omp parallel for
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

void scatter(struct Matrix* W, int tab_size, int rank, int nbr_procs, int next) {
    MPI_Status status;
    
    if (rank == 0) {
        //for (int i = 1; i < nbr_procs; i++) {
        for (int i = 0; i < nbr_procs; i++) {
            MPI_Send(W->data[i], tab_size, MPI_INT, next, TAG_SCATTER, MPI_COMM_WORLD); // On envoie les lignes de la matrice W à P1 (1 à tab_size)
        }
    } else {
        //for (int i = 0; i < (nbr_procs - rank); i++) {
        for (int i = 0; i < nbr_procs; i++) {
            MPI_Recv(W->data[i], tab_size, MPI_INT, ((rank - 1 + nbr_procs) % nbr_procs), TAG_SCATTER, MPI_COMM_WORLD, &status); // Réception des lignes du 'previous'
        }
        //for (int i = 1; i < (nbr_procs - rank); i++) { // On envoie que les lignes de 1 à tab_size
        for (int i = 0; i < nbr_procs; i++) {
            MPI_Send(W->data[i], tab_size, MPI_INT, next, TAG_SCATTER, MPI_COMM_WORLD);
        }
    }
}

void gather(struct Matrix* result, int tab_size, int rank, int nbr_procs, int next, int previous) {
    MPI_Status status;
    
    if (rank == 0) {
        for (int i = 1; i < nbr_procs; i++) {
            MPI_Recv(result->data[i], tab_size, MPI_INT, next, TAG_GATHER, MPI_COMM_WORLD, &status);
        }
    } else {
        MPI_Send(result->data[rank], tab_size, MPI_INT, previous, TAG_GATHER, MPI_COMM_WORLD);
        //MPI_Send(result->data[0], tab_size, MPI_INT, previous, TAG_GATHER, MPI_COMM_WORLD);
        for (int i = 1; i < (nbr_procs - rank); i++) {
            MPI_Recv(result->data[i], tab_size, MPI_INT, next, TAG_GATHER, MPI_COMM_WORLD, &status);
            MPI_Send(result->data[i], tab_size, MPI_INT, previous, TAG_GATHER, MPI_COMM_WORLD);
        }
    }
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
    
    struct Matrix* W;
    struct Matrix* result;
        
    if (rank == 0) {
        // Parse du fichier
        struct Matrix* A = parseFileAndFillMatrix(argv[1]);
        W = transformToW(A);
                
        tab_size = W->rows;
        matrix_size = tab_size * tab_size;
        
        result = allocateMatrix(tab_size);
        
        // Envoie de la taille de la matrice à P1
        MPI_Send(&matrix_size, 1, MPI_INT, next, TAG_SIZES, MPI_COMM_WORLD);
        MPI_Send(&tab_size, 1, MPI_INT, next, TAG_SIZES, MPI_COMM_WORLD);
        
        // Scatter A
        /*scatter(A, tab_size, rank, nbr_procs, next);
        product(A, result, rank,tab_size);
        // Gather result
        gather(result, tab_size, rank, nbr_procs, next, previous);*/
        //printMatrix(result);
        
        // Scatter W
        scatter(W, tab_size, rank, nbr_procs, next);
        
        floyd(W, result, rank, tab_size);
                        
        // Gather result
        gather(result, tab_size, rank, nbr_procs, next, previous);
        
        printMatrix(result);
        
        freeMatrix(A);
        freeMatrix(W);
        freeMatrix(result);
    } else {
        // Réception de la taille de la matrice du processus précedents
        MPI_Recv(&matrix_size, 1, MPI_INT, previous, TAG_SIZES, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&tab_size, 1, MPI_INT, previous, TAG_SIZES, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        // Allocation de la matrice
        W = allocateMatrix(tab_size);
        result = allocateMatrix(tab_size);
        
        // Envoie de la taille de la matrice au processus suivant si le prochain n'est pas 0
        if (rank != nbr_procs - 1) {
            MPI_Send(&matrix_size, 1, MPI_INT, next, TAG_SIZES, MPI_COMM_WORLD);
            MPI_Send(&tab_size, 1, MPI_INT, next, TAG_SIZES, MPI_COMM_WORLD);
        }
        
        // Scatter
        scatter(W, tab_size, rank, nbr_procs, next);
        
        floyd(W, result, rank, tab_size);
        
        //printf("\nRank %d:\n", rank);
        //printMatrix(result);
        
        freeMatrix(W);
        freeMatrix(result);
    }
    
    MPI_Finalize();
    
    return 0;
}
