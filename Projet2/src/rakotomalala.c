#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <omp.h>

#define INF -1

#define TAG_SIZES 11
#define TAG_SCATTER_ROWS 12
#define TAG_SCATTER_COLUMNS 13
#define TAG_CIRCULATE 14
#define TAG_GATHER 15


// Redéfinition du minimum et de l'addition pour l'algorithme de Floyd-Marshall
int min(int a, int b) {
    if (a < 0)
        return b;
    else if (b < 0)
        return a;
    else
        return a < b ? a : b;
}

int add(int a, int b) {
    if (a < 0 || b < 0)
        return -1;
    else
        return a + b;
}

// Structure utilisée pour le projet
struct Matrix {
    int** data;
    int columns;
    int rows;
};

/**
 * Affiche une matrice sur la sortie standard
 */
void printMatrix(struct Matrix *matrix) {
    for (int y = 0; y < matrix->rows; y++) {
        for (int x = 0; x < matrix->columns; x++) {
            if (matrix->data[y][x] == INF)
                printf("i ");
            else
                printf("%d ", matrix->data[y][x]);
        }
        printf("\n");
    }
}

/**
 * Permet d'allouer correctement une matrice
 */
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


/**
 * Permet de free une matrice
 */
void freeMatrix(struct Matrix *matrix) {
    #pragma omp parallel for
    for (int i = 0; i < matrix->rows; i++) {
        free(matrix->data[i]);
    }
    
    free(matrix->data);
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

/**
 * Effectue le produit matriciel entre une ligne et une colonne et stocke le résultat dans le bon indice de la matrice 'result'
 */
void product(struct Matrix* W_row, struct Matrix* W_column, struct Matrix* result, int nbr_tab, int startX) {
    int i = 0;
    for (int z = startX; z < (nbr_tab + startX); z++) {
        #pragma omp parallel for
        for (int y = 0; y < nbr_tab; y++) {
            result->data[y][z] = 0;
            #pragma omp parallel for
            for (int x = 0; x < W_row->columns; x++) {
                result->data[y][z] += (W_row->data[y][x] * W_column->data[i][x]);
            }
        }
        i++;
    }
}

/**
 * Applique l'alogorithme de Floyd-Marshall entre une ligne et une colonne et stocke le résultat dans le bon indice de la matrice 'result'
 */
void floyd(struct Matrix* W_row, struct Matrix* W_column, struct Matrix* result, int nbr_tab, int startX) {
    int i = 0;
    for (int z = startX; z < (nbr_tab + startX); z++) {
        #pragma omp parallel for
        for (int y = 0; y < nbr_tab; y++) {
            result->data[y][z] = INF;
            #pragma omp parallel for
            for (int x = 0; x < W_row->columns; x++) {
                result->data[y][z] = min(result->data[y][z], (add(W_row->data[y][x], W_column->data[i][x])));
            }
        }
        i++;
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

/**
 * Transforme la matrice lu en matrice adjacente W
 */
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
 * Permet à P0 d'envoyer un bout de matricesà P1
 */
void scatterInit(struct Matrix* W, int tab_size, int startY, int endY, int next, int tag) {
    for (int y = startY; y < endY; y++) {
        MPI_Send(W->data[y], tab_size, MPI_INT, next, tag, MPI_COMM_WORLD);
    }
}

/**
 * Permet à un processus de recevoir un bout de matrice de son prédécesseur, de récupérer le bout qui l'intéresse et d'envoyer le reste au suivant
 */
void scatter(struct Matrix* matrix, int previous, int next, int nbr_tab, int tab_size, int nbr_procs_used, int rank, int tag) {
    MPI_Status status;
    
    for (int y = 0; y < nbr_tab; y++) {
        MPI_Recv(matrix->data[y], tab_size, MPI_INT, previous, tag, MPI_COMM_WORLD, &status);
    }
    
    int tmp[tab_size];
    for (int i = rank; i < nbr_procs_used - 1; i++) {
        for (int y = 0; y < nbr_tab; y++) {
            MPI_Recv(&tmp, tab_size, MPI_INT, previous, tag, MPI_COMM_WORLD, &status);
            MPI_Send(&tmp, tab_size, MPI_INT, next, tag, MPI_COMM_WORLD);
        }
    }
}

/**
 * Permet à un processus d'envoyer la matrice résultat à son prédécesseur, et de recevoir du suivant les matrices résultats
 */
void gather(struct Matrix* result, int previous, int next, int nbr_tab, int tab_size, int nbr_procs_used, int rank) {
    MPI_Status status;
    
    for (int y = 0; y < nbr_tab; y++) {
        MPI_Send(result->data[y], tab_size, MPI_INT, previous, TAG_GATHER, MPI_COMM_WORLD);
    }
    
    int tmp[tab_size];
    for (int i = rank; i < nbr_procs_used - 1; i++) {
        for (int y = 0; y < nbr_tab; y++) {
            MPI_Recv(&tmp, tab_size, MPI_INT, next, TAG_GATHER, MPI_COMM_WORLD, &status);
            MPI_Send(&tmp, tab_size, MPI_INT, previous, TAG_GATHER, MPI_COMM_WORLD);
        }
    }
}

/**
 * Permet à P0 de récupérer toutes les lignes et des les placer au bon endroit dans sa matrice
 */
void gatherFinal(struct Matrix* result, int next, int nbr_tab, int tab_size, int nbr_procs_used) {
    MPI_Status status;
    
    for (int i = 0; i < nbr_procs_used - 1; i++) {
        for (int j = 0; j < nbr_tab; j++) {
            MPI_Recv(result->data[nbr_tab + (nbr_tab * i) + j], tab_size, MPI_INT, next, TAG_GATHER, MPI_COMM_WORLD, &status);
        }
    }
}

/**
 * Envoie une matrice au suivant et reçoit une matrice de son prédécesseur
 */
void circulate(struct Matrix* W_column, int nbr_tab, int tab_size, int next, int previous) {
    MPI_Status status;
    
    for (int i = 0; i < nbr_tab; i++) {
        MPI_Send(W_column->data[i], tab_size, MPI_INT, next, TAG_CIRCULATE, MPI_COMM_WORLD);
        MPI_Recv(W_column->data[i], tab_size, MPI_INT, previous, TAG_CIRCULATE, MPI_COMM_WORLD, &status);
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf("Fichier manquant en paramètre\n");
        exit(1);
    }
    
    int rank;
    int nbr_procs;
    int nbr_procs_used;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nbr_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int previous = ((rank - 1 + nbr_procs) % nbr_procs); // TODO: Maybe change with 'nbr_procs_used' when necessary ?
    int next = ((rank + 1) % nbr_procs); // TODO: Maybe change with 'nbr_procs_used' when necessary ?
    
    int tab_size;
    int nbr_tab;
    
    struct Matrix* W_row;
    struct Matrix* W_column;
    
    struct Matrix* result;
        
    if (rank == 0) {
        // Parse du fichier
        struct Matrix* A = parseFileAndFillMatrix(argv[1]);
        struct Matrix* W = transformToW(A);
        
        tab_size = W->rows;
        
        W_row = allocateMatrix(tab_size, nbr_tab);
        W_column = allocateMatrix(tab_size, nbr_tab);
        result = allocateMatrix(tab_size, tab_size);
        
        if ((tab_size / nbr_procs) < 1) {
            nbr_tab = 1;
            nbr_procs_used = tab_size;
        } else {
            nbr_tab = (int) (tab_size / nbr_procs);
            nbr_procs_used = nbr_procs;
        }
                
        // Envoie de la taille de la matrice à P1
        if (nbr_procs_used > 1) {
            MPI_Send(&tab_size, 1, MPI_INT, next, TAG_SIZES, MPI_COMM_WORLD);
            MPI_Send(&nbr_tab, 1, MPI_INT, next, TAG_SIZES, MPI_COMM_WORLD);
        }
        
        for (int i = 0; i < nbr_tab; i++) {
            W_row->data[i] = W->data[i];
            W_column->data[i] = transpose(W)->data[i];
        }
            
        // Scatter W en lignes et en colonnes
        for (int i = 1; i < nbr_procs_used; i++) {
            scatterInit(W, tab_size, (i * nbr_tab), ((i * nbr_tab) + nbr_tab), next, TAG_SCATTER_ROWS);
            scatterInit(transpose(W), tab_size, (i * nbr_tab), ((i * nbr_tab) + nbr_tab), next, TAG_SCATTER_COLUMNS);
        }

        for (int n = 0; n < tab_size - 1; n++) { // Matrice à la puissance N
            if (n != 0) {
                for (int y = 0; y < nbr_tab; y++) {
                    for (int x = 0; x < tab_size; x++)
                        W_row->data[y][x] = result->data[y][x];
                }
            }
            for (int i = nbr_procs_used; i > 0; i--) {
                floyd(W_row, W_column, result, nbr_tab, (nbr_tab * i) % tab_size);
                circulate(W_column, nbr_tab, tab_size, next, previous);
            }
        }
        
        // Récupération de tous les résultats
        gatherFinal(result, next, nbr_tab, tab_size, nbr_procs_used);

        // Affichage du résultat final
        printMatrix(result);
        
        freeMatrix(A);
        freeMatrix(W_row);
        freeMatrix(W_column);
        freeMatrix(result);
    } else {
        MPI_Status status;

        // Réception de la taille de la matrice du processus précedents
        MPI_Recv(&tab_size, 1, MPI_INT, previous, TAG_SIZES, MPI_COMM_WORLD, &status);
        MPI_Recv(&nbr_tab, 1, MPI_INT, previous, TAG_SIZES, MPI_COMM_WORLD, &status);
                
        // Allocation des matrices
        W_row = allocateMatrix(tab_size, nbr_tab);
        W_column = allocateMatrix(tab_size, nbr_tab);
        result = allocateMatrix(tab_size, nbr_tab);
        
        // Envoie de la taille de la matrice au processus suivant si le prochain n'est pas 0
        if (next != 0) {
            MPI_Send(&tab_size, 1, MPI_INT, next, TAG_SIZES, MPI_COMM_WORLD);
            MPI_Send(&nbr_tab, 1, MPI_INT, next, TAG_SIZES, MPI_COMM_WORLD);
        }
        
        if ((tab_size / nbr_procs) < 1)
            nbr_procs_used = tab_size;
        else
            nbr_procs_used = nbr_procs;
        
        scatter(W_row, previous, next, nbr_tab, tab_size, nbr_procs_used, rank, TAG_SCATTER_ROWS);
        scatter(W_column, previous, next, nbr_tab, tab_size, nbr_procs_used, rank, TAG_SCATTER_COLUMNS);
        
        for (int n = 0; n < tab_size - 1; n++) { // Matrice à la puissance N
            if (n != 0) {
                for (int y = 0; y < nbr_tab; y++) {
                    for (int x = 0; x < tab_size; x++)
                        W_row->data[y][x] = result->data[y][x];
                }
            }
            for (int i = nbr_procs_used; i > 0; i--) {
                int startX = ((nbr_tab * (i + rank)) % tab_size);
                floyd(W_row, W_column, result, nbr_tab, startX);
                circulate(W_column, nbr_tab, tab_size, next, previous);
            }
        }

        // Récupération des résultats des suivants et envoie des résultats au précédent
        gather(result, previous, next, nbr_tab, tab_size, nbr_procs_used, rank);
        
        freeMatrix(W_row);
        freeMatrix(W_column);
        freeMatrix(result);
    }
    
    MPI_Finalize();
    
    return 0;
}
