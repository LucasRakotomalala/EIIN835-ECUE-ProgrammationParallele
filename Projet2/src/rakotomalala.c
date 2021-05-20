#include <stdlib.h>
#include <stdio.h>
#include <mpi.h> // MPI
#include <omp.h> // #pragma

// Définitions de macros utilisées tout au long du projet
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

int sum(int a, int b) {
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
    
    free(matrix->data);
    free(matrix);
}

/**
 * Calcule la transposée de la matrice passée en paramètre
 * @param matrix  : matrice dont on calcule la transposée
 * @return la transposée de "matrix"
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
 * Applique l'algorithme de Floyd-Marshall entre une ligne et une colonne (voire plus) et stocke le résultat au bon indice de la matrice "result"
 * @param W_row : la matrice colonne possiblement élevée une puissance quelconque
 * @param W_column : la matrice colonne utilisée pour finalement élever "W_row" à la puissance N
 * @param result : la matrice dans laquelle on va stocker les résultat
 * @param nbr_tab : le nombre de ligne de la matrice
 * @param startX : l'indice à partir duquel on commence l'algorithme pour stocker au bon endroit le résultat
 * @return void
 */
void floyd(struct Matrix* W_row, struct Matrix* W_column, struct Matrix* result, int nbr_tab, int startX) {
    int i = 0;
    for (int z = startX; z < (nbr_tab + startX); z++) {
        #pragma omp parallel for
        for (int y = 0; y < nbr_tab; y++) {
            result->data[y][z] = INF;
            #pragma omp parallel for
            for (int x = 0; x < W_row->columns; x++) {
                result->data[y][z] = min(result->data[y][z], (sum(W_row->data[y][x], W_column->data[i][x])));
            }
        }
        i++;
    }
}

/**
 * Ouverture d'un fichier et lecture de celui-ci pour construire la matrice A
 * @param filePath : le chemin du fichier à ouvrir et lire
 * @return la matrice du fichier
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
            
    struct Matrix* matrix = allocateMatrix(size, size); // Le fichier d'entrée décrit une matrice carrée
    
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
 * Construit la matrice adjacente W de celle passée en paramètre
 * @param A : la matrice dont on doit construire sa matrice adjacente
 * @return la matrice adjacente de celle passée en paramètre
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
 * Permet à un processeur de recevoir de son prédécesseur et d'envoyer à son successeur un entier (si le successeur n'est pas 0)
 * @param value_to_broadcast : le pointeur de la valeur que le souhaite envoyer/récupérer
 * @param rank : le rang du processeur qui appelle la méthode
 * @param previous : le successeur du processeur actuel
 * @param next : le suivant du processeur actuel
 * @return void
 */
void broadcast(int* value_to_broadcast, int rank, int previous, int next) {
    MPI_Status status;

    if (rank != 0)
        MPI_Recv(value_to_broadcast, 1, MPI_INT, previous, TAG_SIZES, MPI_COMM_WORLD, &status);
    if (next != 0)
        MPI_Send(value_to_broadcast, 1, MPI_INT, next, TAG_SIZES, MPI_COMM_WORLD);
}

/**
 * Permet à un processeur de recevoir un bout de matrice de son prédécesseur, de récupérer le bout qui l'intéresse et d'envoyer le reste au suivant
 * @param matrix : la matrice dans laquelle on stocke certaines données
 * @param previous : le prédécesseur du processeur actuel
 * @param next : le successeur du processeur actuel
 * @param nbr_tab : le nombre de ligne de la matrice
 * @param tab_size : le nombre d'éléments par ligne
 * @param nbr_procs_used : le nombre de processeurs utilisés rééllement par le programme
 * @param rank : le rang du processeur qui appelle la méthode
 * @param tag : le tag MPI sur lequel on souhaite envoyer/récupérer les données
 * @return void
 */
void scatter(struct Matrix* matrix, int previous, int next, int nbr_tab, int tab_size, int nbr_procs_used, int rank, int tag) {
    if (rank == 0) {
        for (int i = 1; i < nbr_procs_used; i++) {
            for (int y = (i * nbr_tab); y < ((i * nbr_tab) + nbr_tab); y++) {
                MPI_Send(matrix->data[y], tab_size, MPI_INT, next, tag, MPI_COMM_WORLD);
            }
        }
    } else {
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
}

/**
 * Permet à un processeur (autre que P0) d'envoyer la matrice résultat à son successeur, et de recevoir du prédécesseur la matrice "résultat"
 * @param result : la matrice dans laquelle il y a les résultats de la matrice W élevée à la puissance N
 * @param previous : le prédécesseur du processeur actuel
 * @param next : le successeur du processeur actuel
 * @param nbr_tab : le nombre de ligne de la matrice
 * @param tab_size : le nombre d'éléments par ligne
 * @param nbr_procs_used : le nombre de processeurs utilisés rééllement par le programme
 * @param rank : le rang du processeur qui appelle la méthode
 * @return void
 */
void gather(struct Matrix* result, int previous, int next, int nbr_tab, int tab_size, int rank) {
    MPI_Status status;
    
    for (int y = nbr_tab - 1; y >= 0; y--) {
        MPI_Send(result->data[y], tab_size, MPI_INT, next, TAG_GATHER, MPI_COMM_WORLD);
    }
    
    int tmp[tab_size];
    for (int i = 0; i < rank - 1; i++) {
        for (int y = 0; y < nbr_tab; y++) {
            MPI_Recv(&tmp, tab_size, MPI_INT, previous, TAG_GATHER, MPI_COMM_WORLD, &status);
            MPI_Send(&tmp, tab_size, MPI_INT, next, TAG_GATHER, MPI_COMM_WORLD);
        }
    }
}

/**
 * Permet à P0 de récupérer toutes les lignes et des les placer au bon endroit dans sa matrice
 * @param result : la matrice dans laquelle on va stocker le résultat des matrices résultats des autres processeurs
 * @param previous : le prédécesseur du processeur actuel
 * @param nbr_tab : le nombre de ligne de la matrice
 * @param tab_size : le nombre d'éléments par ligne
 * @param nbr_procs_used : le nombre de processeurs utilisés rééllement par le programme
 * @return void
 */
void gatherFinal(struct Matrix* result, int previous, int nbr_tab, int tab_size, int nbr_procs_used) {
    MPI_Status status;
    
    for (int i = nbr_procs_used - 1; i > 0 ; i--) {
        for (int j = 0; j < nbr_tab; j++) {
            MPI_Recv(result->data[nbr_tab + (nbr_tab * i) - j - 1], tab_size, MPI_INT, previous, TAG_GATHER, MPI_COMM_WORLD, &status);
        }
    }
}

/**
 * Envoie une matrice au successeur et reçoit une matrice de son prédécesseur
 * @param W_column : la matrice à envoyer/recevoir
 * @param nbr_tab : le nombre de ligne de la matrice
 * @param tab_size : le nombre d'éléments par ligne
 * @param previous : le prédeccesseur du processeur actuel
 * @param next : le successeur du processeur actuel
 * @return void
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
    
    int previous = ((rank - 1 + nbr_procs) % nbr_procs);
    int next = ((rank + 1) % nbr_procs);
    
    int tab_size;
    int nbr_tab;
    
    struct Matrix* W_row;
    struct Matrix* W_column;
    
    struct Matrix* result;
        
    if (rank == 0) {
        // Parse du fichier
        struct Matrix* A = parseFileAndFillMatrix(argv[1]);
        
        // Transformation de la matrice A en matrice adjacente W
        struct Matrix* W = transformToW(A);
        
        // Définition des variables
        tab_size = W->rows;
        
        if ((tab_size / nbr_procs) < 1) {
            nbr_tab = 1;
            nbr_procs_used = tab_size;
        } else {
            nbr_tab = (int) (tab_size / nbr_procs);
            nbr_procs_used = nbr_procs;
        }
        
        // Allocation mémoire des matrices
        W_row = allocateMatrix(tab_size, nbr_tab);
        W_column = allocateMatrix(tab_size, nbr_tab);
        result = allocateMatrix(tab_size, tab_size);
        
        // Broadcast sur anneau du nombre de ligne(s)/colonne(s) à traiter par chaque processeur et de la taille d'une ligne/colonne
        broadcast(&tab_size, rank, previous, next);
        broadcast(&nbr_tab, rank, previous, next);
        
        // Création des matrices lignes et colonnes sur lesquelles on va travailler pour éviter d'utiliser la matrice W
        #pragma omp parallel for
        for (int i = 0; i < nbr_tab; i++) {
            W_row->data[i] = W->data[i];
            W_column->data[i] = transpose(W)->data[i];
        }
            
        // Scatter W en lignes et en colonnes
        scatter(W, previous, next, nbr_tab, tab_size, nbr_procs_used, rank, TAG_SCATTER_ROWS);
        scatter(transpose(W), previous, next, nbr_tab, tab_size, nbr_procs_used, rank, TAG_SCATTER_COLUMNS);

        // On élève la matrice ligne (W_row) à la puissance N
        for (int n = 0; n < tab_size - 1; n++) {
            if (n != 0) {
                #pragma omp parallel for
                for (int y = 0; y < nbr_tab; y++) {
                    #pragma omp parallel for
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
        gatherFinal(result, previous, nbr_tab, tab_size, nbr_procs_used);

        // Affichage du résultat final
        printMatrix(result);
        
        freeMatrix(A);
        freeMatrix(W_row);
        freeMatrix(W_column);
        freeMatrix(result);
    } else {
        // Broadcast sur anneau du nombre de ligne(s)/colonne(s) à traiter par chaque processeur et de la taille d'une ligne/colonne
        broadcast(&tab_size, rank, previous, next);
        broadcast(&nbr_tab, rank, previous, next);
        
        // Définition des variables
        if ((tab_size / nbr_procs) < 1) {
            nbr_procs_used = tab_size;
        }
        else {
            nbr_procs_used = nbr_procs;
        }

        // Allocation des matrices
        W_row = allocateMatrix(tab_size, nbr_tab);
        W_column = allocateMatrix(tab_size, nbr_tab);
        result = allocateMatrix(tab_size, nbr_tab);
        
        // Scatter W_row et W_column
        scatter(W_row, previous, next, nbr_tab, tab_size, nbr_procs_used, rank, TAG_SCATTER_ROWS);
        scatter(W_column, previous, next, nbr_tab, tab_size, nbr_procs_used, rank, TAG_SCATTER_COLUMNS);
        
        // On élève la matrice ligne (W_row) à la puissance N
        for (int n = 0; n < tab_size - 1; n++) {
            if (n != 0) {
                #pragma omp parallel for
                for (int y = 0; y < nbr_tab; y++) {
                    #pragma omp parallel for
                    for (int x = 0; x < tab_size; x++)
                        W_row->data[y][x] = result->data[y][x];
                }
            }
            for (int i = nbr_procs_used; i > 0; i--) {
                floyd(W_row, W_column, result, nbr_tab, ((nbr_tab * (i + rank)) % tab_size));
                circulate(W_column, nbr_tab, tab_size, next, previous);
            }
        }

        // Récupération des résultats des suivants et envoie des résultats au précédent
        gather(result, previous, next, nbr_tab, tab_size, rank);
        
        freeMatrix(W_row);
        freeMatrix(W_column);
        freeMatrix(result);
    }
    
    MPI_Finalize();
    
    return 0;
}
