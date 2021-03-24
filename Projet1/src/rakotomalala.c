#include <stdlib.h>
#include <stdio.h>
#include <tgmath.h> // fmaxl()
#include <limits.h> // Élément neutre du max des entiers long
#include <math.h> // pow(), log2()
#include <omp.h> // #pragma

struct tablo {
    long *tab;
    int size;
};

void printTablo(struct tablo *tmp) {
    printf("[%i] :", tmp->size);
    int size = tmp->size;
    
    for (int i = 0; i < size; ++i) {
        printf(" %ld", tmp->tab[i]);
    }
    
    printf("\n");
}

struct tablo *allocateTablo(int size) {
    struct tablo *tmp = malloc(sizeof(struct tablo));
    tmp->size = size;
    tmp->tab = malloc(size * sizeof(long));
    
    return tmp;
}

void freeTablo(struct tablo *tmp) {
    free(tmp->tab);
    free(tmp);
}

struct tablo *parseFileAndFillTablo(FILE *file) {
    struct tablo *source;
    
    int size = 0;
    long nb = 0;
    
    // Recherche du nombre d'entiers (long) présent dans le fichier pour éviter de faire des realloc
    fscanf(file, "%ld", &nb);
    while (!feof(file)) {
        fscanf(file, "%ld", &nb);
        size++;
    }
    
    source = allocateTablo(size);
    
    // Retour au début du fichier
    rewind(file);
    
    // Remplissage du tableau
    fscanf(file, "%ld", &nb);
    for (int i = 0; !feof(file); ++i) {
        source->tab[i] = nb;
        fscanf(file, "%ld", &nb);
    }
        
    return source;
}

void up(struct tablo *source, struct tablo *dest) {
    // Copie du tableau initial à la fin du nouveau tableau
    #pragma omp parallel for
    for (int i = 0; i < source->size; i++) {
        dest->tab[i + source->size] = source->tab[i];
    }
    
    // Algorithme de montée
    for (int i = log2(source->size) - 1; i >= 0; i--) {
        #pragma omp parallel for
        for (int j = pow(2, i); j <= (int) pow(2, i + 1) - 1; j++) {
            dest->tab[j] = dest->tab[2 * j] + dest->tab[2 * j + 1];
        }
    }
}

void down(struct tablo *a, struct tablo *b) {
    b->tab[1] = 0; // Élément neutre de la somme
    
    // Algorithme de descente préfixe
    for (int i = 1; i <= log2(a->size / 2); i++) {
        
        #pragma omp parallel for
        for (int j = pow(2, i); j <= (int) pow(2, i + 1) - 1; j++) {
            if (j % 2 == 0) // Si j est pair
                b->tab[j] = b->tab[j / 2];
            else
                b->tab[j] = b->tab[(j - 1) / 2] + a->tab[j - 1];
        }
    }
}

void downSuffix(struct tablo *a, struct tablo *b) {
    b->tab[1] = 0;
    
    // Algorithme de descente suffixe
    for (int i = 1; i <= log2(a->size / 2); i++) {
        
        #pragma omp parallel for
        for (int j = pow(2, i + 1) - 1; j >= (int) pow(2, i); j--) {
            if (j % 2 == 0) // Si j est pair
                b->tab[j] = b->tab[j / 2] + a->tab[j + 1];
            else
                b->tab[j] = b->tab[(j - 1) / 2];
        }
    }
}

void upMax(struct tablo *source, struct tablo *dest) {
    // Copie du tableau initial à la fin du nouveau tableau
    #pragma omp parallel for
    for (int i = 0; i < source->size; i++) {
        dest->tab[i + source->size] = source->tab[i];
    }
    
    // Algorithme de montée
    for (int i = log2(source->size) - 1; i >= 0; i--) {
        #pragma omp parallel for
        for (int j = pow(2, i); j <= (int) pow(2, i + 1) - 1; j++) {
            dest->tab[j] = fmaxl(dest->tab[2 * j], dest->tab[2 * j + 1]);
        }
    }
}

void downMax(struct tablo *a, struct tablo *b) {
    b->tab[1] = LONG_MIN; // Élément neutre du max des entiers long
    
    // Algorithme de descente préfixe max
    for (int i = 1; i <= log2(a->size / 2); i++) {
        #pragma omp parallel for
        for (int j = pow(2, i); j <= (int) pow(2, i + 1) - 1; j++) {
            if (j % 2 == 0) // Si j est pair
                b->tab[j] = b->tab[j / 2];
            else
                b->tab[j] = fmaxl(b->tab[(j - 1) / 2], a->tab[j - 1]);
        }
    }
}

void downMaxSuffix(struct tablo *a, struct tablo *b) {
    b->tab[1] = LONG_MIN; // Élément neutre du max des entiers long
    
    // Algorithme de descente suffixe max
    for (int i = 1; i <= log2(a->size / 2); i++) {
        #pragma omp parallel for
        for (int j = pow(2, i + 1) - 1; j >= (int) pow(2, i); j--) {
            if (j % 2 == 0) // Si j est pair
                b->tab[j] = fmaxl(b->tab[j / 2], a->tab[j + 1]);
            else
                b->tab[j] = b->tab[(j - 1) / 2];
        }
    }
}

void final(struct tablo *a, struct tablo *b) {
    #pragma omp parallel for
    for (int i = pow(2, log2(a->size / 2)); i < (int) pow(2, log2(a->size / 2) + 1); i++) {
        b->tab[i] = b->tab[i] + a->tab[i];
    }
}

void finalMax(struct tablo *a, struct tablo *b) {
    #pragma omp parallel for
    for (int i = pow(2, log2(a->size / 2)); i < (int) pow(2, log2(a->size / 2) + 1); i++) {
        b->tab[i] = fmaxl(b->tab[i], a->tab[i]);
    }
}

void prefixSum(struct tablo *source, struct tablo *dest) {
    struct tablo *a = allocateTablo(source->size * 2);
    
    up(source, a);
    
    struct tablo *b = allocateTablo(source->size * 2);
    
    down(a, b);
    final(a, b);
    
    // Copie du résultat dans le tablo 'dest'
    #pragma omp parallel for
    for (int i = 0; i < dest->size; i++) {
        dest->tab[i] = b->tab[i + source->size];
    }
    
    freeTablo(a);
    freeTablo(b);
}

void suffixSum(struct tablo *source, struct tablo *dest) {
    struct tablo *a = allocateTablo(source->size * 2);
    
    up(source, a);
    
    struct tablo *b = allocateTablo(source->size * 2);
    
    downSuffix(a, b);
    final(a, b);
    
    // Copie du résultat dans le tablo 'dest'
    #pragma omp parallel for
    for (int i = 0; i < dest->size; i++) {
        dest->tab[i] = b->tab[i + source->size];
    }
    
    freeTablo(a);
    freeTablo(b);
}

void suffixMax(struct tablo *source, struct tablo *dest) {
    struct tablo *a = allocateTablo(source->size * 2);
    
    upMax(source, a);
    
    struct tablo *b = allocateTablo(source->size * 2);
    
    downMaxSuffix(a, b);
    finalMax(a, b);
    
    // Copie du résultat dans le tablo 'dest'
    #pragma omp parallel for
    for (int i = 0; i < dest->size; i++) {
        dest->tab[i] = b->tab[i + source->size];
    }
    
    freeTablo(a);
    freeTablo(b);
}

void prefixMax(struct tablo *source, struct tablo *dest) {
    struct tablo *a = allocateTablo(source->size * 2);
    
    upMax(source, a);
    
    struct tablo *b = allocateTablo(source->size * 2);
    
    downMax(a, b);
    finalMax(a, b);
    
    // Copie du résultat dans le tablo 'dest'
    #pragma omp parallel for
    for (int i = 0; i < dest->size; i++) {
        dest->tab[i] = b->tab[i + source->size];
    }
    
    freeTablo(a);
    freeTablo(b);
}

void displayResult(struct tablo *M, struct tablo *source) {
    long max = 0;
    int index = 0;
    
    // On cherche la valeur maximale et son index
    for (int i = 0; i < source->size; ++i) {
        if (max < M->tab[i]) {
            max = M->tab[i];
            index = i;
        }
    }
    
    // On affiche le max en premier
    printf("%ld", max);
    // Tant que cette valeur est présente dans M, on écrit l'entier correspondant du tableau source
    while (index < M->size) {
        if (!(M->tab[index] == max)) {
            index--;
            break;
        }
        else {
            printf(" %ld", source->tab[index]);
            index++;
        }
    }
    printf("\n");
}

int main(int argc, char **argv) {
    if (argc < 2) {
        printf("Fichier manquant en paramètre\n");
        exit(1);
    }
    
    FILE *file = fopen(argv[1], "r");
    
    struct tablo *Q = parseFileAndFillTablo(file);
    
    fclose(file);
    
    struct tablo *PSUM = allocateTablo(Q->size);
    struct tablo *SSUM = allocateTablo(Q->size);
    struct tablo *SMAX = allocateTablo(Q->size);
    struct tablo *PMAX = allocateTablo(Q->size);
    
    prefixSum(Q, PSUM);
    suffixSum(Q, SSUM);
    suffixMax(PSUM, SMAX);
    prefixMax(SSUM, PMAX);
    
    /*printTablo(Q);
     printTablo(PSUM);
     printTablo(SSUM);
     printTablo(SMAX);
     printTablo(PMAX);*/
    
    struct tablo *M = allocateTablo(Q->size);
    
    // Étape 5
    #pragma omp parallel for
    for (int i = 0; i < Q->size; i++) {
        M->tab[i] = PMAX->tab[i] - SSUM->tab[i] + SMAX->tab[i] - PSUM->tab[i] + Q->tab[i];
    }
    
    //printTablo(M);
    
    displayResult(M, Q);
    
    freeTablo(Q);
    freeTablo(PSUM);
    freeTablo(SSUM);
    freeTablo(SMAX);
    freeTablo(PMAX);
    freeTablo(M);
    
    return 0;
}
