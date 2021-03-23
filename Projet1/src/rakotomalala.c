#include <stdlib.h>
#include <stdio.h>
#include <tgmath.h> // fmaxl()
#include <limits.h> // Élément neutre du max des entiers long
#include <math.h> // pow(), log()
#include <omp.h> // #pragma

struct tablo {
    long *tab;
    int size;
};

void printTablo(struct tablo *tmp) {
    printf("[%i] : ", tmp->size);
    int size = tmp->size;
    
    for (int i = 0; i < size; ++i)
        printf("%ld ", tmp->tab[i]);
    
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

struct tablo *parseFile(char *filePath) {
    struct tablo *source;
    FILE *file = fopen(filePath, "r");
    
    if (file == NULL) {
        printf("Impossible d'ouvrir le fichier\n");
        exit(1);
    }
    
    int size = 0;
    long nb = 0;
    
    // Trouve le nombre d'entier présent dans le fichier
    fscanf(file, "%ld", &nb);
    while (!feof (file)) {
        fscanf(file, "%ld", &nb);
        size += 1;
    }
        
    source = allocateTablo(size);
    
    rewind(file);
    
    // Remplit le tableau
    fscanf(file, "%ld", &nb);
    for (int i = 0; !feof (file); ++i) {
        source->tab[i] = nb;
        fscanf(file, "%ld", &nb);
    }
        
    fclose(file);
    
    return source;
}

void montee(struct tablo *source, struct tablo *dest) {
    // Copie du tableau initial à la fin du nouveau tableau
    #pragma omp parallel for
    for (int i = 0; i < source->size; i++) {
        dest->tab[i + source->size] = source->tab[i];
    }
    
    // Algorithme de montée
    for (int i = (int) log2(source->size) - 1; i >= 0; i--) {
        int beginIndex = pow(2, i);
        int endIndex = pow(2, i + 1) - 1;
        
        #pragma omp parallel for
        for (int j = beginIndex; j <= endIndex; j++) {
            dest->tab[j] = dest->tab[2 * j] + dest->tab[2 * j + 1];
        }
    }
}

void descente(struct tablo *a, struct tablo *b) {
    b->tab[1] = 0; // Élément neutre de la somme
    
    int m = (int) log2(a->size / 2);
    
    // Algorithme de descente préfixe
    for (int i = 1; i <= m; i++) {
        int beginIndex = pow(2, i);
        int endIndex = pow(2, i + 1) - 1;
        
        #pragma omp parallel for
        for (int j = beginIndex; j <= endIndex; j++) {
            if (j % 2 == 0) // Si j est pair
                b->tab[j] = b->tab[j / 2];
            else
                b->tab[j] = b->tab[(j - 1) / 2] + a->tab[j - 1];
        }
    }
}

void descenteSuffixe(struct tablo *a, struct tablo *b) {
    b->tab[1] = 0;
    
    int m = (int) log2(a->size / 2);
    
    // Algorithme de descente suffixe
    for (int i = 1; i <= m; i++) {
        int beginIndex = pow(2, i + 1) - 1;
        int endIndex = pow(2, i);
        
        #pragma omp parallel for
        for (int j = beginIndex; j >= endIndex; j--) {
            if (j % 2 == 0) // Si j est pair
                b->tab[j] = b->tab[j / 2] + a->tab[j + 1];
            else
                b->tab[j] = b->tab[(j - 1) / 2];
        }
    }
}

void monteeMax(struct tablo *source, struct tablo *dest) {
    // Copie du tableau initial à la fin du nouveau tableau
    #pragma omp parallel for
    for (int i = 0; i < source->size; i++) {
        dest->tab[i + source->size] = source->tab[i];
    }
    
    // Algorithme de montée
    for (int i = (int) log2(source->size) - 1; i >= 0; i--) {
        int beginIndex = pow(2, i);
        int endIndex = pow(2, i + 1) - 1;
        
        #pragma omp parallel for
        for (int j = beginIndex; j <= endIndex; j++) {
            dest->tab[j] = fmaxl(dest->tab[2 * j], dest->tab[2 * j + 1]);
        }
    }
}

void descenteMax(struct tablo *a, struct tablo *b) {
    b->tab[1] = LONG_MIN; // Élément neutre du max des entiers long
    
    int m = (int) log2(a->size / 2);
    
    // Algorithme de descente préfixe max
    for (int i = 1; i <= m; i++) {
        int beginIndex = pow(2, i);
        int endIndex = pow(2, i + 1) - 1;
        
        #pragma omp parallel for
        for (int j = beginIndex; j <= endIndex; j++) {
            if (j % 2 == 0) // Si j est pair
                b->tab[j] = b->tab[j / 2];
            else
                b->tab[j] = fmaxl(b->tab[(j - 1) / 2], a->tab[j - 1]);
        }
    }
}

void descenteSuffixeMax(struct tablo *a, struct tablo *b) {
    b->tab[1] = LONG_MIN; // Élément neutre du max des entiers long
    
    int m = (int) log2(a->size / 2);
    
    // Algorithme de descente suffixe max
    for (int i = 1; i <= m; i++) {
        int beginIndex = pow(2, i + 1) - 1;
        int endIndex = pow(2, i);
        
        #pragma omp parallel for
        for (int j = beginIndex; j >= endIndex; j--) {
            if (j % 2 == 0) // Si j est pair
                b->tab[j] = fmaxl(b->tab[j / 2], a->tab[j + 1]);
            else
                b->tab[j] = b->tab[(j - 1) / 2];
        }
    }
}

void final(struct tablo *a, struct tablo *b) {
    int beginIndex = pow(2, log2(a->size / 2));
    int endIndex = pow(2, log2(a->size / 2) + 1);
    
    #pragma omp parallel for
    for (int i = beginIndex; i < endIndex; i++) {
        b->tab[i] = b->tab[i] + a->tab[i];
    }
}

void finalMax(struct tablo *a, struct tablo *b) {
    int beginIndex = pow(2, log2(a->size / 2));
    int endIndex = pow(2, log2(a->size / 2) + 1);
    
    #pragma omp parallel for
    for (int i = beginIndex; i < endIndex; i++) {
        b->tab[i] = fmaxl(b->tab[i], a->tab[i]);
    }
}

void prefixSum(struct tablo *source, struct tablo *dest) {
    struct tablo *a = allocateTablo(source->size * 2);
    
    montee(source, a);
    
    struct tablo *b = allocateTablo(source->size * 2);
    
    descente(a, b);
    final(a, b);
    
    #pragma omp parallel for
    for (int i = 0; i < dest->size; i++) {
        dest->tab[i] = b->tab[i + source->size]; // Copie du résultat dans le tablo 'dest'
    }
    
    freeTablo(a);
    freeTablo(b);
}

void suffixSum(struct tablo *source, struct tablo *dest) {
    struct tablo *a = allocateTablo(source->size * 2);
    
    montee(source, a);
    
    struct tablo *b = allocateTablo(source->size * 2);
    
    descenteSuffixe(a, b);
    final(a, b);
    
    #pragma omp parallel for
    for (int i = 0; i < dest->size; i++) {
        dest->tab[i] = b->tab[i + source->size]; // Copie du résultat dans le tablo 'dest'
    }
    
    freeTablo(a);
    freeTablo(b);
}

void suffixMax(struct tablo *source, struct tablo *dest) {
    struct tablo *a = allocateTablo(source->size * 2);
    
    monteeMax(source, a);
    
    struct tablo *b = allocateTablo(source->size * 2);
    
    descenteSuffixeMax(a, b);
    finalMax(a, b);
    
    #pragma omp parallel for
    for (int i = 0; i < dest->size; i++) {
        dest->tab[i] = b->tab[i + source->size]; // Copie du résultat dans le tablo 'dest'
    }
    
    freeTablo(a);
    freeTablo(b);
}

void prefixMax(struct tablo *source, struct tablo *dest) {
    struct tablo *a = allocateTablo(source->size * 2);
    
    monteeMax(source, a);
    
    struct tablo *b = allocateTablo(source->size * 2);
    
    descenteMax(a, b);
    finalMax(a, b);
    
    #pragma omp parallel for
    for (int i = 0; i < dest->size; i++) {
        dest->tab[i] = b->tab[i + source->size]; // Copie du résultat dans le tablo 'dest'
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
    struct tablo *Q = parseFile(argv[1]);
        
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
    
    #pragma omp parallel for
    for (int i = 0; i < Q->size; i++) {
        long Ms = PMAX->tab[i] - SSUM->tab[i] + Q->tab[i];
        long Mp = SMAX->tab[i] - PSUM->tab[i] + Q->tab[i];
        M->tab[i] = Ms + Mp - Q->tab[i];
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
