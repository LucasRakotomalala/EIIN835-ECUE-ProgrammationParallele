# Recherche du sous tableau maximal

* `void printTablo(struct tablo *tmp)` : non parallele
* `struct tablo *allocateTablo(int size)` : non parallele
* `void freeTablo(struct tablo *tmp)` : non parallele
* `struct tablo *parseFileAndFillTablo(char *filePath)` : non parallele
* `void up(struct tablo *source, struct tablo *dest)` : parallele
* `void down(struct tablo *a, struct tablo *b)` : parallele
* `void downSuffix(struct tablo *a, struct tablo *b)` : parallele
* `void upMax(struct tablo *source, struct tablo *dest)` : parallele
* `void downMax(struct tablo *a, struct tablo *b)` : parallele
* `void downMaxSuffix(struct tablo *a, struct tablo *b)` : parallele
* `void final(struct tablo *a, struct tablo *b)` : parallele
* `void finalMax(struct tablo *a, struct tablo *b)` : parallele
* `void prefixSum(struct tablo *source, struct tablo *dest)` : parallele
* `void suffixSum(struct tablo *source, struct tablo *dest)` : parallele
* `void suffixMax(struct tablo *source, struct tablo *dest)` : parallele
* `void prefixMax(struct tablo *source, struct tablo *dest)` : parallele
* `void displayResult(struct tablo *M, struct tablo *source)` : non parallele
* `int main(int argc, char **argv)` : parallele

## Compilation du fichier source

```c
gcc -Wall -std=c99 -o rakotomalala rakotomalala.c -lm -fopenmp
```