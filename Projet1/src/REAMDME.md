# Recherche du sous tableau maximal

* `void printTablo(struct tablo *tmp)` : non parallele
* `struct tablo *allocateTablo(int size)` : non parallele
* `void freeTablo(struct tablo *tmp)` : non parallele
* `struct tablo *parseFile(char *filePath)` : non parallele
* `void montee(struct tablo *source, struct tablo *dest)` : parallele
* `void descente(struct tablo *a, struct tablo *b)` : parallele
* `void descenteSuffixe(struct tablo *a, struct tablo *b)` : parallele
* `void monteeMax(struct tablo *source, struct tablo *dest)` : parallele
* `void descenteMax(struct tablo *a, struct tablo *b)` : parallele
* `void descenteSuffixeMax(struct tablo *a, struct tablo *b)` : parallele
* `void final(struct tablo *a, struct tablo *b)` : parallele
* `void finalMax(struct tablo *a, struct tablo *b)` : parallele
* `void prefixSum(struct tablo *source, struct tablo *dest)` : parallele
* `void suffixSum(struct tablo *source, struct tablo *dest)` : parallele
* `void suffixMax(struct tablo *source, struct tablo *dest)` : parallele
* `void prefixMax(struct tablo *source, struct tablo *dest)` : parallele
* `void displayResult(struct tablo *M, struct tablo *source)` : non parallele
* `int main(int argc, char **argv)` : parallele