#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

int* tableauAlea(int taille) {
  return (int*) malloc(sizeof(int) * taille);
}

void rempliAleatoire(int* tab, int taille, int max)
{
  #pragma omp parallel for
  for (int i = 0; i < taille; i++)
  {
    tab[i] = (int) rand() % max;
  }
}

// ------------------------------------

int isAPrime(int i)
{
   int divider = 2;
   while (divider < i) {
     if (i % divider == 0) {
       return 0;
     }
     divider++;
   }
   return 1;
}

void simpleLoop(int n)
{
   for (int i = 0; i < n; i++) {
     if (i == n/2) {
       printf("I'm halfway there !\n");
     }
   }
    printf("I'm done !\n");
}

void cpuIntensiveLoop(int n)
{
  struct timeval start, end;
  gettimeofday(&start, NULL);

  #pragma omp parallel for
  for (int i = 0; i < n; i++) {
      isAPrime(i);
  }

  gettimeofday(&end, NULL);
  printf("%ld\n", ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)));
}

int* testPrime(int* tab, int size) {
  int* tabBool = (int*) malloc(sizeof(int) * size);

  #pragma omp parallel for
  for (int i = 0; i < size; i++)
  {
    tabBool[i] = isAPrime(tab[i]);
  }
  return tabBool;
}

int main(int argc, char* argv[])
{
  int taille = (argc == 2) ? atoi(argv[1]) : 1000;
  int size = 0;

  for (int i = 0; i < taille; i = i + 10) {
    struct timeval start, end;
    gettimeofday(&start, NULL);

    size = i;

    int* tab = tableauAlea(size);
    rempliAleatoire(tab, size, 1000000);
    /*
    for (int i = 0; i < size; i++)
      printf("%i_", tab[i]);
    printf("\n");
    */
    //cpuIntensiveLoop(size);
    int* tabBool = testPrime(tab, size);
    /*
    for (int i = 0; i < size; i++)
      printf("%i ", tabBool[i]);
    printf("\n");
    */
    free(tab);
    free(tabBool);

    gettimeofday(&end, NULL);
    printf("%i %ld\n", size, ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)));
}

  return 0;
}