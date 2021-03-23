# Un peu de C

## Question 1

```c
int* tableauAlea(int taille) {
  return (int*) malloc(sizeof(int) * taille);
}



int main()
{
  int size = 10;

  int* tab = tableauAlea(size);
}
```

## Question 2

* Sous la forme `T[]`, avec `T` un objet quelconque.
* Sous la forme d'un pointeur.

## Question 3

* Taille maximum d'un tableau alloué sur la pile :
* Taille maximum d'un tableau alloué sur le tas :

## Question 4

```c
void rempliAleatoire(int* tab, int taille, int max)
{
  for (int i = 0; i < taille; i++)
    tab[i] = (int) rand() % max;
}
```