gcc -Wall -o base-omp base.c -lm -fopenmp
gcc -Wall -o base base.c

./base 3000 > example-3000.data
./base-omp 3000 > example-3000-omp.data

./base 10000 > example-10000.data
./base-omp 10000 > example-10000-omp.data