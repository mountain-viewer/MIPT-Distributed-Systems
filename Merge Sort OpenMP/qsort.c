#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

int less(const void *a, const void *b) {
    return (*(int*)a - *(int*)b);
}

int main(int argc, const char * argv[]) {
    if (argc != 2) {
        printf("Usage: %s n\n", argv[0]);
        return EXIT_FAILURE;
    }
    
    int n = atoi(argv[1]);
    
    // Creating an array with random values of size n
    int *array = (int*) malloc(n * sizeof(int));
    srand(time(NULL));
    for (int i = 0; i < n; i++) {
        array[i] = rand();
    }
    
    double startTime = omp_get_wtime();
    qsort(array, n, sizeof(int), less);
    double endTime = omp_get_wtime();
    
    printf("%.5fs\n", endTime - startTime);
    
    free(array);
    
    return EXIT_SUCCESS;
}
