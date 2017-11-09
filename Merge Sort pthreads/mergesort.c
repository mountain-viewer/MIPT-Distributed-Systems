#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>
#include <time.h>
#include <pthread.h>

// Represent a task for an array to be sorted
typedef struct {
    int n;
    int m;
    int numberOfThreads;
    
    int *array;
    int *buffer;
    pthread_t *threads;
    
} SortTask;

// Represent a task for two parts to be merged
typedef struct {
    int leftSize;
    int rightSize;
    int m;
    int numberOfThreads;
    
    int *left;
    int *right;
    
    int *target;
    int *buffer;
    
    pthread_t *threads;
} MergeTask;

// Simplify the work with command line arguments
void parseInputData(int* arraySize, int* baseChunkSize, int* numberOfThreads, const char * argv[]) {
    *arraySize = atoi(argv[1]);
    *baseChunkSize = atoi(argv[2]);
    *numberOfThreads = atoi(argv[3]);
}

// Comparator function for sorting
int less(const void *a, const void *b) {
    return (*(int*)a - *(int*)b);
}

// Implement a sequential algorithm for merging
void sequentialMerge(MergeTask *task) {
    int i = 0;
    int j = 0;
    
    while (i < task->leftSize && j < task->rightSize) {
        if (task->left[i] < task->right[j]) {
            task->buffer[i + j] = task->left[i];
            ++i;
        } else {
            task->buffer[i + j] = task->right[j];
            ++j;
        }
    }
    
    while (i < task->leftSize) {
        task->buffer[i + j] = task->left[i];
        ++i;
    }
    while (j < task->rightSize) {
        task->buffer[i + j] = task->right[j];
        ++j;
    }
    
    if (task->target != task->buffer) {
        int n = task->leftSize + task->rightSize;
        memcpy(task->target, task->buffer, n * sizeof(int));
    }
}

// Implement a sequential algorithm for merge sort (P = 1)
void sequentialMergeSort(SortTask *task) {
    if (task->m >= task->n) {
        // Base recursion case -- using quick sort for small chunks
        qsort(task->array, task->n, sizeof(int), less);
        return;
    }
    
    int middle = task->n / 2;
    
    // Left subtask
    SortTask leftTask;
    leftTask.n = middle;
    leftTask.m = task->m;
    leftTask.array = task->array;
    leftTask.buffer = task->buffer;
    leftTask.numberOfThreads = 0;
    leftTask.threads = NULL;
    
    // Right subtask
    SortTask rightTask;
    rightTask.n = task->n - middle;
    rightTask.m = task->m;
    rightTask.array = task->array + middle;
    rightTask.buffer = task->buffer + middle;
    rightTask.numberOfThreads = 0;
    rightTask.threads = NULL;
    
    // Implement divide and conquer method -- sorting subarrays recursively
    sequentialMergeSort(&leftTask);
    sequentialMergeSort(&rightTask);
    
    MergeTask mergeTask;
    mergeTask.leftSize = middle;
    mergeTask.rightSize = task->n - middle;
    mergeTask.m = task->m;
    mergeTask.left = task->array;
    mergeTask.right = task->array + middle;
    mergeTask.target = task->array;
    mergeTask.buffer = task->buffer;
    mergeTask.numberOfThreads = 0;
    mergeTask.threads = NULL;
    
    sequentialMerge(&mergeTask);
}

// Implement a binary search
int binarySearch(int *array, int left, int right, int value) {
    int middle = (left + right) / 2;
    
    if (left == right) {
        return right;
    } else if (array[middle] == value) {
        return middle;
    } else if (array[middle] > value) {
        return binarySearch(array, left, middle, value);
    } else {
        return binarySearch(array, middle + 1, right, value);
    }
}

// Implement a parallel merge algorithm
// Signature: void *(*start_routine) (void *) -- because of pthread_create arguments
void* parallelMerge(void *arg) {
    MergeTask *task = (MergeTask*)arg;
    
    if (task->numberOfThreads <= 1 ||  task->leftSize <= task->m || task->rightSize <= task->m) {
        sequentialMerge(task);
        return NULL;
    }
    
    int leftMiddle = task->leftSize / 2;
    int value = task->left[leftMiddle];
    int rightIndex = binarySearch(task->right, 0, task->rightSize, value);
    
    // Left subtask
    MergeTask leftTask;
    leftTask.leftSize = leftMiddle;
    leftTask.rightSize = rightIndex;
    leftTask.m = task->m;
    leftTask.left = task->left;
    leftTask.right = task->right;
    leftTask.target = task->target;
    leftTask.buffer = task->buffer;
    leftTask.threads = task->threads + 2;
    leftTask.numberOfThreads = (task->numberOfThreads - 2) / 2;
    
    
    // Right subtask
    MergeTask rightTask;
    rightTask.leftSize = task->leftSize - leftMiddle;
    rightTask.rightSize = task->rightSize - rightIndex;
    rightTask.m = task->m;
    rightTask.left = task->left + leftMiddle;
    rightTask.right = task->right + rightIndex;
    rightTask.target = task->target + leftMiddle + rightIndex;
    rightTask.buffer = task->buffer + leftMiddle + rightIndex;
    rightTask.threads = task->threads + 2 + (task->numberOfThreads - 2) / 2;
    rightTask.numberOfThreads = task->numberOfThreads - 2 - (task->numberOfThreads - 2) / 2;

    pthread_t leftThread = *(task->threads);
    pthread_t rightThread = *(task->threads + 1);
    
    // Thread calls
    pthread_create(&leftThread, NULL, parallelMerge, &leftTask);
    pthread_create(&rightThread, NULL, parallelMerge, &rightTask);
    
    // Sync
    pthread_join(leftThread, NULL);
    pthread_join(rightThread, NULL);
    
    return NULL;
}

// Implement a parallel merge sort algorithm
// Signature: void *(*start_routine) (void *) -- because of pthread_create arguments
void* parallelMergeSort(void *arg) {
    SortTask *task = (SortTask*)arg;
    
    if (task->numberOfThreads <= 1) {
        // Single thread
        sequentialMergeSort(task);
        return NULL;
    }
    
    if (task->m >= task->n) {
        // Base recursion case -- using quick sort for small chunks
        qsort(task->array, task->n, sizeof(int), less);
        return NULL;
    }
    
    int middle = task->n / 2;
    
    // Left subtask
    SortTask leftTask;
    leftTask.n = middle;
    leftTask.m = task->m;
    leftTask.array = task->array;
    leftTask.buffer = task->buffer;
    leftTask.threads = task->threads + 2;
    leftTask.numberOfThreads = (task->numberOfThreads - 2) / 2;
    
    // Right subtask
    SortTask rightTask;
    rightTask.n = task->n - middle;
    rightTask.m = task->m;
    rightTask.array = task->array + middle;
    rightTask.buffer = task->buffer + middle;
    rightTask.threads = task->threads + 2 + (task->numberOfThreads - 2) / 2;
    rightTask.numberOfThreads = task->numberOfThreads - 2 - (task->numberOfThreads - 2) / 2;
    
    pthread_t leftThread = *(task->threads);
    pthread_t rightThread = *(task->threads + 1);
    
    // Thread call
    pthread_create(&leftThread, NULL, parallelMergeSort, &leftTask);
    pthread_create(&rightThread, NULL, parallelMergeSort, &rightTask);
    
    // Sync
    pthread_join(leftThread, NULL);
    pthread_join(rightThread, NULL);
    
    MergeTask mergeTask;
    mergeTask.leftSize = middle;
    mergeTask.rightSize = task->n - middle;
    mergeTask.m = task->m;
    mergeTask.left = task->array;
    mergeTask.right = task->array + middle;
    mergeTask.target = task->buffer;
    mergeTask.buffer = task->buffer;
    mergeTask.threads = task->threads;
    mergeTask.numberOfThreads = task->numberOfThreads;
    
    parallelMerge(&mergeTask);
    memcpy(task->array, mergeTask.target, task->n * sizeof(int));
    
    return NULL;
}

int main(int argc, const char * argv[]) {
    if (argc != 4) {
        printf("Usage: %s n m P\n", argv[0]);
        return EXIT_FAILURE;
    }
    
    int n, m, P;
    parseInputData(&n, &m, &P, argv);
    
    int *array = (int*)malloc(n * sizeof(int));
    int *buffer = (int*)malloc(n * sizeof(int));
    pthread_t threads[P];
    
    FILE *data = fopen("data.txt", "w");
    srand(time(NULL));
    
    // Filling the array with random numbers
    for (int i = 0; i < n; ++i) {
        array[i] = rand();
    }
    
    // Writing an array to file
    for (int i = 0; i < n; ++i) {
        fprintf(data, "%d ", array[i]);
    }
    fprintf(data, "\n");
    
    SortTask task;
    task.n = n;
    task.m = m;
    task.array = array;
    task.buffer = buffer;
    task.threads = threads;
    task.numberOfThreads = P;

    double startTime = omp_get_wtime();
    
    if (P == 1) {
        sequentialMergeSort(&task);
    } else {
        parallelMergeSort(&task);
    }
    double endTime = omp_get_wtime();
    double totalTime = endTime - startTime;
    
    FILE* stats = fopen("stats.txt", "w");
    fprintf(stats, "%.5fs %d %d %d\n", totalTime, n, m, P);
    fclose(stats);
    
    for (int i = 0; i < n; i++) {
        fprintf(data, "%d ", task.array[i]);
    }
    
    fprintf(data, "\n");
    fclose(data);
    
    free(array);
    free(buffer);
    
    return EXIT_SUCCESS;
}

