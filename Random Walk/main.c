//
//  main.c
//  random_walk
//
//  Created by whoami on 10/16/17.
//  Copyright © 2017 Mountain Viewer. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>

void parseInputData(int *a, int *b, int *x, int *N,
                    float *p,
                    int *P,
                    const char *argv[]) {
    *a = atoi(argv[1]);
    *b = atoi(argv[2]);
    *x = atoi(argv[3]);
    *N = atoi(argv[4]);
    *p = atof(argv[5]);
    *P = atoi(argv[6]);
}

typedef struct {
    int time;
    int isRightEndPoint;
} RandomWalk;

// Returns the 1 or -1 with respect to probability
int randomWalkStep(float p, unsigned int *seed) {
    int threshold = p * RAND_MAX;
    int result = rand_r(seed);
    
    if (result < threshold) {
        return 1;
    }
    
    return -1;
}

RandomWalk runRandomWalk(int a, int b, int x, float p, unsigned int *seed) {
    RandomWalk randomWalk;
    int currentPoint = x;
    randomWalk.time = 0;

    while (currentPoint != a && currentPoint != b) {
        currentPoint += randomWalkStep(p, seed);
        randomWalk.time += 1;
    }

    if (currentPoint == b) {
        randomWalk.isRightEndPoint = 1;
    } else {
        randomWalk.isRightEndPoint = 0;
    }

    return randomWalk;
}

float calculateDeltaTime(struct timeval begin, struct timeval end) {
    return ((end.tv_sec  - begin.tv_sec) * 1000000u +
            end.tv_usec - begin.tv_usec) / 1.e6;
}

int main(int argc, const char * argv[]) {
    int a, b, x, N;
    float p;
    int P;
    
    if (argc != 7) {
        printf("Usage: %s a b x N p P\n", argv[0]);
        return EXIT_FAILURE;
    }

    parseInputData(&a, &b, &x, &N, &p, &P, argv);
    printf("%d %d %d %d %.2f %d\n", a, b, x, N, p, P);
    // Seeds the pseudo random generator for rand()
    srand(time(NULL));
   
    int rightEndPointCount = 0;
    int overallTime = 0;
    struct timeval start, finish;
    
    gettimeofday(&start, NULL);
    
    omp_set_num_threads(P);
    
    #pragma omp parallel
    {
        unsigned int seed = rand();
        
        #pragma omp for reduction(+: rightEndPointCount, overallTime)
        for (int i = 0; i < N; ++i) {
            RandomWalk randomWalk = runRandomWalk(a, b, x, p, &seed);
            rightEndPointCount += randomWalk.isRightEndPoint;
            overallTime += randomWalk.time;
        }
    }

    gettimeofday(&finish, NULL);
    
    float delta = calculateDeltaTime(start, finish);
    
    FILE *file = fopen("stats.txt", "w");
    fprintf(file, "%.2f %.1f %.4fs %d %d %d %d %.2f %d\n",
           rightEndPointCount * 1.0 / N,
           overallTime * 1.0 / N,
           delta,
           a, b, x, N, p, P);
    fclose(file);

    return EXIT_SUCCESS;
}
