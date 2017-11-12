#include <stdio.h>
#include <mpi.h>
#include <pthread.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <string.h>
#include <time.h>

int factor = 2;
int timeInterval = 200;

// Simplify the work with command line arguments
void parseInputData(int *l, int *a, int *b, int *n, int *N,
                    float *pl, float *pr, float *pu, float *pd, char *argv[]) {
    *l = atoi(argv[1]);
    *a = atoi(argv[2]);
    *b = atoi(argv[3]);
    *n = atoi(argv[4]);
    *N = atoi(argv[5]);
    
    *pl = atof(argv[6]);
    *pr = atof(argv[7]);
    *pu = atof(argv[8]);
    *pd = atof(argv[9]);
}

// Work with boolean values in order to make code more readable
typedef int bool;

bool true() {
    return 1;
}

bool false() {
    return 0;
}

// Direction type
typedef enum {
    LEFT = 1,
    RIGHT = 2,
    UP = 3,
    DOWN = 4
} Direction;

// Settings for the process simulation
typedef struct {
    int l;
    int a;
    int b;
    int n;
    int N;
    
    float pl;
    float pr;
    float pu;
    float pd;
    
    int rank;
    int size;
} TaskSetting;

// Point type
typedef struct {
    int x;
    int y;
} Point;

// Particle type
typedef struct {
    Point point;
    int n;
    int seed;
} Particle;

// Calculate the point by given rank and a pair (a, b)
Point getPoint(int rank, int a, int b) {
    Point point;
    point.x = rank % a;
    point.y = rank / a;
    return point;
}

// Calculate rank by given point and a pair (a, b)
int getRank(Point point, int a, int b) {
    return point.y * a + point.x;
}

// Calculate rank of nearest neighbours
int getNeighbourRank(Point point, int a, int b, Direction direction) {
    
    // Left direction
    if (direction == LEFT) {
        point.x -= 1;
        if (point.x < 0) {
            point.x = a - 1;
        }
    }
    
    // Right direction
    if (direction == RIGHT) {
        point.x += 1;
        if (point.x >= a) {
            point.x = 0;
        }
    }
    
    // Up direction
    if (direction == UP) {
        point.y -= 1;
        if (point.y < 0) {
            point.y = b - 1;
        }
    }
    
    // Down direction
    if (direction == DOWN) {
        point.y += 1;
        if (point.y >= b) {
            point.y = 0;
        }
    }
    return getRank(point, a, b);
}

// Append an element to the end of array
// Amortized complexity: O(1)
void append(Particle **array, int *n, int *capacity, Particle *element) {
    if (*n == *capacity) {
        *array = realloc(*array, *capacity * factor * sizeof(Particle));
        *capacity *= factor;
        (*array)[*n] = *element;
        (*n)++;
        return;
    }
    
    (*array)[*n] = *element;
    (*n)++;
}

// Pop the last element from the array
void pop(Particle **array, int *n, int index) {
    (*array)[index] = (*array)[(*n) - 1];
    (*n)--;
}

// Selecting the most probable direction
Direction getDirection(float leftProbability, float rightProbability, float upProbability, float downProbability) {
    if (leftProbability >= rightProbability && leftProbability >= upProbability && leftProbability >= downProbability) {
        return LEFT;
    } else if (rightProbability >= leftProbability && rightProbability >= upProbability && rightProbability >= downProbability) {
        return RIGHT;
    } else if (upProbability >= leftProbability && upProbability >= rightProbability && upProbability >= downProbability) {
        return UP;
    } else {
        return DOWN;
    }
}

void* processArea(void *arguments) {
    TaskSetting *setting = (TaskSetting*)arguments;
    
    int l = setting->l;
    int a = setting->a;
    int b = setting->b;
    int n = setting->n;
    int N = setting->N;
    
    float pl = setting->pl;
    float pr = setting->pr;
    float pu = setting->pu;
    float pd = setting->pd;
    
    int rank = setting->rank;
    int size = setting->size;
    
    
    Point point = getPoint(rank, a, b);
    
    int leftRank = getNeighbourRank(point, a, b, LEFT);
    int rightRank = getNeighbourRank(point, a, b, RIGHT);
    int upRank = getNeighbourRank(point, a, b, UP);
    int downRank = getNeighbourRank(point, a, b, DOWN);

    int seed;
    int *seeds = malloc(size * sizeof(int));
    
    // The process with the least rank seeds the random generator for other processes
    if (rank == 0) {
        srand(time(NULL));
        
        for (int i = 0; i < size; ++i) {
            seeds[i] = rand();
        }
    }

    MPI_Scatter(seeds, 1, MPI_INT, &seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    int particlesSize = N;
    int particlesCapacity = N;
    
    Particle *particles = (Particle*)malloc(particlesCapacity * sizeof(Particle));
    srand(seed);
    
    for (int i = 0; i < particlesSize; ++i) {
        particles[i].point.x = rand() % l;
        particles[i].point.y = rand() % l;
        particles[i].n = n;
        particles[i].seed = rand();
    }
    
    // Left case
    int sendLeftSize = 0;
    int sendLeftCapacity = N;
    Particle *sendLeft = (Particle*)malloc(sendLeftCapacity * sizeof(Particle));
    
    // Right case
    int sendRightSize = 0;
    int sendRightCapacity = N;
    Particle *sendRight = (Particle*)malloc(sendRightCapacity * sizeof(Particle));
    
    // Up case
    int sendUpSize = 0;
    int sendUpCapacity = N;
    Particle *sendUp = (Particle*)malloc(sendUpCapacity * sizeof(Particle));
    
    // Down case
    int sendDownSize = 0;
    int sendDownCapacity = N;
    Particle *sendDown = (Particle*)malloc(sendDownCapacity * sizeof(Particle));
    
    int receiveLeftCapacity = 0;
    int receiveRightCapacity = 0;
    int receiveUpCapacity = 0;
    int receiveDownCapacity = 0;
    
    int completedSize = 0;
    int completedCapacity = N;
    
    Particle *completed = (Particle*)malloc(completedCapacity * sizeof(Particle));
    
    struct timeval t0, t1;
    assert(gettimeofday(&t0, NULL) == 0);
    
    while (true()) {
        int particleID = 0;
        while (particleID < particlesSize) {
            
            Particle *particle = particles + particleID;
            
            bool needIncrement = true();
            for (int time = 0; time < timeInterval; ++time) {
                if (particle->n == 0) {
                    append(&completed, &completedSize, &completedCapacity, particle);
                    pop(&particles, &particlesSize, particleID);
                    needIncrement = false();
                    break;
                }
                
                float leftProbability = rand_r((unsigned int*) &particle->seed) * pl;
                float rightProbability = rand_r((unsigned int*) &particle->seed) * pr;
                float upProbability = rand_r((unsigned int*) &particle->seed) * pu;
                float downProbability = rand_r((unsigned int*) &particle->seed) * pd;
                
                Direction direction = getDirection(leftProbability, rightProbability, upProbability, downProbability);
                
                if (direction == LEFT) {
                    particle->point.x -= 1;
                } else if (direction == RIGHT) {
                    particle->point.x += 1;
                } else if (direction == UP) {
                    particle->point.y -= 1;
                } else {
                    particle->point.y += 1;
                }
                
                particle->n -= 1;
                
                if (particle->point.x < 0) {
                    particle->point.x = l - 1;
                    append(&sendLeft, &sendLeftSize, &sendLeftCapacity, particle);
                    pop(&particles, &particlesSize, particleID);
                    needIncrement = false();
                    break;
                }
                
                if (particle->point.x >= l) {
                    particle->point.x = 0;
                    append(&sendRight, &sendRightSize, &sendRightCapacity, particle);
                    pop(&particles, &particlesSize, particleID);
                    needIncrement = false();
                    break;
                }
                if (particle->point.y < 0) {
                    particle->point.y = l - 1;
                    append(&sendUp, &sendUpSize, &sendUpCapacity, particle);
                    pop(&particles, &particlesSize, particleID);
                    needIncrement = false();
                    break;
                }
                if (particle->point.y >= l) {
                    particle->point.y = 0;
                    append(&sendDown, &sendDownSize, &sendDownCapacity, particle);
                    pop(&particles, &particlesSize, particleID);
                    needIncrement = false();
                    break;
                }
            }
            
            if (needIncrement) {
                particleID += 1;
            }
        }
        
        // Sending and receiving basic info
        MPI_Request *metadata = (MPI_Request*)malloc(8 * sizeof(MPI_Request));
        MPI_Isend(&sendLeftSize, 1, MPI_INT, leftRank, 0, MPI_COMM_WORLD, metadata + 0);
        MPI_Isend(&sendRightSize, 1, MPI_INT, rightRank, 1, MPI_COMM_WORLD, metadata + 1);
        MPI_Isend(&sendUpSize, 1, MPI_INT, upRank, 2, MPI_COMM_WORLD, metadata + 2);
        MPI_Isend(&sendDownSize, 1, MPI_INT, downRank, 3, MPI_COMM_WORLD, metadata + 3);

        MPI_Irecv(&receiveLeftCapacity, 1, MPI_INT, leftRank, 1, MPI_COMM_WORLD, metadata + 4);
        MPI_Irecv(&receiveRightCapacity, 1, MPI_INT, rightRank, 0, MPI_COMM_WORLD, metadata + 5);
        MPI_Irecv(&receiveUpCapacity, 1, MPI_INT, upRank, 3, MPI_COMM_WORLD ,metadata + 6);
        MPI_Irecv(&receiveDownCapacity, 1, MPI_INT, downRank, 2, MPI_COMM_WORLD, metadata + 7);

        // Wait for all processes
        MPI_Waitall(8, metadata, MPI_STATUS_IGNORE);
        
        Particle *receiveLeft = (Particle*)malloc(receiveLeftCapacity * sizeof(Particle));
        Particle *receiveRight = (Particle*)malloc(receiveRightCapacity * sizeof(Particle));
        Particle *receiveUp = (Particle*)malloc(receiveUpCapacity * sizeof(Particle));
        Particle *receiveDown = (Particle*)malloc(receiveDownCapacity * sizeof(Particle));
        
        MPI_Request *data = (MPI_Request*)malloc(8 * sizeof(MPI_Request));
        MPI_Issend(sendLeft, sizeof(Particle) * sendLeftSize, MPI_BYTE, leftRank, 0, MPI_COMM_WORLD, data + 0);
        MPI_Issend(sendRight, sizeof(Particle) * sendRightSize, MPI_BYTE, rightRank, 1, MPI_COMM_WORLD, data + 1);
        MPI_Issend(sendUp, sizeof(Particle) * sendUpSize, MPI_BYTE, upRank, 2, MPI_COMM_WORLD, data + 2);
        MPI_Issend(sendDown, sizeof(Particle) * sendDownSize, MPI_BYTE, downRank, 3, MPI_COMM_WORLD, data + 3);

        MPI_Irecv(receiveLeft, sizeof(Particle) * receiveLeftCapacity, MPI_BYTE, leftRank, 1, MPI_COMM_WORLD, data + 4);
        MPI_Irecv(receiveRight, sizeof(Particle) * receiveRightCapacity, MPI_BYTE, rightRank, 0, MPI_COMM_WORLD, data + 5);
        MPI_Irecv(receiveUp, sizeof(Particle) * receiveUpCapacity, MPI_BYTE, upRank, 3, MPI_COMM_WORLD, data + 6);
        MPI_Irecv(receiveDown, sizeof(Particle) * receiveDownCapacity, MPI_BYTE, downRank, 2, MPI_COMM_WORLD, data + 7);

        // Waiting till all the processes reach this line
        MPI_Waitall(8, data, MPI_STATUS_IGNORE);

        free(metadata);
        free(data);
        
        // Add received particles to current particles
        for (int offset = 0; offset < receiveLeftCapacity; ++offset) {
            append(&particles, &particlesSize, &particlesCapacity, receiveLeft + offset);
        }
        for (int offset = 0; offset < receiveRightCapacity; ++offset) {
            append(&particles, &particlesSize, &particlesCapacity, receiveRight + offset);
        }
        for (int offset = 0; offset < receiveUpCapacity; ++offset) {
            append(&particles, &particlesSize, &particlesCapacity, receiveUp + offset);
        }
        for (int offset = 0; offset < receiveDownCapacity; ++offset) {
            append(&particles, &particlesSize, &particlesCapacity, receiveDown + offset);
        }
        
        free(receiveLeft);
        free(receiveRight);
        free(receiveUp);
        free(receiveDown);
        
        // All previous particles have been sent
        sendLeftSize = 0;
        sendRightSize = 0;
        sendUpSize = 0;
        sendDownSize = 0;

        // Sending info to the master process
        int *recvBuffer = (int*)malloc(sizeof(int));
        MPI_Reduce(&completedSize, recvBuffer, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        
        // Synchronizing
        MPI_Barrier(MPI_COMM_WORLD);
        
        bool finished = false();
        if (rank == 0) {
            if (*recvBuffer == size * N) {
                finished = true();
            } else {
                finished = false();
            }
        }
        MPI_Bcast(&finished, 1, MPI_INT, 0, MPI_COMM_WORLD);

        free(recvBuffer);
        
        if (finished) {
            int *total = (int*)malloc(size * sizeof(int));
            MPI_Gather(&completedSize, 1, MPI_INT, total, 1, MPI_INT, 0, MPI_COMM_WORLD);
            
            if (rank == 0) {
                assert(gettimeofday(&t1, NULL) == 0);
                double delta = ((t1.tv_sec - t0.tv_sec) * 1000000 + t1.tv_usec - t0.tv_usec) / 1000000.0;
                
                FILE *stats = fopen("stats.txt", "w");
                fprintf(stats, "%d %d %d %d %d %f %f %f %f %fs\n", l, a, b, n, N, pl, pr, pu, pd, delta);
                for (int rank = 0; rank < size; ++rank) {
                    fprintf(stats, "%d: %d\n", rank, total[rank]);
                }
                fclose(stats);
            }
            free(total);
            break;
        }
        
        // Synchronizing here
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    free(seeds);

    free(sendLeft);
    free(sendRight);
    free(sendUp);
    free(sendDown);

    free(particles);
    free(completed);

    return NULL;
}

int main(int argc, char * argv[]) {
    if (argc != 10) {
        printf("Usage: %s l a b n N pl pr pu pd\n", argv[0]);
        return EXIT_FAILURE;
    }
    
    int l, a, b, n, N;
    float pl, pr, pu, pd;
    parseInputData(&l, &a, &b, &n, &N, &pl, &pr, &pu, &pd, argv);
    
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    TaskSetting setting;
    setting.l = l;
    setting.a = a;
    setting.b = b;
    setting.n = n;
    setting.N = N;
    
    setting.pl = pl;
    setting.pr = pr;
    setting.pu = pu;
    setting.pd = pd;
    
    setting.rank = rank;
    setting.size = size;
    
    
    pthread_t threadID;
    pthread_create(&threadID, NULL, processArea, &setting);
    pthread_join(threadID, NULL);
    
    MPI_Finalize();
    return EXIT_SUCCESS;
}
