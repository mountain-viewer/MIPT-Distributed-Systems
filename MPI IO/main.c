#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

// Simplify the work with command line arguments
void parseInputData(int *l, int *a, int *b, int *N, char *argv[]) {
    *l = atoi(argv[1]);
    *a = atoi(argv[2]);
    *b = atoi(argv[3]);
    *N = atoi(argv[4]);
}

typedef struct {
    int x;
    int y;
    int r;
} Particle;

// Settings for the process for simulation
typedef struct {
    int rank;
    int size;
    
    int l;
    int a;
    int b;
    int N;
} TaskSetting;

// Simulate two-dimensional random walk
void run(TaskSetting setting) {
    double startTime = MPI_Wtime();
    
    // The coordinates of square that owes the current process
    int x = setting.rank % setting.a;
    int y = setting.rank / setting.a;
    
    Particle *particles = (Particle*)malloc(setting.N * sizeof(Particle));
    int seed;
    int *seeds = (int*)malloc(setting.size * sizeof(int));
    
    // The process with the least rank seeds the random generator for other processes
    if (setting.rank == 0) {
        srand(time(NULL));
        
        for (int i = 0; i < setting.size; ++i) {
            seeds[i] = rand();
        }
    }
    
    MPI_Scatter(seeds, 1, MPI_UNSIGNED, &seed, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    free(seeds);
    
    for (int i = 0; i < setting.N; ++i) {
        particles[i].x = rand_r(&seed) % setting.l;
        particles[i].y = rand_r(&seed) % setting.l;
        particles[i].r = rand_r(&seed) % (setting.a * setting.b);
    }
    
    int maskSize = setting.l * setting.l * setting.size;
    int *mask = (int*)malloc(maskSize * sizeof(int));
    
    for (int i = 0; i < maskSize; ++i) {
        mask[i] = 0;
    }
    
    MPI_File binaryFile;
    MPI_File_delete("data.bin", MPI_INFO_NULL);
    MPI_File_open(MPI_COMM_WORLD, "data.bin", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &binaryFile);
    
    for (int i = 0; i < setting.N; ++i) {
        int currentX = particles[i].x;
        int currentY = particles[i].y;
        int currentR = particles[i].r;
        mask[currentY * setting.l * setting.size + currentX * setting.size + currentR] += 1;
    }
    
    MPI_Aint intex;
    MPI_Aint e;
    MPI_Type_get_extent(MPI_INT, &e, &intex);
    
    MPI_Datatype view;
    MPI_Type_vector(setting.l, setting.l * setting.size, setting.l * setting.a * setting.size, MPI_INT, &view);
    MPI_Type_commit(&view);
    
    int offset = (x * setting.l + y * setting.a * setting.l * setting.l) * setting.size;
    MPI_File_set_view(binaryFile, offset * sizeof(int), MPI_INT, view, "native", MPI_INFO_NULL);
    
    MPI_File_write(binaryFile, mask, maskSize, MPI_INT, MPI_STATUS_IGNORE);
    MPI_Type_free(&view);
    
    MPI_File_close(&binaryFile);
    
    free(mask);
    
    double finishTime = MPI_Wtime();
    
    // Output the stats
    if (setting.rank == 0) {
        double overallTime = finishTime - startTime;
        FILE *stats = fopen("stats.txt", "w");
        
        fprintf(stats, "%d %d %d %d %fs\n",
                setting.l,
                setting.a,
                setting.b,
                setting.N,
                overallTime);
        
        fclose(stats);
    }
    
    free(particles);
}


int main(int argc, char * argv[]) {
    if (argc != 5) {
        printf("Usage: %s l a b N\n", argv[0]);
        return EXIT_FAILURE;
    }
    
    int l, a, b, N;
    parseInputData(&l, &a, &b, &N, argv);

    MPI_Init(&argc, &argv);
    double startTime = MPI_Wtime();

    int rank, size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size != a * b) {
        printf("The condition from the statement is not satisfied. (size != a * b)\n");
        return EXIT_FAILURE;
    }
    
    TaskSetting setting;
    setting.rank = rank;
    setting.size = size;
    
    setting.l = l;
    setting.a = a;
    setting.b = b;
    setting.N = N;

    run(setting);

    MPI_Finalize();
    return EXIT_SUCCESS;
}
