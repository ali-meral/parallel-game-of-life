#include "simulation_control.h"
#include "matrix_operations.h"
#include "utilities.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

void run_simulation(int argc, char *argv[])
{
    int rank, size; // MPI rank and size
    int seed = 10;
    int n = 10; // Grid size
    int verbose = 0;
    int debug = 0;
    int density = 10;   // Density of alive cells
    int iterations = 3; // Number of generations

    // Setting up MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Command line options
    static struct option long_options[] = {
        {"number", required_argument, NULL, 'n'},
        {"seed", required_argument, NULL, 's'},
        {"density", required_argument, NULL, 'd'},
        {"verbose", no_argument, NULL, 'v'},
        {"debug", no_argument, NULL, 'd'},
        {"iterations", required_argument, NULL, 'i'},
        {NULL, 0, NULL, 0}};

    int opt;
    while ((opt = getopt_long(argc, argv, "n:s:d:vbi:", long_options, NULL)) != -1)
    {
        switch (opt)
        {
        case 'n':
            n = atoi(optarg);
            break;
        case 's':
            seed = atoi(optarg);
            break;
        case 'd':
            density = atoi(optarg);
            break;
        case 'i':
            iterations = atoi(optarg);
            break;
        case 'v':
            verbose = 1;
            break;
        case 'b':
            debug = 1;
            break;
        default:
            if (rank == 0)
            {
                fprintf(stderr, "Usage: %s -n <size> -s <seed> -d <density> -i <iterations> -v (verbose) -b (debug)\n", argv[0]);
            }
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
    }

    // Simulation parameters
    int n_loc_r = n; // Local rows
    int n_loc_c = n; // Local columns

    if (verbose && rank == 0)
    {
        printf("Starting simulation with %d x %d grid, density %d%%, %d iterations, seed %d\n", n, n, density, iterations, seed);
    }

    // Allocate matrices
    uint8_t(*current)[n_loc_c] = malloc(n_loc_r * sizeof(*current));
    uint8_t(*next)[n_loc_c] = malloc(n_loc_r * sizeof(*next));

    // Initialize the matrix with random data
    srand(seed);
    fill_matrix(n_loc_r, n_loc_c, current, n, density, 0, 0); // m_offset_r and m_offset_c are 0 for single process

    // Simulation without debug
    double start_time = MPI_Wtime();

    // Set the update function depending on the debug flag
    // Simulation loop
    for (int gen = 0; gen < iterations; gen++)
    {
        //update_matrix(n_loc_r, n_loc_c, current, next);
        update_matrix(n_loc_r, n_loc_c, current, next);
        // Swap pointers
        uint8_t(*temp)[n_loc_c] = current;
        current = next;
        next = temp;
    }

    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;

    // If debug is enabled, repeat with debug output
    if (debug == 1)
    {
        printf("Debugging---------------------------\n");

        // Reinitialize the matrix
        srand(seed);
        fill_matrix(n_loc_r, n_loc_c, current, n, density, 0, 0);
        for (int gen = 0; gen < iterations; gen++)
        {
            update_matrix_debug(n_loc_r, n_loc_c, current, next, rank, size, gen);
            uint8_t(*temp)[n_loc_c] = current;
            current = next;
            next = temp;
        }
        printf("------------------------------------\n");
    }

    // Count alive and dead cells
    int local_alive, local_dead, total_alive, total_dead;
    count_cells(n_loc_r, n_loc_c, current, &local_alive, &local_dead);
    MPI_Reduce(&local_alive, &total_alive, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_dead, &total_dead, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        printf("Alive cells: %d, Dead cells: %d\n", total_alive, total_dead);
        printf("Computation time: %f ms\n\n", elapsed_time * 1000);
    }

    // Cleanup
    free(current);
    free(next);

    MPI_Finalize();
}

// mpicc -Wall -I./include -o main main.c src/simulation_control.c src/matrix_operations.c src/utilities.c
// mpirun -np 1 ./main -n 10 -s 2 -d 20 -v -i 5
