#include "simulation_control.h"
#include "matrix_operations.h"
#include "utilities.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

void run_simulation(int argc, char *argv[])
{
    int rank, size; // always

    int seed = 10;
    int n = 10; // default values
    int opt;
    int verbose = 0;
    int density = 10;   // in percent
    int iterations = 3; // default number of iterations

    int n_loc_r, n_loc_c;
    int nprows, npcols;
    int prow_idx, pcol_idx;

    // Cartesian communicator stuff
    int dims[2] = {0, 0};
    int pers[2] = {0, 0};
    int coords[2];
    MPI_Comm cartcomm;

    static struct option long_options[] = {{"number", optional_argument, 0, 'n'},
                                           {"seed", optional_argument, 0, 's'},
                                           {"density", optional_argument, 0, 'd'},
                                           {"verbose", optional_argument, 0, 'v'},
                                           {"iterations", optional_argument, 0, 'i'},
                                           {0, 0, 0, 0}};

    MPI_Init(&argc, &argv);

    while (1)
    {
        int option_index = 0;
        opt = getopt_long(argc, argv, "n:s:d:vi:", long_options, &option_index);

        if (opt == -1)
            break;

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
        default:
            fprintf(stderr, "Usage: %s -n <n> -k <k>\n", argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    if (density <= 0 || density > 100)
    {

        fprintf(stderr, "Density should be between 1 and 100\n");
        exit(EXIT_FAILURE);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Dims_create(size, 2, dims);

    if (verbose == 1 && rank == 0)
    {
        printf("n: %d\n", n);
        printf("s: %d\n", seed);
        printf("Dimensions created: [%d, %d]\n", dims[0], dims[1]);
    }

    int reorder = 0;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, pers, reorder, &cartcomm);
    MPI_Cart_coords(cartcomm, rank, 2, coords);
    prow_idx = coords[0];
    pcol_idx = coords[1];

    nprows = dims[0];
    npcols = dims[1];
    if (n % nprows != 0 || n % npcols != 0)
    {
        if (rank == 0)
        {
            fprintf(stderr, "n should be divisible by nprows and npcols\n");
        }
        exit(EXIT_FAILURE);
    }

    n_loc_r = n / nprows;
    n_loc_c = n / npcols;
    if (verbose == 1 && rank == 0)
    {
        printf("n_loc_r: %d n_loc_c: %d\n", n_loc_r, n_loc_c);
    }

    srand(seed);

    int m_offset_r = prow_idx * n_loc_r;
    int m_offset_c = pcol_idx * n_loc_c;
    if (verbose == 1)
    {
        printf("%d: prow_idx: %d pcol_idx: %d m_offset_r: %d m_offset_c: %d\n", rank, prow_idx, pcol_idx, m_offset_r, m_offset_c);
    }

    fflush(stdout);

    // Allocate two matrices for current and next states
    uint8_t(*current)[n_loc_c] = (uint8_t(*)[n_loc_c])malloc(n_loc_r * n_loc_c * sizeof(uint8_t));
    uint8_t(*next)[n_loc_c] = (uint8_t(*)[n_loc_c])malloc(n_loc_r * n_loc_c * sizeof(uint8_t));

    fill_matrix(n_loc_r, n_loc_c, current, n, density, m_offset_r, m_offset_c);

    // print if verbose
    if (verbose == 1)
    {
    printf("Generation 0:\n");
    print_matrix(n_loc_r, n_loc_c, current, rank, size);
    }

    // Start timing
    double start_time, end_time, elapsed_time;
    MPI_Barrier(MPI_COMM_WORLD); // Synchronize before starting the timer
    start_time = MPI_Wtime();

    for (int gen = 0; gen < iterations; gen++)
    {
        update_matrix(n_loc_r, n_loc_c, current, next);
        // to debug uncomment the next file instead

        // Debugging options
        if (verbose == 1)
        {
            printf("Generation %d:\n", gen + 1);
            print_matrix(n_loc_r, n_loc_c, next, rank, size);
        }
        // Swap pointers for the next iteration
        uint8_t(*temp)[n_loc_c] = current;
        current = next;
        next = temp;
    }

    end_time = MPI_Wtime(); // Stop the timer before post-processing
    elapsed_time = end_time - start_time;

    int local_alive = 0;
    int local_dead = 0;
    count_cells(n_loc_r, n_loc_c, current, &local_alive, &local_dead);

    int total_alive, total_dead;
    MPI_Reduce(&local_alive, &total_alive, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_dead, &total_dead, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // fflush(stdout);
    // MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0)
    {
        printf("Alive cells: %d, Dead cells: %d\n", total_alive, total_dead);
        printf("Computation time: %f ms\n\n", elapsed_time * 1000);
    }

    free(current);
    free(next);

    MPI_Finalize();
}

// mpicc -Wall -I./include -o main main.c src/simulation_control.c src/matrix_operations.c src/utilities.c
// mpirun -np 1 ./main -n 10 -s 2 -d 20 -v -i 5
