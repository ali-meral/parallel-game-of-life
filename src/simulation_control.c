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
    MPI_Comm cartcomm_reorder;

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

    // Reordering of ranks
    // Create Cartesian Communicator with Reordering
    reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, pers, reorder, &cartcomm_reorder);
    // Compare Communicators
    int result;
    MPI_Comm_compare(cartcomm, cartcomm_reorder, &result);

    // TODO Uncomment
    // if (rank == 0) {
    //     printf("After reordering the communicators are ");
    //     if (result == MPI_IDENT) {
    //         printf("identical.\n");
    //     } else if (result == MPI_CONGRUENT) {
    //         printf("congruent.\n");
    //     } else if (result == MPI_SIMILAR) {
    //         printf("similar.\n");
    //     } else if (result == MPI_UNEQUAL) {
    //         printf("unequal.\n");
    //     }
    // }

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

    // get neighbors of the current rank
    // Determine the neighbors of the current process
    int left, right, up, down;
    MPI_Cart_shift(cartcomm, 1, 1, &left, &right);
    MPI_Cart_shift(cartcomm, 0, 1, &up, &down);

    srand(seed);

    int m_offset_r = prow_idx * n_loc_r;
    int m_offset_c = pcol_idx * n_loc_c;
    if (verbose == 1)
    {
        printf("Rank %d: prow_idx: %d pcol_idx: %d m_offset_r: %d m_offset_c: %d\n", rank, prow_idx, pcol_idx, m_offset_r, m_offset_c);
    }

    fflush(stdout);

    // Allocate two matrices for current and next states
    uint8_t(*current)[n_loc_c] = (uint8_t(*)[n_loc_c])malloc(n_loc_r * n_loc_c * sizeof(uint8_t));
    uint8_t(*next)[n_loc_c] = (uint8_t(*)[n_loc_c])malloc(n_loc_r * n_loc_c * sizeof(uint8_t));

    // for sequential implementation
    uint8_t(*global_matrix_seq)[n] = (uint8_t(*)[n])malloc(n * n * sizeof(uint8_t));
    uint8_t(*next_global_matrix_seq)[n] = (uint8_t(*)[n])malloc(n * n * sizeof(uint8_t));

    // for parallel implementation, just to gather the results
    uint8_t(*global_matrix_par)[n] = (uint8_t(*)[n])malloc(n * n * sizeof(uint8_t));

    fill_matrix(n_loc_r, n_loc_c, current, n, density, m_offset_r, m_offset_c);

    // print if verbose
    if (verbose == 1)
    {
        printf("Generation 0:\n");
        print_matrix(n_loc_r, n_loc_c, current, rank, size);
    }

    // SEQUENTIAL IMPLEMENTATION START ==========================================================
    // Allocate a global matrix for the root process
    // fill the global matrix using the same seed with the fill function without gather
    srand(seed);
    fill_matrix(n, n, global_matrix_seq, n, density, 0, 0);
    //  print global matrix
    if (rank == 0 && verbose == 1)
    {
        printf("Sequential Global Generation 0:\n");
        print_matrix(n, n, global_matrix_seq, rank, size);
    }

    // Sequential update loop
    for (int gen = 1; gen < iterations + 1; gen++)
    {

        update_matrix(n, n, global_matrix_seq, next_global_matrix_seq);

        uint8_t(*temp)[n] = global_matrix_seq;
        global_matrix_seq = next_global_matrix_seq;
        next_global_matrix_seq = temp;

        if (rank == 0 && verbose == 1)
        {
            printf("Sequential Global Generation %d:\n", gen);
            print_matrix(n, n, global_matrix_seq, rank, size);
        }
    }

    // SEQUENTIAL IMPLEMENTATION END ==========================================================

    // Start timing
    double start_time, end_time, elapsed_time;
    MPI_Barrier(MPI_COMM_WORLD); // Synchronize before starting the timer
    start_time = MPI_Wtime();

    for (int gen = 1; gen < iterations + 1; gen++)
    {
        update_matrix(n_loc_r, n_loc_c, current, next);

        // uncomment to print for each process
        // if (verbose == 1)
        // {
        //     printf("Generation %d:\n", gen + 1);
        //     print_matrix(n_loc_r, n_loc_c, next, rank, size);
        // }

        // Swap pointers for the next iteration
        uint8_t(*temp)[n_loc_c] = current;
        current = next;
        next = temp;
    }

    end_time = MPI_Wtime(); // Stop the timer before post-processing
    elapsed_time = end_time - start_time;

    // Gather the submatrices to the root process
    MPI_Gather(current, n_loc_r * n_loc_c, MPI_UINT8_T, global_matrix_par, n_loc_r * n_loc_c, MPI_UINT8_T, 0, MPI_COMM_WORLD);

    if (rank == 0 && verbose == 1)
    {
        printf("Parallel Final Generation %d:\n", iterations);
        print_matrix(n, n, (uint8_t(*)[n])global_matrix_par, rank, size);
    }

    // int local_alive = 0;
    // int local_dead = 0;
    // count_cells(n_loc_r, n_loc_c, current, &local_alive, &local_dead);
    // count instead for the global matrix

    // Count cells in the global matrix at the root process
    if (rank == 0)
    {
        int local_alive = 0;
        int local_dead = 0;
        count_cells(n, n, global_matrix_par, &local_alive, &local_dead);

        printf("Global Alive cells: %d, Global Dead cells: %d\n", local_alive, local_dead);
        printf("Computation time: %f ms\n\n", elapsed_time * 1000);
    }

    // fflush(stdout);
    // MPI_Barrier(MPI_COMM_WORLD);


    free(current);
    free(next);
    free(global_matrix_seq);
    free(next_global_matrix_seq);
    free(global_matrix_par);

    MPI_Finalize();
}

// mpicc -Wall -I./include -o main main.c src/simulation_control.c src/matrix_operations.c src/utilities.c
// mpirun -np 1 ./main -n 10 -s 2 -d 20 -v -i 5