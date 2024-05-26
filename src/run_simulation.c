#include "run_simulation.h"
#include "matrix_operations.h"
#include "mpi_communication.h"
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
    int verbose = 0;
    int density = 10;   // in percent
    int iterations = 3; // default number of iterations
    int verify = 0;

    int n_loc_r, n_loc_c;
    int nprows, npcols;
    int prow_idx, pcol_idx;

    // Cartesian communicator stuff
    int dims[2] = {0, 0};
    int pers[2] = {0, 0};
    int coords[2];
    MPI_Comm cartcomm;

    MPI_Init(&argc, &argv);

    // Parse command-line arguments
    parse_arguments(argc, argv, &n, &seed, &density, &iterations, &verbose, &verify);

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

    // Reorder and compare communicators
    reorder_and_compare_communicators(cartcomm, dims, rank, verbose);

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

    ///////NEW////////////////
    // get neighbors of the current rank
    // Determine the neighbors of the current process
    int left, right, up, down;
    MPI_Cart_shift(cartcomm, 1, 1, &left, &right);
    MPI_Cart_shift(cartcomm, 0, 1, &up, &down);

    ///////NEW////////////////

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

    // for parallel implementation, just to gather the results
    uint8_t(*global_matrix_par)[n] = (uint8_t(*)[n])malloc(n * n * sizeof(uint8_t));

    // Allocate extended matrix for boundary communication
    int extended_r = n_loc_r + 4;
    int extended_c = n_loc_c + 4;
    uint8_t(*extended_matrix)[extended_c] = (uint8_t(*)[extended_c])malloc(extended_r * extended_c * sizeof(uint8_t));


    fill_matrix(n_loc_r, n_loc_c, current, n, density, m_offset_r, m_offset_c);

    // print if verbose
    if (verbose == 1)
    {
        printf("Generation 0:\n");
        print_matrix(n_loc_r, n_loc_c, current, rank, size);
    }

    ///////NEW////////////////
    // Create MPI data type for column communication
    MPI_Datatype col_type;
    MPI_Type_vector(n_loc_r, 1, n_loc_c, MPI_UINT8_T, &col_type);
    MPI_Type_commit(&col_type);
    ///////NEW////////////////

    // Start timing
    double start_time, end_time, elapsed_time;
    MPI_Barrier(MPI_COMM_WORLD); // Synchronize before starting the timer
    start_time = MPI_Wtime();

    // Run the simulation
    for (int gen = 1; gen < iterations + 1; gen++)
    {
        exchange_boundaries(n_loc_r, n_loc_c, current, extended_matrix, cartcomm); // Exchange boundaries

        // print extended matrix if verbose
        if (verbose == 1)
        {
            print_matrix_debug(extended_r, extended_c, extended_matrix, rank, "Extended Matrix");
        }
        
        // update_matrix(n_loc_r, n_loc_c, current, next);
        // update_matrix_w_modulus(n_loc_r, n_loc_c, current, next);

        update_matrix_mpi(n_loc_r, n_loc_c, extended_matrix, next);


        uint8_t(*temp)[n_loc_c] = current;
        current = next;
        next = temp;

        // print if verbose
        if (verbose == 1)
        {
            printf("NotGathered %d:\n", gen);
            print_matrix(n_loc_r, n_loc_c, current, rank, size);
        }
    }

    end_time = MPI_Wtime();
    elapsed_time = end_time - start_time;

    // Gather the submatrices to the root process
    // TODO gather the results to the root process

    if (rank == 0 && verbose == 1)
    {
        printf("Parallel Final Generation %d:\n", iterations);
        print_matrix(n, n, (uint8_t(*)[n])global_matrix_par, rank, size);
    }

    // Count cells in the global matrix at the root process
    if (rank == 0)
    {
        int local_alive = 0;
        int local_dead = 0;
        count_cells(n, n, global_matrix_par, &local_alive, &local_dead);

        printf("Global Alive cells: %d, Global Dead cells: %d\n", local_alive, local_dead);
        // rank 0 prints the computation time for all the processes
        printf("Computation time: %f ms\n\n", elapsed_time * 1000);
    }

    // Compare matrices before freeing memory
    if (verify == 1 && rank == 0)
    {
        uint8_t(*final_matrix)[n] = (uint8_t(*)[n])malloc(n * n * sizeof(uint8_t));
        run_sequential_simulation(n, seed, density, iterations, final_matrix);

        if (verbose == 1)
        {
            printf("Sequential Final Generation %d:\n", iterations);
            print_matrix(n, n, final_matrix, rank, size);
        }

        int result = compare_matrices(n, final_matrix, global_matrix_par);
        if (result == 1)
        {
            printf("Parallel and sequential computed matrices match.\n");
        }
        else
        {
            printf("Parallel and sequential computed matrices do not match.\n");
        }

        free(final_matrix);
    }

    // Free memory
    free(current);
    free(next);
    if (rank == 0)
    {
        free(global_matrix_par);
    }
    free(extended_matrix);
    
    MPI_Finalize();
}

// mpicc -Wall -I./include -o main main.c src/simulation_control.c src/matrix_operations.c src/utilities.c
// mpirun -np 1 ./main -n 10 -s 2 -d 20 -v -i 5