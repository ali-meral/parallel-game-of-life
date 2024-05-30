#include "game.h"
#include "matrix_operations.h"
#include "mpi_communication.h"
#include "utilities.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

typedef struct {
    int n, seed, density, iterations, verbose, verify;
    int dims[2], coords[2], pers[2], rank, size, prow_idx, pcol_idx;
    int n_loc_r, n_loc_c, nprows, npcols;
    MPI_Comm cartcomm;
} SimulationParams;

void initialize(int argc, char *argv[], SimulationParams *params, int reorder)
{
    params->seed = 10;
    params->n = 10;
    params->density = 10;
    params->iterations = 2;
    params->verbose = 0;
    params->verify = 0;

    params->pers[0] = params->pers[1] = 1; // makes it periodic

    MPI_Init(&argc, &argv);

    parse_arguments(argc, argv, &params->n, &params->seed, &params->density, &params->iterations, &params->verbose, &params->verify);

    MPI_Comm_rank(MPI_COMM_WORLD, &params->rank);
    MPI_Comm_size(MPI_COMM_WORLD, &params->size);
    MPI_Dims_create(params->size, 2, params->dims);

    // Hardcode some matrix sizes for experiments
    switch (params->size)
    {
    case 32:
        params->dims[0] = 1;
        params->dims[1] = 32;
        break;
    case 256:
        params->dims[0] = 8;
        params->dims[1] = 32;
        break;
    case 512:
        params->dims[0] = 16;
        params->dims[1] = 32;
        break;
    case 1024:
        params->dims[0] = 32;
        params->dims[1] = 32;
        break;
    }

    if (params->rank == 0)
    {
        // print parameters grid size, seed, local grid size and dimensions in both directions
        printf("%d\t%d\t%d\t%d\t%d\t%d\t", params->n, params->seed, params->density, params->iterations, params->dims[0], params->dims[1]);
    }


    MPI_Cart_create(MPI_COMM_WORLD, 2, params->dims, params->pers, reorder, &params->cartcomm);
    MPI_Cart_coords(params->cartcomm, params->rank, 2, params->coords);
    params->prow_idx = params->coords[0];
    params->pcol_idx = params->coords[1];
    reorder_and_compare_communicators(params->cartcomm, params->dims, params->rank, params->verbose);

    params->nprows = params->dims[0];
    params->npcols = params->dims[1];

    if (params->n % params->nprows != 0 || params->n % params->npcols != 0)
    {
        if (params->rank == 0)
        {
            fprintf(stderr, "n should be divisible by nprows and npcols\n");
        }
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    params->n_loc_r = params->n / params->nprows;
    params->n_loc_c = params->n / params->npcols;

    if (params->verbose && params->rank == 0)
    {
        printf("n_loc_r: %d n_loc_c: %d\n", params->n_loc_r, params->n_loc_c);
    }

    int m_offset_r = params->prow_idx * params->n_loc_r;
    int m_offset_c = params->pcol_idx * params->n_loc_c;

    if (params->verbose)
    {
        printf("Rank %d: prow_idx: %d pcol_idx: %d m_offset_r: %d m_offset_c: %d\n",
                params->rank, params->prow_idx, params->pcol_idx, m_offset_r, m_offset_c);
    }
}

void run_simulation(SimulationParams *params, void (*communicate)(int, int, uint8_t(*)[params->n_loc_c], uint8_t(*)[params->n_loc_c], MPI_Comm), const char *computation_label)
{
    uint8_t(*matrix)[params->n_loc_c] = malloc(params->n_loc_r * params->n_loc_c * sizeof(uint8_t));
    uint8_t(*next_matrix)[params->n_loc_c] = malloc(params->n_loc_r * params->n_loc_c * sizeof(uint8_t));
    uint8_t(*global_matrix_par)[params->n] = malloc(params->n * params->n * sizeof(uint8_t));
    int extended_r = params->n_loc_r + 4;
    int extended_c = params->n_loc_c + 4;
    uint8_t(*extended_matrix)[extended_c] = malloc(extended_r * extended_c * sizeof(uint8_t));

    srand(params->seed);
    fill_matrix(params->n_loc_r, params->n_loc_c, matrix, params->n, params->density, params->prow_idx * params->n_loc_r, params->pcol_idx * params->n_loc_c);

    if (params->verbose)
    {
        printf("Generation 0:\n");
        print_matrix(params->n_loc_r, params->n_loc_c, matrix, params->rank, params->size);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    for (int gen = 1; gen <= params->iterations; gen++)
    {
        communicate(params->n_loc_r, params->n_loc_c, matrix, extended_matrix, params->cartcomm);
        update_matrix(params->n_loc_r, params->n_loc_c, extended_matrix, next_matrix);

        uint8_t(*temp)[params->n_loc_c] = matrix;
        matrix = next_matrix;
        next_matrix = temp;
    }

    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;
    double max_time;
    MPI_Reduce(&elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, params->cartcomm);

    if (params->rank == 0)
    {
        // computation method and time
        printf("%s\t%lf", computation_label, max_time * 1000);
    }

    gather_ranks(params->n_loc_r,                   // number of local rows
                params->n_loc_c,                    // number of local columns
                params->n,                         // global number of rows and columns
                params->rank,                       // rank of the process
                params->size,                       // total number of processes
                params->prow_idx * params->n_loc_r, // offset of the local rows
                params->pcol_idx * params->n_loc_c, // offset of the local columns
                matrix,                            // local matrix
                global_matrix_par,                  // global matrix
                params->cartcomm                    // cartesian communicator
                );

    // print final generation
    if (params->rank == 0 && params->verbose)
    {
        printf("Parallel Final Generation %d:\n", params->iterations);
        print_matrix(params->n, params->n, (uint8_t(*)[params->n])global_matrix_par, params->rank, params->size);
    }

    // counalive and dead cells and print the result
    if (params->rank == 0)
    {
        int local_alive = 0, local_dead = 0;
        count_cells(params->n, params->n, global_matrix_par, &local_alive, &local_dead);
        printf("\t%d\t%d\n", local_alive, local_dead);
    }

    // check if matrices match (sequential vs parallel)
    double elapsed_time_seq;
    if (params->verify && params->rank == 0)
    {
        uint8_t(*final_matrix)[params->n] = malloc(params->n * params->n * sizeof(uint8_t));
        run_sequential_simulation(params->n, params->seed, params->density, params->iterations, final_matrix, &elapsed_time_seq);
        int result = matrices_are_equal(params->n, final_matrix, global_matrix_par);
        printf("Matrices %s.\n", result ? "match" : "do not match");
        free(final_matrix);
    }

    free(matrix);
    free(next_matrix);
    free(global_matrix_par);
    free(extended_matrix);
}

void run_collectives(int argc, char *argv[])
{
    SimulationParams params;
    initialize(argc, argv, &params, 0);
    run_simulation(&params, collectives_communicate, "collectives");
    MPI_Finalize();
}

void run_sendrecv(int argc, char *argv[])
{
    SimulationParams params;
    initialize(argc, argv, &params, 0);

    MPI_Datatype col_type;
    MPI_Type_vector(params.n_loc_r, 1, params.n_loc_c, MPI_UINT8_T, &col_type);
    MPI_Type_commit(&col_type);

    run_simulation(&params, sendrecv_communicate, "sendrecv");

    MPI_Type_free(&col_type);
    MPI_Finalize();
}

void run_sequential_simulation(int n, int seed, int density, int iterations, uint8_t (*final_matrix)[n], double *time)
{
    uint8_t(*matrix)[n] = malloc(n * n * sizeof(uint8_t));
    uint8_t(*next_matrix)[n] = malloc(n * n * sizeof(uint8_t));

    uint8_t(*extended_matrix)[n + 4] = malloc((n + 4) * (n + 4) * sizeof(uint8_t));

    srand(seed);
    fill_matrix(n, n, matrix, n, density, 0, 0);

    double start_time_seq = MPI_Wtime();
    for (int gen = 1; gen <= iterations; gen++)
    {
        fill_extended_grid(n, n, matrix, extended_matrix);
        update_matrix(n, n, extended_matrix, next_matrix);

        uint8_t(*temp)[n] = matrix;
        matrix = next_matrix;
        next_matrix = temp;
    }

    double end_time_seq = MPI_Wtime();

    *time = end_time_seq - start_time_seq;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)

        {
            final_matrix[i][j] = matrix[i][j];
        }
    }

    free(matrix);
    free(next_matrix);
}
