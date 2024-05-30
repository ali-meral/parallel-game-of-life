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
    int dims[2], coords[2], rank, size, prow_idx, pcol_idx;
    int n_loc_r, n_loc_c, nprows, npcols;
    MPI_Comm cartcomm;
    int pers[2];
} SimulationParams;

void initialize(int argc, char *argv[], SimulationParams *params, int reorder)
{
    params->n = params->seed = params->density = params->iterations = params->verbose = params->verify = 0;
    params->rank = params->size = params->prow_idx = params->pcol_idx = params->n_loc_r = params->n_loc_c = params->nprows = params->npcols = 0;
    params->dims[0] = params->dims[1] = params->coords[0] = params->coords[1] = 0;
    params->pers[0] = params->pers[1] = 1;

    MPI_Init(&argc, &argv);
    parse_arguments(argc, argv, &params->n, &params->seed, &params->density, &params->iterations, &params->verbose, &params->verify);
    MPI_Comm_rank(MPI_COMM_WORLD, &params->rank);
    MPI_Comm_size(MPI_COMM_WORLD, &params->size);
    MPI_Dims_create(params->size, 2, params->dims);

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
    case 6:
        params->dims[0] = 1;
        params->dims[1] = 6;
        break;
    }

    if (params->verbose && params->rank == 0)
    {
        printf("n: %d\ns: %d\nDimensions created: [%d, %d]\n", params->n, params->seed, params->dims[0], params->dims[1]);
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

    if (params->verbose)
    {
        printf("Rank %d: prow_idx: %d pcol_idx: %d m_offset_r: %d m_offset_c: %d\n", params->rank, params->prow_idx, params->pcol_idx, params->prow_idx * params->n_loc_r, params->pcol_idx * params->n_loc_c);
    }
}

void run_simulation(SimulationParams *params, void (*communicate)(int, int, uint8_t(*)[params->n_loc_c], uint8_t(*)[params->n_loc_c], MPI_Comm), const char *computation_label)
{
    uint8_t(*current)[params->n_loc_c] = malloc(params->n_loc_r * params->n_loc_c * sizeof(uint8_t));
    uint8_t(*next)[params->n_loc_c] = malloc(params->n_loc_r * params->n_loc_c * sizeof(uint8_t));
    uint8_t(*global_matrix_par)[params->n] = malloc(params->n * params->n * sizeof(uint8_t));
    int extended_r = params->n_loc_r + 4, extended_c = params->n_loc_c + 4;
    uint8_t(*extended_matrix)[extended_c] = malloc(extended_r * extended_c * sizeof(uint8_t));

    srand(params->seed);
    fill_matrix(params->n_loc_r, params->n_loc_c, current, params->n, params->density, params->prow_idx * params->n_loc_r, params->pcol_idx * params->n_loc_c);

    if (params->verbose)
    {
        printf("Generation 0:\n");
        print_matrix(params->n_loc_r, params->n_loc_c, current, params->rank, params->size);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    for (int gen = 1; gen <= params->iterations; gen++)
    {
        communicate(params->n_loc_r, params->n_loc_c, current, extended_matrix, params->cartcomm);
        update_matrix(params->n_loc_r, params->n_loc_c, extended_matrix, next);

        uint8_t(*temp)[params->n_loc_c] = current;
        current = next;
        next = temp;
    }

    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;
    double max_time;
    MPI_Reduce(&elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, params->cartcomm);

    if (params->rank == 0)
    {
        printf("%s Computation time: %f ms\n", computation_label, max_time * 1000);
    }

    gather_ranks(params->n_loc_r, params->n_loc_c, params->n, params->rank, params->size, params->prow_idx * params->n_loc_r, params->pcol_idx * params->n_loc_c, current, global_matrix_par, params->cartcomm, params->verbose, params->iterations);

    if (params->rank == 0 && params->verbose)
    {
        printf("Parallel Final Generation %d:\n", params->iterations);
        print_matrix(params->n, params->n, (uint8_t(*)[params->n])global_matrix_par, params->rank, params->size);
    }

    if (params->rank == 0)
    {
        int local_alive = 0, local_dead = 0;
        count_cells(params->n, params->n, global_matrix_par, &local_alive, &local_dead);
        printf("alive: %d, dead: %d\n", local_alive, local_dead);
    }

    if (params->verify && params->rank == 0)
    {
        uint8_t(*final_matrix)[params->n] = malloc(params->n * params->n * sizeof(uint8_t));
        run_sequential_simulation(params->n, params->seed, params->density, params->iterations, final_matrix);
        int result = compare_matrices(params->n, final_matrix, global_matrix_par);
        printf("Parallel and sequential computed matrices %s.\n", result ? "match" : "do not match");
        free(final_matrix);
    }

    free(current);
    free(next);
    free(global_matrix_par);
    free(extended_matrix);
}

void run_collectives(int argc, char *argv[])
{
    SimulationParams params;
    initialize(argc, argv, &params, 0);
    run_simulation(&params, collectives_communicate, "Collectives");
    MPI_Finalize();
}

void run_parallel(int argc, char *argv[])
{
    SimulationParams params;
    initialize(argc, argv, &params, 0);

    MPI_Datatype col_type;
    MPI_Type_vector(params.n_loc_r, 1, params.n_loc_c, MPI_UINT8_T, &col_type);
    MPI_Type_commit(&col_type);

    run_simulation(&params, sendrecv_communicate, "Parallel");

    MPI_Type_free(&col_type);
    MPI_Finalize();
}

void run_sequential_simulation(int n, int seed, int density, int iterations, uint8_t (*final_matrix)[n])
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
    printf("Sequential Computation time: %f ms\n", (end_time_seq - start_time_seq) * 1000);

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
