#include "run_sequential.h"
#include "matrix_operations.h"
#include "utilities.h"
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

void run_sequential_simulation(int n, int seed, int density, int iterations, uint8_t (*final_matrix)[n])
{
    uint8_t(*matrix)[n] = (uint8_t(*)[n])malloc(n * n * sizeof(uint8_t));
    uint8_t(*next_matrix)[n] = (uint8_t(*)[n])malloc(n * n * sizeof(uint8_t));

    srand(seed);
    fill_matrix(n, n, matrix, n, density, 0, 0);


    // print intial matrix if verbose
    printf("Initial matrix:\n");
    print_matrix_seq(n, n, matrix);

    // Time with MPI_Wtime
    double start_time_seq, end_time_seq, elapsed_time_seq;

    start_time_seq = MPI_Wtime();
    for (int gen = 1; gen <= iterations; gen++)
    {
        update_matrix(n, n, matrix, next_matrix);

        uint8_t(*temp)[n] = matrix;
        matrix = next_matrix;
        next_matrix = temp;

        // Uncomment to print matrix each generation
        // printf("Generation %d:\n", gen);
        // print_matrix_seq(n, n, matrix);

    }

    end_time_seq = MPI_Wtime();
    elapsed_time_seq = end_time_seq - start_time_seq;

    printf("Sequential Computation time: %f ms\n", elapsed_time_seq * 1000);

    // Copy final matrix to the provided pointer
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            final_matrix[i][j] = matrix[ i][j];
        }
    }

    free(matrix);
    free(next_matrix);
}

void main_sequential(int argc, char *argv[])
{

    // default values
    int seed = 10;
    int n = 10;
    int verbose = 0;
    int density = 10;
    int iterations = 3;

    // Parse command-line arguments
    // NULL because we don't need to verify
    parse_arguments(argc, argv, &n, &seed, &density, &iterations, &verbose, NULL); 

    uint8_t(*final_matrix)[n] = (uint8_t(*)[n])malloc(n * n * sizeof(uint8_t));


    // Time with MPI_Wtime
    run_sequential_simulation(n, seed, density, iterations, final_matrix);

    // print final matrix if verbose
    if (verbose)
    {
        printf("Final matrix:\n");
        print_matrix_seq(n, n, final_matrix);
    }


    // count cells using count_cells
    int alive = 0;
    int dead = 0;
    count_cells(n, n, final_matrix, &alive, &dead);
    printf("alive: %d, dead: %d\n", alive, dead);
    

    free(final_matrix);

}