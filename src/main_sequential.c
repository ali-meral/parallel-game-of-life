#include "matrix_operations.h"
#include "utilities.h"
#include "game.h"
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

void main_sequential(int argc, char *argv[])
{
    // default values
    int seed = 10;
    int n = 10;
    int verbose = 0;
    int density = 10;
    int iterations = 3;

    // Parse command-line arguments
    parse_arguments(argc, argv, &n, &seed, &density, &iterations, &verbose, NULL);

    uint8_t(*final_matrix)[n] = (uint8_t(*)[n])malloc(n * n * sizeof(uint8_t));

    // time variable that run sequential simulation function will modify
    // to store the time it took to run the simulation
    double seq_time;
    run_sequential_simulation(n, seed, density, iterations, final_matrix, &seq_time);

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
    printf("%d\t%d\t%d\t%d\tsequential\t%lf\t%d\t%d\n", n, seed, density, iterations, seq_time, alive, dead);


    free(final_matrix);
}

int main(int argc, char *argv[])
{
    main_sequential(argc, argv);
    return 0;
}