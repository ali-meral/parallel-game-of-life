#include "utilities.h"
#include "matrix_operations.h"
#include <getopt.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void count_cells(int n_loc_r, int n_loc_c, uint8_t matrix[n_loc_r][n_loc_c], int *alive_count, int *dead_count)
{
    *alive_count = 0;
    *dead_count = 0;
    for (int i = 0; i < n_loc_r; i++)
    {
        for (int j = 0; j < n_loc_c; j++)
        {
            if (matrix[i][j] == 1)
            {
                (*alive_count)++;
            }
            else
            {
                (*dead_count)++;
            }
        }
    }
}

int compare_matrices(int n, uint8_t (*matrix1)[n], uint8_t (*matrix2)[n])
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {

            if (matrix1[i][j] != matrix2[i][j])
            {
                printf("matrix1[%d][%d] = %d, matrix2[%d][%d] = %d\n", i, j, matrix1[i][j], i, j, matrix2[i][j]);
                return 0; // Matrices are not equal
            }
        }
    }
    return 1; // Matrices are equal
}

// Function to parse command-line arguments
void parse_arguments(int argc, char *argv[], int *n, int *seed, int *density, int *iterations, int *verbose, int *verify)
{
    static struct option long_options[] = {
        {"number", optional_argument, 0, 'n'},
        {"seed", optional_argument, 0, 's'},
        {"density", optional_argument, 0, 'd'},
        {"verbose", optional_argument, 0, 'v'},
        {"iterations", optional_argument, 0, 'i'},
        {"verify", optional_argument, 0, 'c'},
        {0, 0, 0, 0}};

    int opt;
    while ((opt = getopt_long(argc, argv, "n:s:d:vi:c", long_options, NULL)) != -1)
    {
        switch (opt)
        {
        case 'n':
            *n = atoi(optarg);
            break;
        case 's':
            *seed = atoi(optarg);
            break;
        case 'd':
            *density = atoi(optarg);
            break;
        case 'i':
            *iterations = atoi(optarg);
            break;
        case 'v':
            *verbose = 1;
            break;
        case 'c':
            *verify = 1;
            break;
        default:
            fprintf(stderr, "Usage: %s [-n size] [-s seed] [-d density] [-i iterations] [-v verbose]\n", argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    if (*density <= 0 || *density > 100)
    {
        fprintf(stderr, "Density should be between 1 and 100\n");
        exit(EXIT_FAILURE);
    }

    // check for the --sequential flag
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "--sequential") == 0)
        {
            continue;
        }
    }
}

void reorder_and_compare_communicators(MPI_Comm cartcomm, int *dims, int rank, int verbose)
{
    int reorder = 1;
    MPI_Comm cartcomm_reorder;
    int pers[2] = {0, 0};

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, pers, reorder, &cartcomm_reorder);

    int result;
    MPI_Comm_compare(cartcomm, cartcomm_reorder, &result);

    if (rank == 0 && verbose == 1)
    {
        printf("After reordering the communicators are ");
        if (result == MPI_IDENT)
        {
            printf("identical.\n");
        }
        else if (result == MPI_CONGRUENT)
        {
            printf("congruent.\n");
        }
        else if (result == MPI_SIMILAR)
        {
            printf("similar.\n");
        }
        else if (result == MPI_UNEQUAL)
        {
            printf("unequal.\n");
        }
    }
}
