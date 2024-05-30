#include "utilities.h"
#include <getopt.h>
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
        printf("Communicator after reordering are ");
        if (result == MPI_IDENT)
        {
            printf("MPI_IDENT.\n");
        }
        else if (result == MPI_CONGRUENT)
        {
            printf("MPI_CONGRUENT.\n");
        }
        else if (result == MPI_SIMILAR)
        {
            printf("MPI_SIMILAR.\n");
        }
        else if (result == MPI_UNEQUAL)
        {
            printf("MPI_UNEQUAL.\n");
        }
    }
}

int wrap(int idx, int limit)
{
  return idx + limit * (idx < 0) - limit * (idx >= limit);
}