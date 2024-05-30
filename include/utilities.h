#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdint.h>
#include <mpi.h>

void count_cells(int n_loc_r, int n_loc_c, uint8_t matrix[n_loc_r][n_loc_c], int *alive_count, int *dead_count);
void parse_arguments(int argc, char *argv[], int *n, int *seed, int *density, int *iterations, int *verbose, int *verify);
void reorder_and_compare_communicators(MPI_Comm cartcomm, int *dims, int rank, int verbose);
int wrap(int idx, int limit);

#endif
