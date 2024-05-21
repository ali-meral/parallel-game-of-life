#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdint.h>

void count_cells(int n_loc_r, int n_loc_c, uint8_t matrix[n_loc_r][n_loc_c], int *alive_count, int *dead_count);
int compare_matrices(int n, uint8_t (*matrix1)[n], uint8_t (*matrix2)[n]);
void parse_arguments(int argc, char *argv[], int *n, int *seed, int *density, int *iterations, int *verbose, int *verify);
void run_sequential_simulation(int n, int seed, int density, int iterations, uint8_t (*final_matrix)[n]);

#endif
