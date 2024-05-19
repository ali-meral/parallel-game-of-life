#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdint.h>

void count_cells(int n_loc_r, int n_loc_c, uint8_t matrix[n_loc_r][n_loc_c], int *alive_count, int *dead_count);

#endif
