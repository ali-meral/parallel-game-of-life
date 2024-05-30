#ifndef MATRIX_OPERATIONS_H
#define MATRIX_OPERATIONS_H

#include <stdint.h>

int matrices_are_equal(int n, uint8_t (*matrix1)[n], uint8_t (*matrix2)[n]);
void fill_matrix(int n_loc_r, int n_loc_c, uint8_t (*matrix)[n_loc_c], int n, int density, int m_offset_r, int m_offset_c);
void fill_extended_grid(int n_loc_r, int n_loc_c, uint8_t (*current)[n_loc_c], uint8_t (*extended)[n_loc_c + 4]);
void print_matrix(int n_loc_r, int n_loc_c, uint8_t (*matrix)[n_loc_c], int rank, int size);
void print_matrix_seq(int n_loc_r, int n_loc_c, uint8_t (*matrix)[n_loc_c]);
void update_matrix(int n_loc_r, int n_loc_c, uint8_t (*extended_matrix)[n_loc_c + 4], uint8_t (*next)[n_loc_c]);

#endif
