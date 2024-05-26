#ifndef MATRIX_OPERATIONS_H
#define MATRIX_OPERATIONS_H

#include <stdint.h>

void fill_matrix(int n_loc_r, int n_loc_c, uint8_t (*matrix)[n_loc_c], int n, int density, int m_offset_r, int m_offset_c);
void print_matrix(int n_loc_r, int n_loc_c, uint8_t (*matrix)[n_loc_c], int rank, int size);
void update_matrix_w_modulus(int n_loc_r, int n_loc_c, uint8_t (*current)[n_loc_c], uint8_t (*next)[n_loc_c]);
void update_matrix_mpi(int n_loc_r, int n_loc_c, uint8_t (*extended_matrix)[n_loc_c + 4], uint8_t (*next)[n_loc_c]);
void update_matrix(int n_loc_r, int n_loc_c, uint8_t (*current)[n_loc_c], uint8_t (*next)[n_loc_c]);
void print_matrix_debug(int n_loc_r, int n_loc_c, uint8_t (*matrix)[n_loc_c], int rank, const char *matrix_name);

#endif
