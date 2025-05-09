#ifndef MPI_COMMUNICATION_H
#define MPI_COMMUNICATION_H
#include <mpi.h>
#include <stdlib.h>
#include <stdint.h>

void collectives_communicate(int n_loc_r, int n_loc_c, uint8_t (*matrix)[n_loc_c], uint8_t (*extended_matrix)[n_loc_c + 4], MPI_Comm cartcomm);
void gather_ranks(int n_loc_r, int n_loc_c, int n, int rank, int size, int m_offset_r, int m_offset_c, 
                    uint8_t (*matrix)[n_loc_c], uint8_t (*global_matrix_par)[n], MPI_Comm cartcomm);
void sendrecv_communicate(int n_loc_r, int n_loc_c, uint8_t (*matrix)[n_loc_c], uint8_t (*extended_matrix)[n_loc_c + 4], MPI_Comm cartcomm);
#endif