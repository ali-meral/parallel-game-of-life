#ifndef MPI_COMMUNICATION_H
#define MPI_COMMUNICATION_H
#include <mpi.h>
#include <stdlib.h>
#include <stdint.h>

void exchange_boundaries(int n_loc_r, int n_loc_c, uint8_t (*current)[n_loc_c], uint8_t (*extended_matrix)[n_loc_c + 4], MPI_Comm cartcomm);
#endif