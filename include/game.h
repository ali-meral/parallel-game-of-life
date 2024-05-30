#ifndef GAME_H
#define GAME_H

#include <stdint.h>

void run_collectives(int argc, char *argv[]);
void run_sendrecv(int argc, char *argv[]);
void run_sequential_simulation(int n, int seed, int density, int iterations, uint8_t (*final_matrix)[n], double *time);

#endif
