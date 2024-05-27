#ifndef RUN_SEQUENTIAL_H
#define RUN_SEQUENTIAL_H

#include <stdint.h>

void run_sequential_simulation(int n, int seed, int density, int iterations, uint8_t (*final_matrix)[n]);
void main_sequential(int argc, char *argv[]);

#endif
