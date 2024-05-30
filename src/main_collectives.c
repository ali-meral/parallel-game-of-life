#include "matrix_operations.h"
#include "game.h"
#include "mpi_communication.h"
#include "utilities.h"
#include <mpi.h>

int main(int argc, char *argv[]) {
    run_collectives(argc, argv);
    return 0;
}
