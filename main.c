#include "matrix_operations.h"
#include "run_parallel.h"
#include "run_sequential.h"
#include "mpi_communication.h"
#include "utilities.h"
#include <mpi.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char *argv[]) {
    int sequential = 0;

    // Check if --sequential argument is provided
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--sequential") == 0) {
            sequential = 1;
            break;
        }
    }

    if (sequential) {
        // Run sequential simulation
        run_sequential(argc, argv);
    } else {
        run_parallel(argc, argv);
    }

    return 0;
}
