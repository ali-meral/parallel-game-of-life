#include "run_collectives.h"
#include "run_sequential.h"
#include "matrix_operations.h"
#include "mpi_communication.h"
#include "utilities.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

int c(int x, int y, int rows)
{
    return y * rows + x;
}

void run_collectives(int argc, char *argv[])
{
    int rank, size;         // rank and size of MPI communicator
    int n_loc_r, n_loc_c;   // local dimensions of the matrix
    int nprows, npcols;     // number of rows and columns in the process grid
    int prow_idx, pcol_idx; // row and column index of the current process in the process grid

    // default values
    int seed = 1;
    int n = 10;
    int density = 30;
    int iterations = 3;
    int verbose = 0;

    int verify = 0;

    MPI_Init(&argc, &argv);

    int dims[2] = {0, 0};
    int pers[2] = {1, 1};
    int coords[2];

    MPI_Comm cartcomm;

    // Parse command-line arguments
    parse_arguments(argc, argv, &n, &seed, &density, &iterations, &verbose, &verify);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Dims_create(size, 2, dims);
    if (verbose == 1 && rank == 0)
    {
        printf("n: %d\n", n);
        printf("s: %d\n", seed);
        printf("Dimensions created: [%d, %d]\n", dims[0], dims[1]);
    }

    int reorder = 0;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, pers, reorder, &cartcomm);
    MPI_Cart_coords(cartcomm, rank, 2, coords);
    prow_idx = coords[0];
    pcol_idx = coords[1];

    // Reorder and compare communicators
    reorder_and_compare_communicators(cartcomm, dims, rank, verbose);

    nprows = dims[0];
    npcols = dims[1];
    if (n % nprows != 0 || n % npcols != 0)
    {
        if (rank == 0)
        {
            fprintf(stderr, "n should be divisible by nprows and npcols\n");
        }
        exit(EXIT_FAILURE);
    }

    n_loc_r = n / nprows;
    n_loc_c = n / npcols;
    if (verbose == 1 && rank == 0)
    {
        printf("n_loc_r: %d n_loc_c: %d\n", n_loc_r, n_loc_c);
    }

    int m_offset_r = prow_idx * n_loc_r;
    int m_offset_c = pcol_idx * n_loc_c;
    if (verbose == 1)
    {
        printf("Rank %d: prow_idx: %d pcol_idx: %d m_offset_r: %d m_offset_c: %d\n", rank, prow_idx, pcol_idx, m_offset_r, m_offset_c);
    }

    fflush(stdout);

    // Allocate two matrices for current and next states
    uint8_t(*current)[n_loc_c] = (uint8_t(*)[n_loc_c])malloc(n_loc_r * n_loc_c * sizeof(uint8_t));
    uint8_t(*next)[n_loc_c] = (uint8_t(*)[n_loc_c])malloc(n_loc_r * n_loc_c * sizeof(uint8_t));

    // for parallel implementation, just to gather the results
    uint8_t(*global_matrix_par)[n] = (uint8_t(*)[n])malloc(n * n * sizeof(uint8_t));

    // Allocate extended matrix for boundary communication
    int extended_r = n_loc_r + 4;
    int extended_c = n_loc_c + 4;
    uint8_t(*extended_matrix)[extended_c] = (uint8_t(*)[extended_c])malloc(extended_r * extended_c * sizeof(uint8_t));

    srand(seed);
    fill_matrix(n_loc_r, n_loc_c, current, n, density, m_offset_r, m_offset_c);

    // print if verbose
    if (verbose == 1)
    {
        printf("Generation 0:\n");
        print_matrix(n_loc_r, n_loc_c, current, rank, size);
    }

    // Prepate the extended matrix by copying the current matrix in its center. Rest should be 0

    for (int i = 0; i < extended_r; i++)
    {
        for (int j = 0; j < extended_c; j++)
        {
            if (i >= 2 && i < n_loc_r + 2 && j >= 2 && j < n_loc_c + 2)
            {
                extended_matrix[i][j] = current[i - 2][j - 2];
            }
            else
            {
                extended_matrix[i][j] = 0;
            }
        }
    }

    // COLLECTIVE ======================================================================

    printf("Rank %d, Coordinates: (%d, %d)\n", rank, coords[0], coords[1]);
    MPI_Datatype rowtype_2_layers, coltype_2_layers;

    // Create datatype for two layers of rows
    MPI_Type_vector(2, n_loc_c, n_loc_c + 4, MPI_UINT8_T, &rowtype_2_layers);
    MPI_Type_commit(&rowtype_2_layers);

    // Create datatype for two layers of columns
    MPI_Type_vector(n_loc_r, 2, n_loc_c + 4, MPI_UINT8_T, &coltype_2_layers);
    MPI_Type_commit(&coltype_2_layers);

MPI_Datatype cornertype_4;
int blocklengths[4] = {1, 1, 1, 1};
MPI_Aint displacements[4];
MPI_Datatype types[4] = {MPI_UINT8_T, MPI_UINT8_T, MPI_UINT8_T, MPI_UINT8_T};

displacements[0] = 0;
displacements[1] = 1 * sizeof(uint8_t);
displacements[2] = n_loc_c * sizeof(uint8_t);
displacements[3] = (n_loc_c + 1) * sizeof(uint8_t);

MPI_Type_create_struct(4, blocklengths, displacements, types, &cornertype_4);
MPI_Type_commit(&cornertype_4);

    int neighbors[8]; // Store ranks of neighbors
    int upper_left_coords[2], upper_right_coords[2], lower_left_coords[2], lower_right_coords[2];

    // Determine neighbor coordinates
    upper_left_coords[0] = (coords[0] - 1 + dims[0]) % dims[0];
    upper_left_coords[1] = (coords[1] - 1 + dims[1]) % dims[1];

    upper_right_coords[0] = (coords[0] - 1 + dims[0]) % dims[0];
    upper_right_coords[1] = (coords[1] + 1) % dims[1];

    lower_left_coords[0] = (coords[0] + 1) % dims[0];
    lower_left_coords[1] = (coords[1] - 1 + dims[1]) % dims[1];

    lower_right_coords[0] = (coords[0] + 1) % dims[0];
    lower_right_coords[1] = (coords[1] + 1) % dims[1];

    // Get ranks of diagonal neighbors
    MPI_Cart_rank(cartcomm, upper_left_coords, &neighbors[0]);
    MPI_Cart_rank(cartcomm, upper_right_coords, &neighbors[1]);
    MPI_Cart_rank(cartcomm, lower_left_coords, &neighbors[2]);
    MPI_Cart_rank(cartcomm, lower_right_coords, &neighbors[3]);

    // Get ranks of direct neighbors
    MPI_Cart_shift(cartcomm, 1, 1, &neighbors[4], &neighbors[5]); // left-right
    MPI_Cart_shift(cartcomm, 0, 1, &neighbors[6], &neighbors[7]); // up-down

    // Create distributed graph communicator
    MPI_Comm dist_graph_comm;
    MPI_Dist_graph_create_adjacent(cartcomm, 8, neighbors, MPI_UNWEIGHTED, 8, neighbors, MPI_UNWEIGHTED, MPI_INFO_NULL, 0, &dist_graph_comm);

    // Set up communication buffers
    int sendcounts[8], recvcounts[8];
    MPI_Aint senddispls[8], recvdispls[8];
    MPI_Datatype sendtypes[8], recvtypes[8];

    // Initialize sendcounts and recvcounts
    for (int i = 0; i < 8; i++)
    {
        sendcounts[i] = recvcounts[i] = 1;
    }

    // Calculate and set send displacements using the current matrix
senddispls[0] = ((MPI_Aint)&current[0][0]) - (MPI_Aint)current; 
sendtypes[0] = cornertype_4;

// Top-right corner (first 2x2 block)
senddispls[1] = ((MPI_Aint)&current[0][n_loc_c - 2]) - (MPI_Aint)current; 
sendtypes[1] = cornertype_4;

// Bottom-left corner (last 2x2 block)
senddispls[2] = ((MPI_Aint)&current[n_loc_r - 2][0]) - (MPI_Aint)current; 
sendtypes[2] = cornertype_4;

// Bottom-right corner (last 2x2 block)
senddispls[3] = ((MPI_Aint)&current[n_loc_r - 2][n_loc_c - 2]) - (MPI_Aint)current; 
sendtypes[3] = cornertype_4;

    senddispls[4] = ((MPI_Aint)&current[0][0]) - (MPI_Aint)current;
    sendtypes[4] = coltype_2_layers;

    senddispls[5] = ((MPI_Aint)&current[0][n_loc_c - 2]) - (MPI_Aint)current;
    sendtypes[5] = coltype_2_layers;

    senddispls[6] = ((MPI_Aint)&current[0][0]) - (MPI_Aint)current;
    sendtypes[6] = rowtype_2_layers;

    senddispls[7] = ((MPI_Aint)&current[n_loc_r - 2][0]) - (MPI_Aint)current;
    sendtypes[7] = rowtype_2_layers;

    // Set the receive displacements using the extended matrix
    recvdispls[0] = ((MPI_Aint)&extended_matrix[0][0]) - (MPI_Aint)extended_matrix;
    recvtypes[0] = cornertype_4;

    recvdispls[1] = ((MPI_Aint)&extended_matrix[0][n_loc_c + 2]) - (MPI_Aint)extended_matrix;
    recvtypes[1] = cornertype_4;

    recvdispls[2] = ((MPI_Aint)&extended_matrix[n_loc_r + 2][0]) - (MPI_Aint)extended_matrix;
    recvtypes[2] = cornertype_4;

    recvdispls[3] = ((MPI_Aint)&extended_matrix[n_loc_r + 2][n_loc_c + 2]) - (MPI_Aint)extended_matrix;
    recvtypes[3] = cornertype_4;

    recvdispls[4] = ((MPI_Aint)&extended_matrix[2][0]) - (MPI_Aint)extended_matrix;
    recvtypes[4] = coltype_2_layers;

    recvdispls[5] = ((MPI_Aint)&extended_matrix[2][n_loc_c + 2]) - (MPI_Aint)extended_matrix;
    recvtypes[5] = coltype_2_layers;

    recvdispls[6] = ((MPI_Aint)&extended_matrix[0][2]) - (MPI_Aint)extended_matrix;
    recvtypes[6] = rowtype_2_layers;

    recvdispls[7] = ((MPI_Aint)&extended_matrix[n_loc_r + 2][2]) - (MPI_Aint)extended_matrix;
    recvtypes[7] = rowtype_2_layers;



    // Print the values to be sent
// Debug print to verify the buffers
for (int i = 0; i < 8; i++) {
    printf("Rank %d, senddispls[%d]=%ld, type=", rank, i, senddispls[i]);
    if (sendtypes[i] == cornertype_4) {
        printf("cornertype_4, Values: ");
        uint8_t* ptr = (uint8_t*)((MPI_Aint)current + senddispls[i]);
        printf("%d %d %d %d\n", ptr[0], ptr[1], ptr[n_loc_c], ptr[n_loc_c + 1]);
    } else if (sendtypes[i] == coltype_2_layers) {
        printf("coltype_2_layers, Values: ");
        for (int j = 0; j < n_loc_r; j++) {
            printf("%d %d ", current[j][senddispls[i] / sizeof(uint8_t)], current[j][(senddispls[i] / sizeof(uint8_t)) + 1]);
        }
        printf("\n");
    } else if (sendtypes[i] == rowtype_2_layers) {
        printf("rowtype_2_layers, Values: ");
        int row_index = senddispls[i] / (n_loc_c * sizeof(uint8_t));
        for (int j = 0; j < n_loc_c; j++) {
            printf("%d %d ", current[row_index][j], current[row_index + 1][j]);
        }
        printf("\n");
    }
}


    // Call MPI_Neighbor_alltoallw
    // MPI_Neighbor_alltoallw(
    //     current, sendcounts, senddispls, sendtypes,
    //     extended_matrix, recvcounts, recvdispls, recvtypes,
    //     dist_graph_comm
    // );

    // Print the extended matrix to verify
    // printf("Rank %d, Extended matrix after communication:\n", rank);
    // for (int i = 0; i < n_loc_r + 4; i++) {
    //     for (int j = 0; j < n_loc_c + 4; j++) {
    //         printf("%d ", extended_matrix[i][j]);
    //     }
    //     printf("\n");
    // }

    // Print the values to be sent
    // printf("Rank %d, Sending values:\n", rank);
    // for (int i = 0; i < 8; i++) {
    //     printf("senddispls[%d]=%ld, type=", i, senddispls[i]);
    //     if (sendtypes[i] == cornertype) {
    //         printf("cornertype, Value: %d\n", *(uint8_t*)((MPI_Aint)extended_matrix + senddispls[i]));
    //     } else if (sendtypes[i] == coltype) {
    //         printf("coltype, Values: ");
    //         for (int j = 0; j < n_loc_r; j++) {
    //             printf("%d ", extended_matrix[j + 1][senddispls[i] / sizeof(uint8_t) % extended_c]);
    //         }
    //         printf("\n");
    //     } else if (sendtypes[i] == rowtype) {
    //         printf("rowtype, Values: ");
    //         for (int j = 0; j < n_loc_c; j++) {
    //             printf("%d ", extended_matrix[senddispls[i] / sizeof(uint8_t) / extended_c][j + 1]);
    //         }
    //         printf("\n");
    //     }
    // }

    // print extended matrix for each process
    // block until all processes have reached this point
    MPI_Barrier(cartcomm);
    // debug prints
    for (int r = 0; r < size; r++)
    {
        if (rank == r)
        {
            printf("Rank %d\n", rank);
            for (int i = 0; i < extended_r; i++)
            {
                for (int j = 0; j < extended_c; j++)
                {
                    printf("%d ", extended_matrix[i][j]);
                }
                printf("\n");
            }
        }
        MPI_Barrier(cartcomm);
    }

    //==================================================================================

    // Create MPI data type for column communication

    MPI_Datatype col_type;
    MPI_Type_vector(n_loc_r, 1, n_loc_c, MPI_UINT8_T, &col_type);
    MPI_Type_commit(&col_type);

    // Start timing
    double start_time, end_time, elapsed_time, max_time;
    MPI_Barrier(MPI_COMM_WORLD); // Synchronize before starting the timer
    start_time = MPI_Wtime();

    // Run the simulation
    for (int gen = 1; gen < iterations + 1; gen++)
    {
        exchange_boundaries(n_loc_r, n_loc_c, current, extended_matrix, cartcomm); // Exchange boundaries
        update_matrix_mpi(n_loc_r, n_loc_c, extended_matrix, next);

        uint8_t(*temp)[n_loc_c] = current;
        current = next;
        next = temp;
    }

    end_time = MPI_Wtime();
    elapsed_time = end_time - start_time;
    MPI_Reduce(&elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, cartcomm);

    // print max time and also time of each process
    if (rank == 0)
    {
        printf("Parallel Computation time: %f ms\n", max_time * 1000);
    }

    // now we need to gather
    if (rank == 0)
    {
        // root process
        for (int i = 0; i < n_loc_r; ++i)
        {
            for (int j = 0; j < n_loc_c; ++j)
            {
                global_matrix_par[m_offset_r + i][m_offset_c + j] = current[i][j];
            }
        }

        // other processes
        for (int r = 1; r < size; ++r)
        {
            uint8_t *recv_buffer = (uint8_t *)malloc(n_loc_r * n_loc_c * sizeof(uint8_t));
            MPI_Recv(recv_buffer, n_loc_r * n_loc_c, MPI_UINT8_T, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // get again the coordinates of the process
            MPI_Cart_coords(cartcomm, r, 2, coords); // to get offset of received matrix
            int recv_offset_r = coords[0] * n_loc_r;
            int recv_offset_c = coords[1] * n_loc_c;

            // place received submatrix in the global matrix
            for (int i = 0; i < n_loc_r; ++i)
            {
                for (int j = 0; j < n_loc_c; ++j)
                {
                    global_matrix_par[recv_offset_r + i][recv_offset_c + j] = recv_buffer[i * n_loc_c + j];
                }
            }
            free(recv_buffer);
        }
    }
    else
    {
        // processes need to send to root
        MPI_Send(&(current[0][0]), n_loc_r * n_loc_c, MPI_UINT8_T, 0, 0, MPI_COMM_WORLD);
    }

    if (rank == 0 && verbose == 1)
    {
        printf("Parallel Final Generation %d:\n", iterations);
        print_matrix(n, n, (uint8_t(*)[n])global_matrix_par, rank, size);
    }

    // Count cells in the global matrix at the root process
    if (rank == 0)
    {
        int local_alive = 0;
        int local_dead = 0;
        count_cells(n, n, global_matrix_par, &local_alive, &local_dead);

        printf("alive: %d, dead: %d\n", local_alive, local_dead);
    }

    // Compare matrices before freeing memory
    if (verify == 1 && rank == 0)
    {
        uint8_t(*final_matrix)[n] = (uint8_t(*)[n])malloc(n * n * sizeof(uint8_t));

        // Time with MPI_Wtime
        run_sequential_simulation(n, seed, density, iterations, final_matrix);
        int result = compare_matrices(n, final_matrix, global_matrix_par);
        if (result == 1)
        {
            printf("Parallel and sequential computed matrices match.\n");
        }
        else
        {
            printf("Parallel and sequential computed matrices do not match.\n");
        }

        free(final_matrix);
    }

    // Free memory
    free(current);
    free(next);
    free(global_matrix_par);
    free(extended_matrix);

    MPI_Finalize();
}