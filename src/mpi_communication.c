#include "matrix_operations.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

void exchange_boundaries(int n_loc_r, int n_loc_c, uint8_t (*current)[n_loc_c], uint8_t (*extended_matrix)[n_loc_c + 4], MPI_Comm cartcomm) {
    MPI_Status status;
    int rank, coords[2], dims[2], periods[2];
    MPI_Comm_rank(cartcomm, &rank);
    MPI_Cart_get(cartcomm, 2, dims, periods, coords);

    int neighbor_coords[8][2];
    int neighbors[8];

    // Define the relative positions of the 8 neighbors
    int rel_pos[8][2] = {
        {-1, -1}, {-1, 0}, {-1, 1},
        { 0, -1},         { 0, 1},
        { 1, -1}, { 1, 0}, { 1, 1}
    };

    // Calculate the coordinates and ranks of the 8 neighbors
    for (int i = 0; i < 8; i++) {
        neighbor_coords[i][0] = (coords[0] + rel_pos[i][0] + dims[0]) % dims[0];
        neighbor_coords[i][1] = (coords[1] + rel_pos[i][1] + dims[1]) % dims[1];
        MPI_Cart_rank(cartcomm, neighbor_coords[i], &neighbors[i]);
    }

    // Allocate buffers for sending and receiving data
    uint8_t *send_top = (uint8_t *)malloc(2 * n_loc_c * sizeof(uint8_t));
    uint8_t *recv_bottom = (uint8_t *)malloc(2 * n_loc_c * sizeof(uint8_t));

    uint8_t *send_bottom = (uint8_t *)malloc(2 * n_loc_c * sizeof(uint8_t));
    uint8_t *recv_top = (uint8_t *)malloc(2 * n_loc_c * sizeof(uint8_t));

    uint8_t *send_left = (uint8_t *)malloc(2 * n_loc_r * sizeof(uint8_t));
    uint8_t *recv_right = (uint8_t *)malloc(2 * n_loc_r * sizeof(uint8_t));

    uint8_t *send_right = (uint8_t *)malloc(2 * n_loc_r * sizeof(uint8_t));
    uint8_t *recv_left = (uint8_t *)malloc(2 * n_loc_r * sizeof(uint8_t));

    uint8_t *send_top_left = (uint8_t *)malloc(4 * sizeof(uint8_t));
    uint8_t *recv_bottom_right = (uint8_t *)malloc(4 * sizeof(uint8_t));

    uint8_t *send_top_right = (uint8_t *)malloc(4 * sizeof(uint8_t));
    uint8_t *recv_bottom_left = (uint8_t *)malloc(4 * sizeof(uint8_t));

    uint8_t *send_bottom_left = (uint8_t *)malloc(4 * sizeof(uint8_t));
    uint8_t *recv_top_right = (uint8_t *)malloc(4 * sizeof(uint8_t));

    uint8_t *send_bottom_right = (uint8_t *)malloc(4 * sizeof(uint8_t));
    uint8_t *recv_top_left = (uint8_t *)malloc(4 * sizeof(uint8_t));

    // Fill send buffers for top and bottom rows
    for (int j = 0; j < n_loc_c; j++) {
        send_top[j] = current[0][j];
        send_top[j + n_loc_c] = current[1][j];
        send_bottom[j] = current[n_loc_r - 2][j];
        send_bottom[j + n_loc_c] = current[n_loc_r - 1][j];
    }

    // Fill send buffers for left and right columns
    for (int i = 0; i < n_loc_r; i++) {
        send_left[i] = current[i][0];
        send_left[i + n_loc_r] = current[i][1];
        send_right[i] = current[i][n_loc_c - 1];
        send_right[i + n_loc_r] = current[i][n_loc_c - 2];
    }

    // Fill send buffers for corners
    send_top_left[0] = current[0][0];
    send_top_left[1] = current[0][1];
    send_top_left[2] = current[1][0];
    send_top_left[3] = current[1][1];

    send_top_right[0] = current[0][n_loc_c - 1];
    send_top_right[1] = current[0][n_loc_c - 2];
    send_top_right[2] = current[1][n_loc_c - 1];
    send_top_right[3] = current[1][n_loc_c - 2];

    send_bottom_left[0] = current[n_loc_r - 1][0];
    send_bottom_left[1] = current[n_loc_r - 2][0];
    send_bottom_left[2] = current[n_loc_r - 1][1];
    send_bottom_left[3] = current[n_loc_r - 2][1];

    send_bottom_right[0] = current[n_loc_r - 1][n_loc_c - 1];
    send_bottom_right[1] = current[n_loc_r - 2][n_loc_c - 1];
    send_bottom_right[2] = current[n_loc_r - 1][n_loc_c - 2];
    send_bottom_right[3] = current[n_loc_r - 2][n_loc_c - 2];

    // Exchange top and bottom rows
    MPI_Sendrecv(send_top, 2 * n_loc_c, MPI_UINT8_T, neighbors[1], 0,
                 recv_bottom, 2 * n_loc_c, MPI_UINT8_T, neighbors[6], 0,
                 cartcomm, &status);

    MPI_Sendrecv(send_bottom, 2 * n_loc_c, MPI_UINT8_T, neighbors[6], 1,
                 recv_top, 2 * n_loc_c, MPI_UINT8_T, neighbors[1], 1,
                 cartcomm, &status);

    // Exchange left and right columns
    MPI_Sendrecv(send_left, 2 * n_loc_r, MPI_UINT8_T, neighbors[3], 2,
                 recv_right, 2 * n_loc_r, MPI_UINT8_T, neighbors[4], 2,
                 cartcomm, &status);

    MPI_Sendrecv(send_right, 2 * n_loc_r, MPI_UINT8_T, neighbors[4], 3,
                 recv_left, 2 * n_loc_r, MPI_UINT8_T, neighbors[3], 3,
                 cartcomm, &status);

    // Exchange corners
    MPI_Sendrecv(send_top_left, 4, MPI_UINT8_T, neighbors[0], 4,
                 recv_bottom_right, 4, MPI_UINT8_T, neighbors[7], 4,
                 cartcomm, &status);

    MPI_Sendrecv(send_top_right, 4, MPI_UINT8_T, neighbors[2], 5,
                 recv_bottom_left, 4, MPI_UINT8_T, neighbors[5], 5,
                 cartcomm, &status);

    MPI_Sendrecv(send_bottom_left, 4, MPI_UINT8_T, neighbors[5], 6,
                 recv_top_right, 4, MPI_UINT8_T, neighbors[2], 6,
                 cartcomm, &status);

    MPI_Sendrecv(send_bottom_right, 4, MPI_UINT8_T, neighbors[7], 7,
                 recv_top_left, 4, MPI_UINT8_T, neighbors[0], 7,
                 cartcomm, &status);

    // // Print received corners
    // printf("Rank %d received top left corner: ", rank);
    // for (int i = 0; i < 4; i++) {
    //     printf("%d ", recv_top_left[i]);
    // }
    // printf("\n");

    // printf("Rank %d received top right corner: ", rank);
    // for (int i = 0; i < 4; i++) {
    //     printf("%d ", recv_top_right[i]);
    // }
    // printf("\n");

    // printf("Rank %d received bottom left corner: ", rank);
    // for (int i = 0; i < 4; i++) {
    //     printf("%d ", recv_bottom_left[i]);
    // }
    // printf("\n");

    // printf("Rank %d received bottom right corner: ", rank);
    // for (int i = 0; i < 4; i++) {
    //     printf("%d ", recv_bottom_right[i]);
    // }
    // printf("\n");

    // Copy the current matrix into the center of the extended matrix
    for (int i = 0; i < n_loc_r; i++) {
        for (int j = 0; j < n_loc_c; j++) {
            extended_matrix[i + 2][j + 2] = current[i][j];
        }
    }

    // Update the extended matrix with received data
    for (int j = 0; j < n_loc_c; j++) {
        extended_matrix[0][j + 2] = recv_top[j];
        extended_matrix[1][j + 2] = recv_top[j + n_loc_c];
        extended_matrix[n_loc_r + 2][j + 2] = recv_bottom[j];
        extended_matrix[n_loc_r + 3][j + 2] = recv_bottom[j + n_loc_c];
    }

    for (int i = 0; i < n_loc_r; i++) {
        extended_matrix[i + 2][0] = recv_left[i + n_loc_r];
        extended_matrix[i + 2][1] = recv_left[i];
        extended_matrix[i + 2][n_loc_c + 2] = recv_right[i];
        extended_matrix[i + 2][n_loc_c + 3] = recv_right[i + n_loc_r];
    }

    // Update the extended matrix with received corner data
    extended_matrix[0][0] = recv_top_left[3];
    extended_matrix[0][1] = recv_top_left[1];
    extended_matrix[1][0] = recv_top_left[2];
    extended_matrix[1][1] = recv_top_left[0];

    extended_matrix[0][n_loc_c + 2] = recv_top_right[1];
    extended_matrix[0][n_loc_c + 3] = recv_top_right[3];
    extended_matrix[1][n_loc_c + 2] = recv_top_right[2];
    extended_matrix[1][n_loc_c + 3] = recv_top_right[0];

    extended_matrix[n_loc_r + 2][0] = recv_bottom_left[1];
    extended_matrix[n_loc_r + 2][1] = recv_bottom_left[3];
    extended_matrix[n_loc_r + 3][0] = recv_bottom_left[0];
    extended_matrix[n_loc_r + 3][1] = recv_bottom_left[2];

    extended_matrix[n_loc_r + 2][n_loc_c + 2] = recv_bottom_right[1];
    extended_matrix[n_loc_r + 2][n_loc_c + 3] = recv_bottom_right[3];
    extended_matrix[n_loc_r + 3][n_loc_c + 2] = recv_bottom_right[2];
    extended_matrix[n_loc_r + 3][n_loc_c + 3] = recv_bottom_right[0];

    free(recv_top);
    free(recv_bottom);
    free(send_left);
    free(recv_right);
    free(send_right);
    free(recv_left);
    free(send_top_left);
    free(recv_bottom_right);
    free(send_top_right);
    free(recv_bottom_left);
    free(send_bottom_left);
    free(recv_top_right);
    free(send_bottom_right);
    free(recv_top_left);
}
