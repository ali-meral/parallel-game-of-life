#include "matrix_operations.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

void collectives_communicate(int n_loc_r, int n_loc_c, uint8_t (*current)[n_loc_c], uint8_t (*extended_matrix)[n_loc_c + 4], MPI_Comm cartcomm)
{
    int extended_r = n_loc_r + 4, extended_c = n_loc_c + 4;

    // Initialize the extended matrix with zeros
    for (int i = 0; i < extended_r; i++)
    {
        for (int j = 0; j < extended_c; j++)
        {
            extended_matrix[i][j] = 0;
        }
    }

    // Copy the current matrix to the center of the extended matrix
    for (int i = 2; i < n_loc_r + 2; i++)
    {
        for (int j = 2; j < n_loc_c + 2; j++)
        {
            extended_matrix[i][j] = current[i - 2][j - 2];
        }
    }

    int coords[2], dims[2], periods[2], reorder;
    MPI_Cart_get(cartcomm, 2, dims, periods, coords);

    int neighbors[8];
    int upper_left_coords[2], upper_right_coords[2], lower_left_coords[2], lower_right_coords[2];

    upper_left_coords[0] = (coords[0] - 1 + dims[0]) % dims[0];
    upper_left_coords[1] = (coords[1] - 1 + dims[1]) % dims[1];

    upper_right_coords[0] = (coords[0] - 1 + dims[0]) % dims[0];
    upper_right_coords[1] = (coords[1] + 1) % dims[1];

    lower_left_coords[0] = (coords[0] + 1) % dims[0];
    lower_left_coords[1] = (coords[1] - 1 + dims[1]) % dims[1];

    lower_right_coords[0] = (coords[0] + 1) % dims[0];
    lower_right_coords[1] = (coords[1] + 1) % dims[1];

    MPI_Cart_rank(cartcomm, upper_left_coords, &neighbors[0]);
    MPI_Cart_rank(cartcomm, upper_right_coords, &neighbors[1]);
    MPI_Cart_rank(cartcomm, lower_left_coords, &neighbors[2]);
    MPI_Cart_rank(cartcomm, lower_right_coords, &neighbors[3]);
    MPI_Cart_shift(cartcomm, 1, 1, &neighbors[4], &neighbors[5]);
    MPI_Cart_shift(cartcomm, 0, 1, &neighbors[6], &neighbors[7]);

    MPI_Comm dist_graph_comm;
    MPI_Dist_graph_create_adjacent(cartcomm, 8, neighbors, MPI_UNWEIGHTED, 8, neighbors, MPI_UNWEIGHTED, MPI_INFO_NULL, 0, &dist_graph_comm);

    int total_send_size = 4 * 4 + 2 * 2 * n_loc_r + 2 * 2 * n_loc_c;
    int total_recv_size = total_send_size;

    uint8_t *sendbuf = (uint8_t *)malloc(total_send_size * sizeof(uint8_t));
    uint8_t *recvbuf = (uint8_t *)malloc(total_recv_size * sizeof(uint8_t));

    int sendcounts[8], recvcounts[8];
    int senddispls[8], recvdispls[8];

    int pos = 0;
    senddispls[0] = pos;
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            sendbuf[pos++] = current[i][j];
        }
    }
    senddispls[1] = pos;
    for (int i = 0; i < 2; i++)
    {
        for (int j = n_loc_c - 2; j < n_loc_c; j++)
        {
            sendbuf[pos++] = current[i][j];
        }
    }
    senddispls[2] = pos;
    for (int i = n_loc_r - 2; i < n_loc_r; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            sendbuf[pos++] = current[i][j];
        }
    }
    senddispls[3] = pos;
    for (int i = n_loc_r - 2; i < n_loc_r; i++)
    {
        for (int j = n_loc_c - 2; j < n_loc_c; j++)
        {
            sendbuf[pos++] = current[i][j];
        }
    }
    senddispls[4] = pos;
    for (int i = 0; i < n_loc_r; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            sendbuf[pos++] = current[i][j];
        }
    }
    senddispls[5] = pos;
    for (int i = 0; i < n_loc_r; i++)
    {
        for (int j = n_loc_c - 2; j < n_loc_c; j++)
        {
            sendbuf[pos++] = current[i][j];
        }
    }
    senddispls[6] = pos;
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < n_loc_c; j++)
        {
            sendbuf[pos++] = current[i][j];
        }
    }
    senddispls[7] = pos;
    for (int i = n_loc_r - 2; i < n_loc_r; i++)
    {
        for (int j = 0; j < n_loc_c; j++)
        {
            sendbuf[pos++] = current[i][j];
        }
    }

    for (int i = 0; i < 8; i++)
    {
        sendcounts[i] = recvcounts[i] = (i < 4) ? 4 : ((i < 6) ? n_loc_r * 2 : n_loc_c * 2);
        recvdispls[i] = senddispls[i];
    }

    MPI_Neighbor_alltoallv(sendbuf, sendcounts, senddispls, MPI_UINT8_T,
                           recvbuf, recvcounts, recvdispls, MPI_UINT8_T,
                           dist_graph_comm);

    // Unpack the received data into the extended matrix
    pos = 0;

    // Fill the corners

    // Lower right
    for (int i = n_loc_r + 2; i < n_loc_r + 4; i++)
    {
        for (int j = n_loc_c + 2; j < n_loc_c + 4; j++)
        {
            extended_matrix[i][j] = recvbuf[pos++];
        }
    }
    // Lower left
    for (int i = n_loc_r + 2; i < n_loc_r + 4; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            extended_matrix[i][j] = recvbuf[pos++];
        }
    }
    // Upper right
    for (int i = 0; i < 2; i++)
    {
        for (int j = n_loc_c + 2; j < n_loc_c + 4; j++)
        {
            extended_matrix[i][j] = recvbuf[pos++];
        }
    }
    // Upper left
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            extended_matrix[i][j] = recvbuf[pos++];
        }
    }

    // Fill the edges
    // Lower edge
    for (int i = 2; i < n_loc_r + 2; i++)
    {
        for (int j = n_loc_c + 2; j < n_loc_c + 4; j++)
        {
            extended_matrix[i][j] = recvbuf[pos++];
        }
    }
    // Upper edge
    for (int i = 2; i < n_loc_r + 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            extended_matrix[i][j] = recvbuf[pos++];
        }
    }

    // Right edge
    for (int i = n_loc_r + 2; i < n_loc_r + 4; i++)
    {
        for (int j = 2; j < n_loc_c + 2; j++)
        {
            extended_matrix[i][j] = recvbuf[pos++];
        }
    }

    // Left edge
    for (int i = 0; i < 2; i++)
    {
        for (int j = 2; j < n_loc_c + 2; j++)
        {
            extended_matrix[i][j] = recvbuf[pos++];
        }
    }

    free(sendbuf);
    free(recvbuf);
}

void sendrecv_communicate(int n_loc_r, int n_loc_c, uint8_t (*current)[n_loc_c], uint8_t (*extended_matrix)[n_loc_c + 4], MPI_Comm cartcomm)
{
    MPI_Status status;
    int rank, coords[2], dims[2], periods[2];
    MPI_Comm_rank(cartcomm, &rank);
    MPI_Cart_get(cartcomm, 2, dims, periods, coords);

    int neighbor_coords[8][2];
    int neighbors[8];

    // relative positions of the 8 neighbors
    int rel_pos[8][2] = {
        {-1, -1}, {-1, 0}, {-1, 1}, {0, -1}, {0, 1}, {1, -1}, {1, 0}, {1, 1}};

    // get coordinates and ranks of the 8 neighbors
    for (int i = 0; i < 8; i++)
    {
        neighbor_coords[i][0] = (coords[0] + rel_pos[i][0] + dims[0]) % dims[0];
        neighbor_coords[i][1] = (coords[1] + rel_pos[i][1] + dims[1]) % dims[1];
        MPI_Cart_rank(cartcomm, neighbor_coords[i], &neighbors[i]);
    }

    // reserve space for the buffers
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

    // buffers for top and bottom rows
    for (int j = 0; j < n_loc_c; j++)
    {
        send_top[j] = current[0][j];
        send_top[j + n_loc_c] = current[1][j];
        send_bottom[j] = current[n_loc_r - 2][j];
        send_bottom[j + n_loc_c] = current[n_loc_r - 1][j];
    }

    // buffers for left and right columns
    for (int i = 0; i < n_loc_r; i++)
    {
        send_left[i] = current[i][0];
        send_left[i + n_loc_r] = current[i][1];
        send_right[i] = current[i][n_loc_c - 1];
        send_right[i + n_loc_r] = current[i][n_loc_c - 2];
    }

    // buffers for corners
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

    // top and bottom rows
    MPI_Sendrecv(send_top, 2 * n_loc_c, MPI_UINT8_T, neighbors[1], 0,
                 recv_bottom, 2 * n_loc_c, MPI_UINT8_T, neighbors[6], 0,
                 cartcomm, &status);

    MPI_Sendrecv(send_bottom, 2 * n_loc_c, MPI_UINT8_T, neighbors[6], 1,
                 recv_top, 2 * n_loc_c, MPI_UINT8_T, neighbors[1], 1,
                 cartcomm, &status);

    // left and right columns
    MPI_Sendrecv(send_left, 2 * n_loc_r, MPI_UINT8_T, neighbors[3], 2,
                 recv_right, 2 * n_loc_r, MPI_UINT8_T, neighbors[4], 2,
                 cartcomm, &status);

    MPI_Sendrecv(send_right, 2 * n_loc_r, MPI_UINT8_T, neighbors[4], 3,
                 recv_left, 2 * n_loc_r, MPI_UINT8_T, neighbors[3], 3,
                 cartcomm, &status);

    // corners
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

    // place the original matrix into the extended matrix
    for (int i = 0; i < n_loc_r; i++)
    {
        for (int j = 0; j < n_loc_c; j++)
        {
            extended_matrix[i + 2][j + 2] = current[i][j];
        }
    }

    // copy received into extended matrix
    for (int j = 0; j < n_loc_c; j++)
    {
        extended_matrix[0][j + 2] = recv_top[j];
        extended_matrix[1][j + 2] = recv_top[j + n_loc_c];
        extended_matrix[n_loc_r + 2][j + 2] = recv_bottom[j];
        extended_matrix[n_loc_r + 3][j + 2] = recv_bottom[j + n_loc_c];
    }

    for (int i = 0; i < n_loc_r; i++)
    {
        extended_matrix[i + 2][0] = recv_left[i + n_loc_r];
        extended_matrix[i + 2][1] = recv_left[i];
        extended_matrix[i + 2][n_loc_c + 2] = recv_right[i];
        extended_matrix[i + 2][n_loc_c + 3] = recv_right[i + n_loc_r];
    }

    extended_matrix[0][0] = recv_top_left[3]; // farthest point
    extended_matrix[0][1] = recv_top_left[1]; // next layer
    extended_matrix[1][0] = recv_top_left[2]; // next layer
    extended_matrix[1][1] = recv_top_left[0]; // closest point

    extended_matrix[0][n_loc_c + 2] = recv_top_right[1];
    extended_matrix[0][n_loc_c + 3] = recv_top_right[3];
    extended_matrix[1][n_loc_c + 2] = recv_top_right[0];
    extended_matrix[1][n_loc_c + 3] = recv_top_right[2];

    extended_matrix[n_loc_r + 2][0] = recv_bottom_left[1];
    extended_matrix[n_loc_r + 2][1] = recv_bottom_left[0];
    extended_matrix[n_loc_r + 3][0] = recv_bottom_left[3];
    extended_matrix[n_loc_r + 3][1] = recv_bottom_left[2];

    extended_matrix[n_loc_r + 2][n_loc_c + 2] = recv_bottom_right[0];
    extended_matrix[n_loc_r + 2][n_loc_c + 3] = recv_bottom_right[1];
    extended_matrix[n_loc_r + 3][n_loc_c + 2] = recv_bottom_right[2];
    extended_matrix[n_loc_r + 3][n_loc_c + 3] = recv_bottom_right[3];

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

void gather_ranks(int n_loc_r, int n_loc_c, int n, int rank, int size, int m_offset_r, int m_offset_c,
                  uint8_t (*current)[n_loc_c], uint8_t (*global_matrix_par)[n], MPI_Comm cartcomm,
                  int verbose, int iterations)
{
    if (rank == 0)
    {
        for (int i = 0; i < n_loc_r; ++i)
        {
            for (int j = 0; j < n_loc_c; ++j)
            {
                global_matrix_par[m_offset_r + i][m_offset_c + j] = current[i][j];
            }
        }

        for (int r = 1; r < size; ++r)
        {
            uint8_t *recv_buffer = (uint8_t *)malloc(n_loc_r * n_loc_c * sizeof(uint8_t));
            MPI_Recv(recv_buffer, n_loc_r * n_loc_c, MPI_UINT8_T, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            int coords[2];
            MPI_Cart_coords(cartcomm, r, 2, coords);
            int recv_offset_r = coords[0] * n_loc_r;
            int recv_offset_c = coords[1] * n_loc_c;

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
        MPI_Send(&(current[0][0]), n_loc_r * n_loc_c, MPI_UINT8_T, 0, 0, MPI_COMM_WORLD);
    }
}
