#include "matrix_operations.h"
#include <stdio.h>
#include <stdlib.h>

void fill_matrix(int n_loc_r, int n_loc_c, uint8_t (*matrix)[n_loc_c], int n, int density, int m_offset_r, int m_offset_c)
{
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      int r = rand() % 100;
      if (r < density)
      {
        r = 1;
      }
      else
      {
        r = 0;
      }

      if (i >= m_offset_r && i < m_offset_r + n_loc_r &&
          j >= m_offset_c && j < m_offset_c + n_loc_c)
      {
        matrix[i - m_offset_r][j - m_offset_c] = r;
      }
    }
  }
}

void print_matrix(int n_loc_r, int n_loc_c, uint8_t (*matrix)[n_loc_c], int rank, int size)
{
  for (int i = 0; i < size; i++)
  {
    if (rank == i)
    {
      for (int i = 0; i < n_loc_r; i++)
      {
        for (int j = 0; j < n_loc_c; j++)
        {
          printf("%d ", matrix[i][j]);
        }
        printf("\n");
      }
    }
    fflush(stdout);
  }
}

void update_matrix(int n_loc_r, int n_loc_c, uint8_t (*current)[n_loc_c], uint8_t (*next)[n_loc_c])
{
  int neighbors;
  // Use extended indices for seamless wrap-around with an extra row/column at each boundary
  int extended_r = n_loc_r + 2;
  int extended_c = n_loc_c + 2;
  uint8_t extended[extended_r][extended_c];

  // Prepare extended matrix with wrap-around
  for (int i = 0; i < n_loc_r; i++)
  {
    for (int j = 0; j < n_loc_c; j++)
    {
      extended[i + 1][j + 1] = current[i][j];
    }
  }
  // Fill in wrap-around cells
  for (int i = 1; i <= n_loc_r; i++)
  {
    extended[i][0] = current[i - 1][n_loc_c - 1];
    extended[i][n_loc_c + 1] = current[i - 1][0];
  }
  for (int j = 1; j <= n_loc_c; j++)
  {
    extended[0][j] = current[n_loc_r - 1][j - 1];
    extended[n_loc_r + 1][j] = current[0][j - 1];
  }
  // Corners
  extended[0][0] = current[n_loc_r - 1][n_loc_c - 1];
  extended[0][n_loc_c + 1] = current[n_loc_r - 1][0];
  extended[n_loc_r + 1][0] = current[0][n_loc_c - 1];
  extended[n_loc_r + 1][n_loc_c + 1] = current[0][0];

  // Update cells without conditional statements
  for (int i = 1; i <= n_loc_r; i++)
  {
    for (int j = 1; j <= n_loc_c; j++)
    {
      neighbors = extended[i - 1][j - 1] + extended[i - 1][j] + extended[i - 1][j + 1] +
                  extended[i][j - 1] + extended[i][j + 1] +
                  extended[i + 1][j - 1] + extended[i + 1][j] + extended[i + 1][j + 1];
      // Compute next state using arithmetic
      next[i - 1][j - 1] = neighbors == 3 || (neighbors == 2 && extended[i][j]);
    }
  }
}

void update_matrix_debug(int n_loc_r, int n_loc_c, uint8_t (*current)[n_loc_c], uint8_t (*next)[n_loc_c], int rank, int size, int generation)
{
  // Call the standard update function first (assuming it's already implemented)
  update_matrix(n_loc_r, n_loc_c, current, next);

  // Print the 'next' matrix with generation info
  if (rank == 0)
  { // Assuming rank 0 is responsible for output
    printf("Generation %d:\n", generation);
    print_matrix(n_loc_r, n_loc_c, current, rank, size);
  }
}