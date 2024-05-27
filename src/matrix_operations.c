#include "matrix_operations.h"
#include "utilities.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

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

int wrap(int idx, int limit)
{
  return idx + limit * (idx < 0) - limit * (idx >= limit);
}
void update_matrix(int n_loc_r, int n_loc_c, uint8_t (*current)[n_loc_c], uint8_t (*next)[n_loc_c])
{
  for (int i = 0; i < n_loc_r; i++)
  {
    for (int j = 0; j < n_loc_c; j++)
    {
      // precompute indices
      int im2 = wrap(i - 2, n_loc_r);
      int im1 = wrap(i - 1, n_loc_r);
      int ip2 = wrap(i + 2, n_loc_r);
      int jm2 = wrap(j - 2, n_loc_c);
      int jp1 = wrap(j + 1, n_loc_c);
      int jp2 = wrap(j + 2, n_loc_c);

      // calculate alive neighbors
      int neighbors =
          current[im2][jm2] + // C[i ⊖ 2][j ⊖ 2]
          current[im1][j] +   // C[i ⊖ 1][j]
          current[im1][jp1] + // C[i ⊖ 1][j ⊕ 1]
          current[i][jm2] +   // C[i][j ⊖ 2]
          current[i][jp1] +   // C[i][j ⊕ 1]
          current[ip2][jm2] + // C[i ⊕ 2][j ⊖ 2]
          current[ip2][j] +   // C[i ⊕ 2][j]
          current[ip2][jp2];  // C[i ⊕ 2][j ⊕ 2]

      // determine who lives and who dies
      next[i][j] = (neighbors == 3 || (neighbors == 2 && current[i][j]));
    }
  }
}

void print_matrix_debug(int n_loc_r, int n_loc_c, uint8_t (*matrix)[n_loc_c], int rank, const char *matrix_name)
{
  printf("Rank %d - %s:\n", rank, matrix_name);
  for (int i = 0; i < n_loc_r; i++)
  {
    for (int j = 0; j < n_loc_c; j++)
    {
      printf("%d ", matrix[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

void update_matrix_mpi(int n_loc_r, int n_loc_c, uint8_t (*extended_matrix)[n_loc_c + 4], uint8_t (*next)[n_loc_c])
{
  for (int i = 2; i < n_loc_r + 2; i++)
  {
    for (int j = 2; j < n_loc_c + 2; j++)
    {
      int neighbors =
          extended_matrix[i - 2][j - 2] +
          extended_matrix[i - 1][j] +
          extended_matrix[i - 1][j + 1] +
          extended_matrix[i][j - 2] +
          extended_matrix[i][j + 1] +
          extended_matrix[i + 2][j - 2] +
          extended_matrix[i + 2][j] +
          extended_matrix[i + 2][j + 2];

      next[i - 2][j - 2] = (neighbors == 3 || (neighbors == 2 && extended_matrix[i][j]));
    }
  }
}

