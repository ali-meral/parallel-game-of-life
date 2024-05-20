#include "utilities.h"

void count_cells(int n_loc_r, int n_loc_c, uint8_t matrix[n_loc_r][n_loc_c], int *alive_count, int *dead_count)
{
  *alive_count = 0;
  *dead_count = 0;
  for (int i = 0; i < n_loc_r; i++)
  {
    for (int j = 0; j < n_loc_c; j++)
    {
      if (matrix[i][j] == 1)
      {
        (*alive_count)++;
      }
      else
      {
        (*dead_count)++;
      }
    }
  }
}

int compare_matrices(int n, uint8_t (*matrix1)[n], uint8_t (*matrix2)[n]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (matrix1[i][j] != matrix2[i][j]) {
                return 0; // Matrices are not equal
            }
        }
    }
    return 1; // Matrices are equal
}