#include "matrix_operations.h"
#include "utilities.h"
#include <stdio.h>
#include <stdlib.h>

int compare_matrices(int n, uint8_t (*matrix1)[n], uint8_t (*matrix2)[n])
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {

            if (matrix1[i][j] != matrix2[i][j])
            {
                printf("matrix1[%d][%d] = %d, matrix2[%d][%d] = %d\n", i, j, matrix1[i][j], i, j, matrix2[i][j]);
                return 0; // Matrices are not equal
            }
        }
    }
    return 1; // Matrices are equal
}

void fill_extended_grid(int n_loc_r, int n_loc_c, uint8_t (*current)[n_loc_c], uint8_t (*extended)[n_loc_c + 4]) {
    // Fill the center of the extended grid with the original grid
    for (int i = 0; i < n_loc_r; i++) {
        for (int j = 0; j < n_loc_c; j++) {
            extended[i + 2][j + 2] = current[i][j];
        }
    }

    // Fill the wrap-around rows
    for (int j = 0; j < n_loc_c; j++) {
        extended[0][j + 2] = current[n_loc_r - 2][j];
        extended[1][j + 2] = current[n_loc_r - 1][j];
        extended[n_loc_r + 2][j + 2] = current[0][j];
        extended[n_loc_r + 3][j + 2] = current[1][j];
    }

    // Fill the wrap-around columns
    for (int i = 0; i < n_loc_r; i++) {
        extended[i + 2][0] = current[i][n_loc_c - 2];
        extended[i + 2][1] = current[i][n_loc_c - 1];
        extended[i + 2][n_loc_c + 2] = current[i][0];
        extended[i + 2][n_loc_c + 3] = current[i][1];
    }

    // Fill the corners
    extended[0][0] = current[n_loc_r - 2][n_loc_c - 2];
    extended[0][1] = current[n_loc_r - 2][n_loc_c - 1];
    extended[0][n_loc_c + 2] = current[n_loc_r - 2][0];
    extended[0][n_loc_c + 3] = current[n_loc_r - 2][1];
    
    extended[1][0] = current[n_loc_r - 1][n_loc_c - 2];
    extended[1][1] = current[n_loc_r - 1][n_loc_c - 1];
    extended[1][n_loc_c + 2] = current[n_loc_r - 1][0];
    extended[1][n_loc_c + 3] = current[n_loc_r - 1][1];
    
    extended[n_loc_r + 2][0] = current[0][n_loc_c - 2];
    extended[n_loc_r + 2][1] = current[0][n_loc_c - 1];
    extended[n_loc_r + 2][n_loc_c + 2] = current[0][0];
    extended[n_loc_r + 2][n_loc_c + 3] = current[0][1];
    
    extended[n_loc_r + 3][0] = current[1][n_loc_c - 2];
    extended[n_loc_r + 3][1] = current[1][n_loc_c - 1];
    extended[n_loc_r + 3][n_loc_c + 2] = current[1][0];
    extended[n_loc_r + 3][n_loc_c + 3] = current[1][1];
}


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

void print_matrix_seq(int n_loc_r, int n_loc_c, uint8_t (*matrix)[n_loc_c]) {
    for (int i = 0; i < n_loc_r; i++) {
        for (int j = 0; j < n_loc_c; j++) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
    fflush(stdout);
}

void update_matrix(int n_loc_r, int n_loc_c, uint8_t (*extended_matrix)[n_loc_c + 4], uint8_t (*next)[n_loc_c])
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