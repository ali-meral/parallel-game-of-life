#include <assert.h>
#include <getopt.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

void fill_matrix(int n_loc_r, int n_loc_c, uint8_t (*matrix)[n_loc_c], int n, int density, int m_offset_r, int m_offset_c);
void print_matrix(int n_loc_r, int n_loc_c, uint8_t (*matrix)[n_loc_c], int rank, int size);
void update_matrix(int n_loc_r, int n_loc_c, uint8_t (*current)[n_loc_c], uint8_t (*next)[n_loc_c]);

/**
 * @brief
 *
 * @param n_loc_r number of rows of local matrix
 * @param n_loc_c number of columns of local matrix
 * @param n       number of rows/columns of global matrix
 * @param density density of 1s in the matrix
 * @param iterations number of generations
 * @param m_offset_r matrix offset of rank i in rows
 * @param m_offset_c matrix offset of rank i in columns
 */
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
    MPI_Barrier(MPI_COMM_WORLD);
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
    for (int i = 0; i < n_loc_r; i++) {
        for (int j = 0; j < n_loc_c; j++) {
            extended[i + 1][j + 1] = current[i][j];
        }
    }
    // Fill in wrap-around cells
    for (int i = 1; i <= n_loc_r; i++) {
        extended[i][0] = current[i - 1][n_loc_c - 1];
        extended[i][n_loc_c + 1] = current[i - 1][0];
    }
    for (int j = 1; j <= n_loc_c; j++) {
        extended[0][j] = current[n_loc_r - 1][j - 1];
        extended[n_loc_r + 1][j] = current[0][j - 1];
    }
    // Corners
    extended[0][0] = current[n_loc_r - 1][n_loc_c - 1];
    extended[0][n_loc_c + 1] = current[n_loc_r - 1][0];
    extended[n_loc_r + 1][0] = current[0][n_loc_c - 1];
    extended[n_loc_r + 1][n_loc_c + 1] = current[0][0];

    // Update cells without conditional statements
    for (int i = 1; i <= n_loc_r; i++) {
        for (int j = 1; j <= n_loc_c; j++) {
            neighbors = extended[i - 1][j - 1] + extended[i - 1][j] + extended[i - 1][j + 1] +
                        extended[i][j - 1]                     + extended[i][j + 1] +
                        extended[i + 1][j - 1] + extended[i + 1][j] + extended[i + 1][j + 1];
            // Compute next state using arithmetic
            next[i - 1][j - 1] = neighbors == 3 || (neighbors == 2 && extended[i][j]);
        }
    }
}

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

int main(int argc, char *argv[])
{
  int rank, size; // always

  int seed = 10;
  int n = 10; // default values
  int opt;
  int verbose = 0;
  int density = 10;   // in percent
  int iterations = 3; // default number of iterations

  int n_loc_r, n_loc_c;
  int nprows, npcols;
  int prow_idx, pcol_idx;

  // Cartesian communicator stuff
  int dims[2] = {0, 0};
  int pers[2] = {0, 0};
  int coords[2];
  MPI_Comm cartcomm;

  static struct option long_options[] = {{"number", optional_argument, 0, 'n'},
                                         {"seed", optional_argument, 0, 's'},
                                         {"density", optional_argument, 0, 'd'},
                                         {"verbose", optional_argument, 0, 'v'},
                                         {"iterations", optional_argument, 0, 'i'},
                                         {0, 0, 0, 0}};

  MPI_Init(&argc, &argv);

  while (1)
  {
    int option_index = 0;
    opt = getopt_long(argc, argv, "n:s:d:vi:", long_options, &option_index);

    if (opt == -1)
      break;

    switch (opt)
    {
    case 'n':
      n = atoi(optarg);
      break;
    case 's':
      seed = atoi(optarg);
      break;
    case 'd':
      density = atoi(optarg);
      break;
    case 'i':
      iterations = atoi(optarg);
      break;
    case 'v':
      verbose = 1;
      break;
    default:
      fprintf(stderr, "Usage: %s -n <n> -k <k>\n", argv[0]);
      exit(EXIT_FAILURE);
    }
  }

  if (density <= 0 || density > 100)
  {
    fprintf(stderr, "Density should be between 1 and 100\n");
    exit(EXIT_FAILURE);
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  MPI_Dims_create(size, 2, dims);

  if (verbose == 1 && rank == 0)
  {
    printf("Dimensions created: [%d, %d]\n", dims[0], dims[1]);
  }

  int reorder = 0;
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, pers, reorder, &cartcomm);
  MPI_Cart_coords(cartcomm, rank, 2, coords);
  prow_idx = coords[0];
  pcol_idx = coords[1];

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

  srand(seed);

  // matrix with npcols; this declaration is enough to use matrix[i][j]
  uint8_t(*matrix)[n_loc_c];
  matrix = (uint8_t(*)[n_loc_c])malloc(n_loc_r * n_loc_c * sizeof(uint8_t));

  int m_offset_r = prow_idx * n_loc_r;
  int m_offset_c = pcol_idx * n_loc_c;
  if (verbose == 1)
  {
    printf("%d: prow_idx: %d pcol_idx: %d m_offset_r: %d m_offset_c: %d\n", rank, prow_idx, pcol_idx, m_offset_r, m_offset_c);
  }

  fill_matrix(n_loc_r, n_loc_c, matrix, n, density, m_offset_r, m_offset_c);

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  // Debugging options
  // printf("Generation 0:\n");
  // print_matrix(n_loc_r, n_loc_c, matrix, rank, size);


  free(matrix);

  // Allocate two matrices for current and next states
  uint8_t(*current)[n_loc_c] = (uint8_t(*)[n_loc_c])malloc(n_loc_r * n_loc_c * sizeof(uint8_t));
  uint8_t(*next)[n_loc_c] = (uint8_t(*)[n_loc_c])malloc(n_loc_r * n_loc_c * sizeof(uint8_t));

  fill_matrix(n_loc_r, n_loc_c, current, n, density, m_offset_r, m_offset_c);


  // Start timing
  double start_time, end_time, elapsed_time;
  MPI_Barrier(MPI_COMM_WORLD); // Synchronize before starting the timer
  start_time = MPI_Wtime();

  for (int gen = 0; gen < iterations; gen++)
  {
    update_matrix(n_loc_r, n_loc_c, current, next);
    // Swap pointers for the next iteration
    uint8_t(*temp)[n_loc_c] = current;
    current = next;
    next = temp;

    // Debugging options
    // printf("Generation %d:\n", gen + 1);
    // print_matrix(n_loc_r, n_loc_c, current, rank, size);
  }

  end_time = MPI_Wtime(); // Stop the timer before post-processing
  elapsed_time = end_time - start_time;

  int local_alive = 0, local_dead = 0;
  count_cells(n_loc_r, n_loc_c, current, &local_alive, &local_dead);

  int total_alive, total_dead;
  MPI_Reduce(&local_alive, &total_alive, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_dead, &total_dead, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0)
  {
    printf("Alive cells: %d, Dead cells: %d\n", total_alive, total_dead);

    // if (verbose == 1)
    // {
      printf("Computation time: %f ms\n\n", elapsed_time * 1000);
    // }
  }

  free(current);
  free(next);

  MPI_Finalize();

  return 0;
}

// how to compile and run sequentially
// mpicc -o main main.c && mpirun -np 1 ./main -n 10 -s 1 -d 20 -v -i 100
// how to compile and run sequentially and save the results for plotting
// mpicc -o main main.c && mpirun -np 1 ./main -n 10 -s 1 -d 20 -v -i 100  > input.txt
// how to run repetitions
// make clean
// make
// make run REPS=10 N=10 SEED=1 DENSITY=20 VERBOSE= ITERATIONS=100