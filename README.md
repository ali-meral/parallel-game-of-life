### Overview

HPC Assignment 2

### Compilation

To compile the program, run:

```sh
make
```

To clean up the compiled binaries and object files, run:

```sh
make clean
```

### Running the Program

#### Parallel Run

To run the parallel program, use the following command:

```sh
mpirun -np <number_of_processes> ./bin/main_parallel -n <grid_size> -s <seed> -d <density> -i <iterations> [-v] [-c]
```

Example:

```sh
mpirun -np 4 ./bin/main_parallel -n 1000 -s 2 -d 20 -i 10 -v -c
```

- `-np <number_of_processes>`: Number of processes (e.g., 4)
- `-n <grid_size>`: Size of the grid (e.g., 1000)
- `-s <seed>`: Seed for random number generation (e.g., 2)
- `-d <density>`: Density percentage (e.g., 20)
- `-i <iterations>`: Number of iterations (e.g., 10)
- `-v`: Verbose mode (optional)
- `-c`: Verify the results by comparing with the sequential computation (optional)

#### Sequential Run

To run the sequential program, use the following command:

```sh
./bin/main_sequential -n <grid_size> -s <seed> -d <density> -i <iterations> [-v] [-c]
```

Example:

```sh
./bin/main_sequential -n 1000 -s 2 -d 20 -i 10 -v -c
```

- `-n <grid_size>`: Size of the grid (e.g., 1000)
- `-s <seed>`: Seed for random number generation (e.g., 2)
- `-d <density>`: Density percentage (e.g., 20)
- `-i <iterations>`: Number of iterations (e.g., 10)
- `-v`: Verbose mode (optional)
- `-c`: Verify the results by comparing with the parallel computation (optional)

#### Repeated Runs

To run the parallel program multiple times with customizable parameters, use the `make run_parallel` target:

```sh
make run_parallel P=<number_of_processes> N=<grid_size> SEED=<seed> DENSITY=<density> ITERATIONS=<iterations> REPS=<repetitions>
```

Example:

```sh
make run_parallel P=4 N=1000 SEED=2 DENSITY=20 ITERATIONS=10 REPS=3
```

To run the sequential program multiple times with customizable parameters, use the `make run_sequential` target:

```sh
make run_sequential N=<grid_size> SEED=<seed> DENSITY=<density> ITERATIONS=<iterations> REPS=<repetitions>
```

Example:

```sh
make run_sequential N=1000 SEED=2 DENSITY=20 ITERATIONS=10 REPS=3
```

- `P=<number_of_processes>`: Number of processes (default: 1, only for parallel)
- `N=<grid_size>`: Size of the grid (default: 1000)
- `SEED=<seed>`: Seed for random number generation (default: 2)
- `DENSITY=<density>`: Density percentage (default: 20)
- `ITERATIONS=<iterations>`: Number of iterations (default: 10)
- `REPS=<repetitions>`: Number of repetitions (default: 1)

### Example Usage

```sh
make clean
make
mpirun -np 4 ./main_parallel -n 8 -s 1 -d 30 -i 10 -v -c
./main_sequential -n 8 -s 1 -d 30 -i 10 -v
make run_parallel P=4 N=1000 SEED=1 DENSITY=30 ITERATIONS=10 REPS=3
make run_sequential N=1000 SEED=1 DENSITY=30 ITERATIONS=10 REPS=3
```

This `README.md` file provides an overview of the compilation and execution process for both the parallel and sequential versions of your program. The commands and examples are tailored to reflect the current Makefile and the directory structure where the binaries and object files are stored in the `bin` directory.