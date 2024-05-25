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

#### Single Run

To run the program once, use the following command:

```sh
mpirun -np <number_of_processes> ./main -n <grid_size> -s <seed> -d <density> -i <iterations> [-v] [-c]
```

Example:

```sh
mpirun -np 1 ./main -n 10 -s 1 -d 30 -v -i 10 -c
```
```sh
mpirun -np 2 ./main -n 10 -s 1 -d 30 -v -i 10 -c
```

- `-np <number_of_processes>`: Number of processes (e.g., 1)
- `-n <grid_size>`: Size of the grid (e.g., 10)
- `-s <seed>`: Seed for random number generation (e.g., 1)
- `-d <density>`: Density percentage (e.g., 30)
- `-i <iterations>`: Number of iterations (e.g., 10)
- `-v`: Verbose mode (optional)
- `-c`: Verify the results by comparing with the sequential computation (optional)


#### Repeated Runs

To run the program multiple times with customizable parameters, use the `make run` target:

```sh
make run P=<number_of_processes> N=<grid_size> SEED=<seed> DENSITY=<density> ITERATIONS=<iterations> REPS=<repetitions>
```

Example:

```sh
make run P=4 N=100 SEED=1 DENSITY=30 ITERATIONS=2 REPS=10
```

- `P=<number_of_processes>`: Number of processes (default: 1)
- `N=<grid_size>`: Size of the grid (default: 100)
- `SEED=<seed>`: Seed for random number generation (default: 2)
- `DENSITY=<density>`: Density percentage (default: 20)
- `ITERATIONS=<iterations>`: Number of iterations (default: 5)
- `REPS=<repetitions>`: Number of repetitions (default: 1)

### Example Usage

```sh
make clean
make
mpirun -np 1 ./main -n 10 -s 1 -d 30 -v -i 10 -c
make run P=4 N=100 SEED=1 DENSITY=30 ITERATIONS=2 REPS=10
```
