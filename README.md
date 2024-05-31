## HPC Assignment 2

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

#### Sequential Run

To run the sequential program, use the following command:

```sh
./main_sequential -n <grid_size> -s <seed> -d <density> -i <iterations> [-v]
```

Example:

```sh
./main_sequential -n 8 -s 2 -d 20 -i 10
```

#### Parallel Run

To run the parallel program, use the following command:

```sh
mpirun -np <number_of_processes> ./main_{sendrecv,collectives} -n <grid_size> -s <seed> -d <density> -i <iterations> [-v] [-c]
```

Example:

```sh
mpirun -np 4 ./main_sendrecv -n 8 -s 2 -d 20 -i 10
```

- `-np <number_of_processes>`: Number of processes
- `-n <grid_size>`: Size of the grid
- `-s <seed>`: Seed for random number generation
- `-d <density>`: Density percentage
- `-i <iterations>`: Number of iterations
- `-v`: Verbose mode (optional)
- `-c`: Verify the results by comparing with the sequential computation (optional)


#### Repeated Runs

To run the sequential program multiple times, use the `make run_sequential` target:

```sh
make run_sequential N=<grid_size> SEED=<seed> DENSITY=<density> ITERATIONS=<iterations> REPS=<repetitions>
```

Example:

```sh
make run_sequential N=1000 SEED=2 DENSITY=20 ITERATIONS=10 REPS=3
```

To run the parallel program multiple times use:

```sh
make {run_sendrecv,run_collectives} P=<number_of_processes> N=<grid_size> SEED=<seed> DENSITY=<density> ITERATIONS=<iterations> REPS=<repetitions>
```

Example:

```sh
make run_sendrecv P=4 N=1000 SEED=2 DENSITY=20 ITERATIONS=10 REPS=3
```

- `P=<number_of_processes>`: Number of processes (default: 4, only for parallel)
- `N=<grid_size>`: Size of the grid (default: 100)
- `SEED=<seed>`: Seed for random number generation (default: 1)
- `DENSITY=<density>`: Density percentage (default: 30)
- `ITERATIONS=<iterations>`: Number of iterations (default: 10)
- `REPS=<repetitions>`: Number of repetitions (default: 10)

### Output

#### Sequential

```sh
make run_sequential N=1000 SEED=1 DENSITY=30 ITERATIONS=10 REPS=3
```

```raw
reps    n       seed    density iters   implementation  time (ms)       alive   dead
1       1000    1       30      10      sequential      152.763877      365156  634844
2       1000    1       30      10      sequential      149.611482      365156  634844
3       1000    1       30      10      sequential      149.281325      365156  634844
```

#### Parallel SendRecv

```sh
make run_sendrecv P=4 N=1000 SEED=1 DENSITY=30 ITERATIONS=10 REPS=3
```
```raw
make run_sendrecv P=4 N=1000 SEED=1 DENSITY=30 ITERATIONS=10 REPS=3
reps    np      n       seed    density iters   dimx    dimy    implementation  time (ms)       alive   dead
1       4       1000    1       30      10      2       2       sendrecv        34.868822       365156  634844
2       4       1000    1       30      10      2       2       sendrecv        44.094889       365156  634844
3       4       1000    1       30      10      2       2       sendrecv        35.395157       365156  634844
```
