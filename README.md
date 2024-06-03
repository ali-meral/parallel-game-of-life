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
mpirun -np <number_of_processes> ./main_{sendrecv,collectives} -n <grid_size> -s <seed> -d <density> -i <iterations> -r <repetitions> [-v] [-c] 
```

Example:

```sh
mpirun -np 4 ./main_sendrecv -n 8 -s 2 -d 20 -i 10 -r 1
```

- `-np <number_of_processes>`: Number of processes
- `-n <grid_size>`: Size of the grid
- `-s <seed>`: Seed for random number generation
- `-d <density>`: Density percentage
- `-i <iterations>`: Number of iterations
- `-v`: Verbose mode (optional)
- `-c`: Verify the results by comparing with the sequential computation (optional)
- `-r`: Repetitions


### Output

#### Sequential

```sh
./main_sequential -n 8 -s 2 -d 20 -i 10 -r 5
```

```raw
n       seed    density iters   implementation  time (ms)       alive   dead
8       2       20      10      sequential      0.007739        0       64
8       2       20      10      sequential      0.001500        0       64
8       2       20      10      sequential      0.001250        0       64
```

#### Parallel SendRecv

```sh
mpirun -np 4 ./main_sendrecv -n 8 -s 2 -d 20 -i 10 -r 3
```

```raw
np      n       seed    density iters   dimx    dimy    implementation  time (ms)       alive   dead
4       8       2       20      10      2       2       sendrecv        0.053482        0       64
4       8       2       20      10      2       2       sendrecv        0.028771        0       64
4       8       2       20      10      2       2       sendrecv        0.028081        0       64
```
