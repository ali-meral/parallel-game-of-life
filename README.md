# hpc
HPC Project 2024SS

### 1. Sequential


make clean
make

To run with repetitions:
make run N=100 SEED=1 DENSITY=30 ITERATIONS=2 REPS=10

To run once:
mpirun -np 1 ./main -n 10 -s 1 -d 30 -v -i 10