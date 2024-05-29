# Default parameters, these can be overridden at runtime
P ?= 4
N ?= 8
SEED ?= 1
DENSITY ?= 30
ITERATIONS ?= 1
REPS ?= 1
VERIFY ?= 1
VERBOSE ?= 0

# Compile the program
all: main_parallel main_sequential main_collectives

# Parallel target
main_parallel: main_parallel.o game.o  matrix_operations.o utilities.o mpi_communication.o
	mpicc -Wall -I./include -o $@ $^

# Collective target
main_collectives: main_collectives.o game.o matrix_operations.o utilities.o mpi_communication.o
	mpicc -Wall -I./include -o $@ $^

# Sequential target
main_sequential: main_sequential.o game.o matrix_operations.o utilities.o mpi_communication.o
	mpicc -Wall -I./include -o $@ $^

# Compile main_parallel.c
main_parallel.o: src/main_parallel.c
	mpicc -Wall -I./include -c $< -o $@

# Compile main_collectives.c
main_collectives.o: src/main_collectives.c
	mpicc -Wall -I./include -c $< -o $@

# Compile main_sequential.c
main_sequential.o: src/main_sequential.c
	mpicc -Wall -I./include -c $< -o $@

# Compile game.c
game.o: src/game.c
	mpicc -Wall -I./include -c $< -o $@

# Compile run_sequential.c
run_sequential.o: src/run_sequential.c
	mpicc -Wall -I./include -c $< -o $@

# Compile matrix_operations.c
matrix_operations.o: src/matrix_operations.c
	mpicc -Wall -I./include -c $< -o $@

# Compile utilities.c
utilities.o: src/utilities.c
	mpicc -Wall -I./include -c $< -o $@

# Compile mpi_communication.c
mpi_communication.o: src/mpi_communication.c
	mpicc -Wall -I./include -c $< -o $@

# Run the parallel program with custom parameters
run_parallel: main_parallel
	@echo "Running parallel simulation $(REPS) times with the following settings:"
	@echo "P: $(P) , Grid size: $(N), Seed: $(SEED), Density: $(DENSITY)%, Iterations: $(ITERATIONS), Verbose: $(VERBOSE)"
	@for i in $$(seq 1 $(REPS)); do \
		echo "Repetition $$i:"; \
		mpirun -np $(P) ./main_parallel -n $(N) -s $(SEED) -d $(DENSITY) -i $(ITERATIONS); \
	done

# Run the collectives program
run_collectives: main_collectives
	@echo "Running with collectives $(REPS) times:"
	@echo "P: $(P) , Grid size: $(N), Seed: $(SEED), Density: $(DENSITY)%, Iterations: $(ITERATIONS), Verbose: $(VERBOSE)"
	@for i in $$(seq 1 $(REPS)); do \
		echo "Repetition $$i:"; \
		mpirun -np $(P) ./main_collectives -n $(N) -s $(SEED) -d $(DENSITY) -i $(ITERATIONS); \
	done

# Run the sequential program with custom parameters
run_sequential: main_sequential
	@echo "Running sequential simulation $(REPS) times with the following settings:"
	@echo "Grid size: $(N), Seed: $(SEED), Density: $(DENSITY)%, Iterations: $(ITERATIONS), Verbose: $(VERBOSE)"
	@for i in $$(seq 1 $(REPS)); do \
		echo "Repetition $$i:"; \
		./main_sequential -n $(N) -s $(SEED) -d $(DENSITY) -i $(ITERATIONS) -v $(VERBOSE); \
	done

# Clean up binary files and objects
clean:
	rm -f main_parallel main_collectives main_sequential s*.o