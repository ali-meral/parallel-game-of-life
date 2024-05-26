# Default parameters, these can be overridden at runtime
P ?= 1
N ?= 100
SEED ?= 2
DENSITY ?= 20
ITERATIONS ?= 5
REPS ?= 1

# Compile the program
all: main

# Link objects and build the main executable
main: main.o run_simulation.o matrix_operations.o utilities.o mpi_communication.o
	mpicc -Wall -I./include -o $@ $^

# Compile main.c
main.o: main.c
	mpicc -Wall -I./include -c $< -o $@

# Compile run_simulation.c
run_simulation.o: src/run_simulation.c
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



# Run the program with custom parameters
run:
	@echo "Running simulation $(REPS) times with the following settings:"
	@echo "C: $(P) ,Grid size: $(N), Seed: $(SEED), Density: $(DENSITY)%, Iterations: $(ITERATIONS)"
	@for i in $$(seq 1 $(REPS)); do \
		echo "Repetition $$i:"; \
		mpirun -np $(P) ./main -n $(N) -s $(SEED) -d $(DENSITY) -i $(ITERATIONS); \
	done

# Clean up binary files and objects
clean:
	rm -f main *.o

.PHONY: all clean run
