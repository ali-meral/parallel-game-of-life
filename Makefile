# Default parameters, these can be overridden at runtime
C ?= 1
N ?= 100
SEED ?= 2
DENSITY ?= 20
ITERATIONS ?= 5
REPS ?= 1

# Compile the program
all: main

# Link objects and build the main executable
main: main.o simulation_control.o matrix_operations.o utilities.o
	mpicc -Wall -I./include -o $@ $^

# Compile main.c
main.o: main.c
	mpicc -Wall -I./include -c $< -o $@

# Compile simulation_control.c
simulation_control.o: src/simulation_control.c
	mpicc -Wall -I./include -c $< -o $@

# Compile matrix_operations.c
matrix_operations.o: src/matrix_operations.c
	mpicc -Wall -I./include -c $< -o $@

# Compile utilities.c
utilities.o: src/utilities.c
	mpicc -Wall -I./include -c $< -o $@

# Run the program with custom parameters
run:
	@echo "Running simulation $(REPS) times with the following settings:"
	@echo "C: $(C) ,Grid size: $(N), Seed: $(SEED), Density: $(DENSITY)%, Iterations: $(ITERATIONS)"
	@for i in $$(seq 1 $(REPS)); do \
		echo "Repetition $$i:"; \
		mpirun -np $(C) ./main -n $(N) -s $(SEED) -d $(DENSITY) -i $(ITERATIONS); \
	done

# Clean up binary files and objects
clean:
	rm -f main *.o

.PHONY: all clean run
