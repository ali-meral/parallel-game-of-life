# input parameters
P ?= 4
N ?= 8
SEED ?= 1
DENSITY ?= 30
ITERATIONS ?= 2
REPS ?= 1
VERIFY ?= 1
VERBOSE ?= 0

# object files
OBJS = game.o matrix_operations.o utilities.o mpi_communication.o

# targets
all: main_parallel main_sequential main_collectives

# linkings
%: %.o $(OBJS)
	mpicc -I./include -o $@ $^

# Compilation pattern rule
%.o: src/%.c
	mpicc -I./include -c $< -o $@

# Rules
run_parallel: main_parallel
	@echo "SendRecv running $(REPS) times:"
	@echo "P: $(P) , Grid size: $(N), Seed: $(SEED), Density: $(DENSITY)%, Iterations: $(ITERATIONS), Verbose: $(VERBOSE)"
	@for i in $$(seq 1 $(REPS)); do \
		echo "Repetition $$i:"; \
		mpirun -np $(P) ./main_parallel -n $(N) -s $(SEED) -d $(DENSITY) -i $(ITERATIONS); \
	done

run_collectives: main_collectives
	@echo "Collectives running $(REPS) times:"
	@echo "P: $(P) , Grid size: $(N), Seed: $(SEED), Density: $(DENSITY)%, Iterations: $(ITERATIONS), Verbose: $(VERBOSE)"
	@for i in $$(seq 1 $(REPS)); do \
		echo "Repetition $$i:"; \
		mpirun -np $(P) ./main_collectives -n $(N) -s $(SEED) -d $(DENSITY) -i $(ITERATIONS); \
	done

run_sequential: main_sequential
	@echo "Sequential running $(REPS) times:"
	@echo "Grid size: $(N), Seed: $(SEED), Density: $(DENSITY)%, Iterations: $(ITERATIONS), Verbose: $(VERBOSE)"
	@for i in $$(seq 1 $(REPS)); do \
		echo "Repetition $$i:"; \
		./main_sequential -n $(N) -s $(SEED) -d $(DENSITY) -i $(ITERATIONS) -v $(VERBOSE); \
	done

clean:
	rm -f main_parallel main_collectives main_sequential *.o
