CC=mpicc
CFLAGS=-Wall
EXECUTABLE=main
SRC=main.c

# Default values if not provided
REPS ?= 1
N ?= 100
SEED ?= 2
DENSITY ?= 20
VERBOSE ?= v
ITERATIONS ?= 5

all: $(EXECUTABLE)

$(EXECUTABLE): $(SRC)
	$(CC) $(CFLAGS) -o $@ $^

run:
	@echo "Running $(EXECUTABLE) $(REPS) times with the following parameters:"
	@echo "Grid size: $(N), Seed: $(SEED), Density: $(DENSITY)%, Verbose: $(VERBOSE), Iterations: $(ITERATIONS)"
	@for i in $$(seq 1 $(REPS)); do \
		echo "Repetition $$i:"; \
		mpirun -np 1 ./$(EXECUTABLE) -n $(N) -s $(SEED) -d $(DENSITY) -$(VERBOSE) -i $(ITERATIONS); \
	done | tee times.txt

clean:
	rm -f $(EXECUTABLE)

.PHONY: all run clean
