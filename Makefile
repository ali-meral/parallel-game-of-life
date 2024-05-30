# input parameters
P ?= 4
N ?= 100
SEED ?= 1
DENSITY ?= 30
ITERATIONS ?= 10
REPS ?= 10

# Compiler and flags
MPICC = mpicc
CFLAGS = -I./include -w  # Add -w to suppress all warnings

# object files
OBJS = game.o matrix_operations.o utilities.o mpi_communication.o

# targets
all: main_sendrecv main_sequential main_collectives

# linkings
%: %.o $(OBJS)
	$(MPICC) $(CFLAGS) -o $@ $^

# Compilation pattern rule
%.o: src/%.c
	$(MPICC) $(CFLAGS) -c $< -o $@

# Rules
run_sendrecv: main_sendrecv
	@echo "reps\tnp\tn\tseed\tdensity\titers\tdimx\tdimy\timplementation\ttime (ms)\talive\tdead"
	@for i in $$(seq 1 $(REPS)); do \
		echo -n "$$i\t$(P)\t"; \
		mpirun -np $(P) ./main_sendrecv -n $(N) -s $(SEED) -d $(DENSITY) -i $(ITERATIONS); \
	done

run_collectives: main_collectives
	@echo "reps\tnp\tn\tseed\tdensity\titers\tdimx\tdimy\timplementation\ttime (ms)\talive\tdead"
	@for i in $$(seq 1 $(REPS)); do \
		echo -n "$$i\t$(P)\t"; \
		mpirun -np $(P) ./main_collectives -n $(N) -s $(SEED) -d $(DENSITY) -i $(ITERATIONS); \
	done

run_sequential: main_sequential
	@echo "reps\tn\tseed\tdensity\titers\timplementation\ttime (ms)\talive\tdead"
	@for i in $$(seq 1 $(REPS)); do \
		echo -n "$$i\t"; \
		./main_sequential -n $(N) -s $(SEED) -d $(DENSITY) -i $(ITERATIONS); \
	done

clean:
	rm -f main_sendrecv main_collectives main_sequential *.o
