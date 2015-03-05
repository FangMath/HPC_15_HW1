EXECUTABLES = int_ring jacobi-mpi.c
COMPILER = mpicc
FLAGS = -O3 -Wall -lrt

all: $(EXECUTABLES)

int_ring: int_ring.c
	$(COMPILER) $(FLAGS) int_ring.c -o int_ring 

jacobi-mpi: jacobi-mpi.c
	$(COMPILER) $(FLAGS) jacobi-mpi.c -o jacobi-mpi

clean:
	rm -rf $(EXECUTABLES)
