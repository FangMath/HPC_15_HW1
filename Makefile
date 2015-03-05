EXECUTABLES = sr sr2 pingpong pingpong_array pingpong_nonblocking
COMPILER = mpicc
FLAGS = -O3 -Wall

all: $(EXECUTABLES)

sr: sr.c
	$(COMPILER) $(FLAGS) sr.c -o sr

sr2: sr2.c
	$(COMPILER) $(FLAGS) sr2.c -o sr2

pingpong: pingpong.c
	$(COMPILER) $(FLAGS) $^ -o pingpong

pingpong_array: pingpong_array.c
	$(COMPILER) $(FLAGS) $^ -o pingpong_array

pingpong_nonblocking: pingpong_nonblocking.c
	$(COMPILER) $(FLAGS) $^ -o pingpong_nonblocking

clean:
	rm -rf $(EXECUTABLES)
