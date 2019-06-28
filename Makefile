MPI_PROGRAMS = nbody_mpi_collective.x
MPI_OBJECTS = $(MPI_PROGRAMS:.x=.o)
TEST_PROGRAMS = xstring.x

PROGRAMS = $(MPI_PROGRAMS) $(TEST_PROGRAMS)
OBJECTS = $(MPI_OBJECTS) xdmf_write.o xstring.o

CC = gcc
CC_ARGS = -std=c99 -I./
CC_ARGS += -Ofast -march=native -Wall -Werror -Wpedantic -Wno-unknown-pragmas

MPICC = mpicc
MPICC_ARGS = $(CC_ARGS) -DTEST_SCRATCH

all: $(PROGRAMS)

$(PROGRAMS): %.x: xdmf_write.o xstring.o

$(MPI_PROGRAMS): %.x: %.o
	$(MPICC) $(MPICC_ARGS) -o $@ $^ -lm

$(MPI_OBJECTS): %.o: %.c
	$(MPICC) -c $(MPICC_ARGS) -o $@ $<

$(TEST_PROGRAMS): %.x: %.c
	$(CC) -DTEST_MAIN $(CC_ARGS) -o $@ $< -lm

xdmf_write.o: %.o : %.c
	$(CC) -c $(CC_ARGS) -o $@ $<

xstring.o: %.o : %.c
	$(CC) -c $(CC_ARGS) -o $@ $<

exec:
	mpirun -np 3 $(MPI_PROGRAMS) kaplan_1000.bin

.PHONY: clean

clean:
	-rm $(PROGRAMS) $(OBJECTS)
