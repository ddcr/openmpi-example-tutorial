SHELL=/bin/bash

MODULE_LOAD_COMPILER = \
  module load $(1) || true; \
  echo === Compiler and version ===; \
  $(2) --version

MODULE_UNLOAD_COMPILER = \
  module unload $(1) || true; \
  echo === Compiler $(1) unloaded ===

MODULE_LOAD_PACKAGE = \
  module load $(1) || true; \
  echo === $(2) ===; \
  echo $${$(strip $(2))}

MODULE_UNLOAD_PACKAGE = \
  module unload $(1) || true; \
  echo === $(1) unloaded ===

COMPILER=
COMPILERNAME=intel
# COMPILERNAME=gnu       # gnu/5.4.0    (cluster veredas)
# COMPILERNAME=gnu7      # gnu7/7.2.0   (cluster veredas)
# COMPILERNAME=intel     # intel/14.0.1 (cluster veredas)

# MPINAME=mpi
# MPINAME=openmpi
# MPINAME=openmpi3
# MPINAME=mvapich2
MPINAME=impi

CC = gcc
CXX = g++
F77 = gfortran
FC = gfortran

MPICC = mpicc
MPICXX = mpicxx
MPIF77 = mpifort
MPIFC = mpifort

ifeq ("$(COMPILERNAME)", "intel")
  CC = icc
  CXX = icpc
  F77 = ifort
  FC = ifort
  MPICC = mpiicc
  MPICXX = mpiicpc
  MPIF77 = mpiifort
  MPIFC = mpiifort
  CC_ARGS = -std=c99 -I./
  CC_ARGS += -O3 -xHost
else
  CC_ARGS = -std=c99 -fgnu89-inline -I./
  CC_ARGS += -Ofast -march=native -Wall -Werror -Wpedantic -Wno-unknown-pragmas
endif
MPICC_ARGS = $(CC_ARGS) -DTEST_SCRATCH


ifeq ("$(COMPILERNAME)", "docker")
  # Use docker (my laptop)
  CROSSCOMPILE=./centos5-devtoolset2-gcc4-local
else ifeq ("$(COMPILERNAME)", "devtoolset-2")
  # Use RedHat software collections (cluster veredas)
  CROSSCOMPILE=scl enable devtoolset-2 --
else
  # Use modulefiles
  MODULE_LOAD_CC = $(call MODULE_LOAD_COMPILER, $(COMPILERNAME), $(CC))
  MODULE_UNLOAD_CC = $(call MODULE_UNLOAD_COMPILER, $(COMPILERNAME))
endif


MODULE_LOAD_MPICC = $(call MODULE_LOAD_PACKAGE, $(MPINAME), MPIHOME)
MODULE_UNLOAD_MPICC = $(call MODULE_LOAD_PACKAGE, $(MPINAME), MPIHOME)

MPI_PROGRAMS = nbody_mpi_collective.x
MPI_OBJECTS = $(MPI_PROGRAMS:.x=.o)
TEST_PROGRAMS = xstring.x

PROGRAMS = $(MPI_PROGRAMS) $(TEST_PROGRAMS)
OBJECTS = $(MPI_OBJECTS) xdmf_write.o xstring.o

all:
	( \
	  $(MODULE_LOAD_CC); \
	  $(MODULE_LOAD_MPICC); \
	  $(MAKE) $(PROGRAMS); \
	  $(MODULE_UNLOAD_MPICC); \
	  $(MODULE_UNLOAD_CC); \
	)

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
