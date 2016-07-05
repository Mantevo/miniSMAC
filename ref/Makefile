EXE = smac2d_mpi_and_threads_with_AVX

FC=mpif90
#FC=ifort

# for MPI only, double precision
#FFLAGS=-O3 -c -g -DD_PRECISION -check bounds -traceback
# for MPI + OpenMP, double precision
FFLAGS=-O3 -xAVX -c -g -DD_PRECISION -check bounds -traceback -openmp \
      -I/opt/intel-12.1/compiler/include

# for linux machines
#CPP=/lib/cpp
# for Mac OSX
#CPP=/opt/local/bin/cpp
# for TLCC2 linux clusters like glory, chama
CPP=/usr/bin/cpp

# preprocessor flags, with double precision
# ... for MPI only
#CPPFLAGS=-P -traditional -DD_PRECISION
# ... for MPI + OpenMP, double precision
CPPFLAGS=-P -traditional -DD_PRECISION -D_OPENMP \
      -I/opt/intel-12.1/compiler/include

# linker
LD=mpif90

# for MPI only
#LDFLAGS=-O3
# for MPI + OpenMP
LDFLAGS=-O3 -openmp \
      -I/opt/intel-12.1/compiler/include

INCDIR=./include

# for linux machines
#MPIDIR=-I/home/projects/openmpi/1.7.3/gnu/4.7.2/include
# for Mac OSX
#MPIDIR=-I/Users/dwbarne/local/packages/OpenMPI-1.7.2/include
# for redsky
#MPIDIR=-I/apps/x86_64/mpi/openmpi/intel-12.1-2011.13.367/openmpi-1.4.3_oobpr/include
# for glory, chama
MPIDIR=-I/apps/x86_64/mpi/openmpi/intel-12.1-2011.13.367/openmpi-1.4.3/include

# delete all previous suffixes
.SUFFIXES:
# assign new suffixes for defining transformation rule '.f.o'
.SUFFIXES: .f .o

OBJ =   bc.o \
	btlu3.o \
        flux.o \
        geom.o \
        initia.o \
        lreltur.o \
        metric.o \
        rbc.o \
        smac2d.o \
        tursa.o \
        blas.o \
        conup.o \
        fm.o \
	ilusol.o \
        lrelax.o \
        lusgs.o \
        preproc.o \
        rgrid.o \
        step.o \
        visc.o

#---------------- compilation rules ----------------------
.f.o:
	$(CPP) $(CPPFLAGS) $*.f -I$(INCDIR) $(MPIDIR) > $*F.f
	$(FC) $(FFLAGS) $(*F)F.f -o $*.o
#	rm $*F.f

$(EXE): $(OBJ)
	$(LD) $(LDFLAGS) -o $(EXE) $(OBJ)

# clean up
clean:
	rm -f *F.f *.o *.f~ *F.f.*
