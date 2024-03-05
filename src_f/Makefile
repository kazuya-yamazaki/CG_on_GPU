F90	  =  mpif90
F90LINKER = $(F90)
LIB_DIR   = 
INC_DIR   = 
OPTFLAGS  = -fast -acc -gpu=cc80,O3 --ptxas-options=-v -cuda -Minfo=accel
FFLAGS = -DOPENACC -cpp $(OPTFLAGS)
FLIBS = 
F90LFLAGS= -lnvhpcwrapnvtx
#
TARGET = ../run/prog_gpu
default: $(TARGET)
OBJS =\
pfem_util.o nvtxwrap.o \
solver_SR.o solver_CG.o \
solver11.o test1.o util.o pfem_init.o input_cntl.o input_grid.o define_file_name.o\
mat_con0.o mat_con1.o mat_ass_main.o mat_ass_bc.o pfem_finalize.o output_ucd.o

$(TARGET):  $(OBJS)
	$(F90LINKER) $(OPTFLAGS) -o $(TARGET) $(OBJS) $(F90LFLAGS)
clean:
	/bin/rm -f *.o $(TARGET) *~ *.mod
.f.o:
	$(F90) $(FFLAGS) $(INC_DIR) -c  $*.f
.f90.o:
	$(F90) $(FFLAGS) $(INC_DIR) -c  $*.f90
.SUFFIXES: .f90 .f
