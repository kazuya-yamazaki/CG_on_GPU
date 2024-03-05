CC	  =  mpicc
LIB_DIR  = 
INC_DIR  = 
OPTFLAGS = -fast -acc -gpu=cc80 -cuda -Minfo=accel -DOPENACC
LFLAGS   = 
#
TARGET = ../run/prog_gpu
default: $(TARGET)
OBJS =\
	test1.o\
	pfem_init.o\
	input_cntl.o\
	input_grid.o\
	define_file_name.o\
	mat_con0_atomic.o\
	mat_con1.o\
	mat_ass_main.o\
	mat_ass_bc.o\
	solve11.o\
	solver_CG.o\
	solver_SR.o\
	output_ucd.o\
	pfem_finalize.o\
	allocate.o\
	util.o

$(TARGET):  $(OBJS)
	$(CC) $(OPTFLAGS) -o $@ $(OBJS) $(LFLAGS)
.c.o:
	$(CC) $(OPTFLAGS) -c  $*.c
clean:
	/bin/rm -f *.o $(TARGET) *~ *.mod
