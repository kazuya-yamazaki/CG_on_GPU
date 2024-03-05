#!/bin/sh
#PJM -L rscgrp=short-a
#PJM -L node=1
#PJM --mpi proc=1
#PJM -L elapse=00:15:00
#PJM -g your_group_number
#PJM -j

module load nvidia nvmpi cuda

mpirun -n 1 ./prog_gpu
# To use profiler: mpirun -n 1 nsys profile -o report_%q{OMPI_COMM_WORLD_RANK} --trace=cuda,nvtx --sample=none --cpuctxsw=none --force-overwrite true ./prog_gpu

