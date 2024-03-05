#!/bin/sh
#PJM -L rscgrp=short-a
#PJM -L node=1
#PJM --mpi proc=1
#PJM -L elapse=00:15:00
#PJM -g your_group_name
#PJM -j

module load nvidia nvmpi cuda

mpirun -n 1 ./pmesh

rm wk.*
