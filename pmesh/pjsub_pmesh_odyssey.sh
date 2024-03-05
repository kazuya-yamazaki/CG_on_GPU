#!/bin/sh
#PJM -L rscgrp=short-o
#PJM -L node=1
#PJM --mpi proc=1
#PJM -L elapse=00:15:00
#PJM -g your_group_name
#PJM -j

module load fj
module load fjmpi

mpiexec ./pmesh

rm wk.*
