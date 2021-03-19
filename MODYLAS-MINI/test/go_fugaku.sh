#!/bin/sh
#PJM -L "node=8"
#PJM -L "elapse=00:15:00"
#PJM -S

MODYLAS=../src/modylas_mini

export OMP_NUM_THREADS=2

(mpiexec $MODYLAS ./wat111 \
  && python2 ./check.py wat111.mdmntr.8n_8t wat111.mdmntr) 2>&1 | tee log
