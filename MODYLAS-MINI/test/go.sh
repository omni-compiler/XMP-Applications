#!/bin/sh

MODYLAS=../src/modylas_mini

export OMP_NUM_THREADS=2

(mpiexec -n 8 $MODYLAS ./wat111 \
  && ./check.py wat111.mdmntr.8n_8t wat111.mdmntr) 2>&1 | tee log
