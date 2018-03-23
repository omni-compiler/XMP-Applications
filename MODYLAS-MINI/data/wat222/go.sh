#!/bin/sh
#PJM --rsc-list "node=4x4x4"
#PJM --rsc-list "elapse=0:15:00"
#PJM --rsc-list="rscgrp=small"
#PJM --mpi "use-rankdir"
#PJM --stgin  "rank=* ../../src/modylas_mini %r:./"
#PJM --stgin  "rank=* ./wat222.mdconf        %r:./"
#PJM --stgin  "rank=* ./wat222.mdff.bin      %r:./"
#PJM --stgin  "rank=* ./wat222.mdxyz.bin     %r:./"
#PJM --stgout "rank=0 %r:./wat222.mdmntr      ./job%j/"
#PJM --stgout "rank=0 %r:./wat222.restart.bin ./job%j/"
#PJM --stgout "rank=0 %r:./wat222.mdtrj.bin   ./job%j/"
#PJM -S

. /work/system/Env_base

NTHREADS=8
export PARALLEL=${NTHREADS}
export OMP_NUM_THREADS=${NTHREADS}
LPG="/opt/FJSVxosmmm/sbin/lpgparm -t 4MB -s 4MB -h 4MB -d 4MB -p 4MB"


# Endian
export FORT90L='-Wl,-T'

# MCA Parameters
MCA_PARAM="--mca common_tofu_fastmode_threshold 0"
MCA_PARAM="${MCA_PARAM} --mca common_tofu_max_fastmode_procs 40"

LD="./modylas_mini ./wat222"

mpiexec ${MCA_PARAM} ${LPG} ${LD}
