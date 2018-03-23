#!/bin/sh
#PJM --rsc-list "node=8x8x8"
#PJM --rsc-list "elapse=0:15:00"
#PJM --rsc-list="rscgrp=large"
#PJM --mpi "use-rankdir"
#PJM --stgin  "rank=* ../../src/modylas_mini %r:./"
#PJM --stgin  "rank=* ./wat444.mdconf        %r:./"
#PJM --stgin  "rank=* ./wat444.mdff.bin      %r:./"
#PJM --stgin  "rank=* ./wat444.mdxyz.bin     %r:./"
#PJM --stgout "rank=0 %r:./wat444.mdmntr      ./job%j/"
#PJM --stgout "rank=0 %r:./wat444.restart.bin ./job%j/"
#PJM --stgout "rank=0 %r:./wat444.mdtrj.bin   ./job%j/"
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

LD="./modylas_mini ./wat444"

mpiexec ${MCA_PARAM} ${LPG} ${LD}
