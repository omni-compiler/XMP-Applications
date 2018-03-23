#!/bin/sh
#PJM --rsc-list "node=2x2x2"
#PJM --rsc-list "elapse=10:00"
#PJM --mpi "use-rankdir"
#PJM --stgin  "rank=* ../../src/modylas_mini %r:./"
#PJM --stgin  "rank=* ./wat111.mdconf        %r:./"
#PJM --stgin  "rank=* ./wat111.mdff.bin      %r:./"
#PJM --stgin  "rank=* ./wat111.mdxyz.bin     %r:./"
#PJM --stgout "rank=0 %r:./wat111.mdmntr      ./job%j/"
#PJM --stgout "rank=0 %r:./wat111.restart.bin ./job%j/"
#PJM --stgout "rank=0 %r:./wat111.mdtrj.bin   ./job%j/"
#PJM -S

. /work/system/Env_base

export OMP_NUM_THREADS=8
export PARALLEL=8
export FORT90L='-Wl,-T'

mpiexec ./modylas_mini wat111
