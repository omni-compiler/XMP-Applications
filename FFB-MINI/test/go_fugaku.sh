#!/bin/sh
#
#  pjsub --interact --sparam "wait-time=600" ./go_fugaku.sh
#
#PJM -L "node=8"
#PJM -L "elapse=00:05:00"
#PJM -S

PROG=../bin/ffb_mini
export PARALLEL=8
export XMP_ONESIDED_HEAP_SIZE=8192M

mpiexec -np 8 ${PROG} 2 2 2 46 color_partsize=2000 reorder_ndiv=10 unroll=on \
 | tee les3x.log.P0001 
