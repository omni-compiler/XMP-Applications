#!/bin/bash -x
#
#PJM -o "c6h6_rimp2.out"
#PJM -e "c6h6_rimp2.err"
#PJM --rsc-list "node=2"
#PJM --rsc-list "elapse=0:10:00"
#PJM -s
#

export FLIB_FASTOMP=FALSE
export FLIB_CNTL_BARRIER_ERR=FALSE

mpiexec ./rimp2.exe
