#!/bin/sh

PROG=../bin/ffb_mini

mpiexec -np 8 -x XMP_ONESIDED_HEAP_SIZE=8192M $PROG 2 2 2 46 | tee les3x.log.P0001 \
 && ./check.py master/les3x.log les3x.log
