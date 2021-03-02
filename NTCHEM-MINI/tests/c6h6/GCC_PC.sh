#!/bin/bash
date
hostname
MOLECULE=c6h6
pwd
NPROCS=1
#export OMP_NUM_THREADS=8
export OMP_NUM_THREADS=1
time mpirun -np $NPROCS ../../bin/rimp2.exe | tee ${MOLECULE}_rimp2.out
ls -go
#python ./check.py

