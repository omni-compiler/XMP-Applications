#! /bin/bash -x

GLEV=${1}
RLEV=${2}
NMPI=${3}
ZL=${4}
VGRID=${5}
TOPDIR=${6}
BINNAME=${7}
RUNCONF=${8}

# System specific
MPIEXEC="mpiexec"

GL=`printf %02d ${GLEV}`
RL=`printf %02d ${RLEV}`
if   [ ${NMPI} -ge 10000 ]; then
	NP=`printf %05d ${NMPI}`
elif [ ${NMPI} -ge 1000 ]; then
	NP=`printf %04d ${NMPI}`
elif [ ${NMPI} -ge 100 ]; then
	NP=`printf %03d ${NMPI}`
else
	NP=`printf %02d ${NMPI}`
fi

dir2d=gl${GL}rl${RL}pe${NP}
dir3d=gl${GL}rl${RL}z${ZL}pe${NP}
res2d=GL${GL}RL${RL}
res3d=GL${GL}RL${RL}z${ZL}

MNGINFO=rl${RL}-prc${NP}.info

if [ ${NMPI} -gt 36864 ]; then
   rscgrp="huge"
elif [ ${NMPI} -gt 384 ]; then
   rscgrp="large"
else
   rscgrp="small"
fi

outdir=${dir3d}
cd ${outdir}
HERE=${PWD}

ln -s ${TOPDIR}/bin/${BINNAME} .
ln -s ${TOPDIR}/data/mnginfo/${MNGINFO} .
ln -s ${TOPDIR}/data/grid/vgrid/${VGRID} .

for f in $( ls ${TOPDIR}/data/grid/boundary/${dir2d} )
do
   ln -s ${TOPDIR}/data/grid/boundary/${dir2d}/${f} .
done


cat << EOF1 > run.sh
#! /bin/bash -x
################################################################################
#
# for Fugaku computer
#
################################################################################
#PJM --rsc-list "node=${NMPI}"
#PJM --rsc-list "elapse=02:00:00"
#PJM --mpi "use-rankdir"
#PJM -j
#PJM -s
#
export PARALLEL=8
export OMP_NUM_THREADS=8

# run
${MPIEXEC} ./${BINNAME} || exit

################################################################################
EOF1

exit
