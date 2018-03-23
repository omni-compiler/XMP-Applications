#
export NTQC_TOP=/data/ra000002/a03358/FY2017-2H/xmp/ntchem-mini.xmp.Rev3
export HOSTTYPE=linux64_xmp_omp_k_fx10
export LAPACK=
export BLAS=
export ATLAS=-SSL2BLAMP
export SCRATCH=/home/a03358/scr/ntchem
export PARALLEL=mpiomp
#
export TARGET=LINUX64
unset USE_MPI

# if you want to use MPICH, you can set the environmental variables as
# follos (see ./GA/README)
# 
# export MPI_USE=yes
# export MPI_INCLUDE=/usr/include
# export MPI_LIB=/usr/lib
# export LIBMPI=-lmpi

export LARGE_FILES=yes

