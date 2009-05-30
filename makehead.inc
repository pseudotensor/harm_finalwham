#
#
#echo "Make sure MPICH/GM params set in both global.h and makefile!"
#echo "Make sure MPICH/GM params set in both global.h and makefile!"
USEMPI=1

LAPACKLIB=lapack
BLASLIB=blas
F2CLIB=g2c  #-l$(F2CLIB)

PREP = prep
CMD=grmhd
CMD2=grmhdldouble
