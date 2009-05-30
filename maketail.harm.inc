
#SRCS = \
#bounds.c boundsint.c coord.c diag.c dudp_calc_3vel.c dudp_calc.c dump.c dump_ener.c fixup.c gaussj.c \
#image.c initbase.c init.c interppoint.c lubksb.c ludcmp.c main.c metric.c \
#mnewt.c nrutil.c phys.c ranc.c restart.c set_arrays.c set_grid.c \
#step_ch.c tensor.c utoprim.c utoprim_2d.c utoprim_1d.c \
#utoprim_1d_opt.c vchar.c transforms.c fail.c utoprim_ldz.c \
#utoprim_jon.c \
#mpi_init.c mpi_fileio.c boundmpi.c boundmpiint.c sources.c utoprim_1d_final.c \
#utoprim_2d_final.c utoprim_5d2_final.c phys.ffde.c \
#advance.c reconstructeno.c utoprimgen.c flux.c dumpgen.c boundsflux.c \
#interpline.c monointerp.c fluxct.c fluxcompute.c smonointerp.c \
#jmonointerp.c wavespeeds.c phys.coldgrmhd.c eos.c

#OBJS = \
#bounds.o boundsint.o coord.o diag.o dudp_calc_3vel.o dudp_calc.o dump.o dump_ener.o fixup.o gaussj.o \
#image.o initbase.o init.o interppoint.o lubksb.o ludcmp.o main.o metric.o \
#mnewt.o nrutil.o phys.o ranc.o restart.o set_arrays.o set_grid.o \
#step_ch.o tensor.o utoprim.o utoprim_2d.o utoprim_1d.o \
#utoprim_1d_opt.o vchar.o transforms.o fail.o utoprim_ldz.o \
#utoprim_jon.o \
#mpi_init.o mpi_fileio.o boundmpi.o boundmpiint.o sources.o utoprim_1d_final.o \
#utoprim_2d_final.o utoprim_5d2_final.o phys.ffde.o \
#advance.o reconstructeno.o utoprimgen.o flux.o dumpgen.o boundsflux.o \
#interpline.o monointerp.o fluxct.o fluxcompute.o smonointerp.o \
#jmonointerp.o wavespeeds.o phys.coldgrmhd.o eos.o

SRCSMANYPOINTS = \
init.c initbase.c restart.c set_arrays.c set_grid.c \
advance.c reconstructeno.c flux.c wavespeeds.c boundsflux.c \
fluxct.c fluxcompute.c \
bounds.c boundsint.c \
diag.c dump.c dumpgen.c dump_ener.c image.c \
\
fail.c fixup.c \
mpi_init.c mpi_fileio.c boundmpi.c boundmpiint.c \
interppoint.c interpline.c monointerp.c jmonointerp.c smonointerp.c


OBJSMANYPOINTS = \
init.o initbase.o restart.o set_arrays.o set_grid.o \
advance.o reconstructeno.o flux.o wavespeeds.o boundsflux.o \
fluxct.o fluxcompute.o \
bounds.o boundsint.o \
diag.o dump.o dumpgen.o dump_ener.o image.o \
\
fail.o fixup.o \
mpi_init.o mpi_fileio.o boundmpi.o boundmpiint.o \
interppoint.o interpline.o monointerp.o jmonointerp.o smonointerp.o \

# phys.c does currents across multiple points, but otherwise one pointed
# wavespeeds kinda at more than one point

# for some reason main.c really matters alot.  Needs to be compiled
# withOUT -mp -pc64 for precision in inversion to be ok  ODD!  GODMARK

SRCSONEPOINT = \
\
main.c metric.c coord.c phys.c step_ch.c \
phys.ffde.c phys.coldgrmhd.c eos.c \
vchar.c transforms.c sources.c \
gaussj.c lubksb.c ludcmp.c mnewt.c nrutil.c ranc.c tensor.c \
utoprimgen.c dudp_calc_3vel.c dudp_calc.c \
utoprim.c utoprim_2d.c utoprim_1d.c utoprim_1d_opt.c utoprim_ldz.c \
utoprim_jon.c utoprim_1d_final.c utoprim_2d_final.c utoprim_5d2_final.c

OBJSONEPOINT = \
\
main.o metric.o coord.o phys.o step_ch.o \
phys.ffde.o phys.coldgrmhd.o eos.o \
vchar.o transforms.o sources.o \
gaussj.o lubksb.o ludcmp.o mnewt.o nrutil.o ranc.o tensor.o \
utoprimgen.o dudp_calc_3vel.o dudp_calc.o \
utoprim.o utoprim_2d.o utoprim_1d.o utoprim_1d_opt.o utoprim_ldz.o \
utoprim_jon.o utoprim_1d_final.o utoprim_2d_final.o utoprim_5d2_final.o







