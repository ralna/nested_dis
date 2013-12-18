# Possible Fortran compilers with options

ifort = /home/opt/intel/Compiler/11.1/075/bin/intel64/ifort
ifort = ifort

f90= nagfor -dcfuns -C=all -C=undefined -gline -f95
f90= $(ifort) -g -C -static -openmp
f90= $(ifort) -C -check noarg_temp_created -u -debug extended -traceback -openmp
#f90 = g95 -O2 -Wall -pedantic -fbounds-check -ftrace=full -g
f90 = g95   -fintrinsic-extensions -std=f2003 -Wimplicit-none -ftrace=full 
# Note: g95 needs LIBS = -L/usr/local/lib/ -lmetis_long

f90 = gfortran  -Wall -pedantic -fbounds-check -march=native -g -Wunused -fbacktrace -fopenmp -fprofile-arcs -ftest-coverage
#f90 = /numerical/local_beale/bin/pgfortran -gopt -Mbounds -Mchkptr -Mchkstk -Mstandard -Mallocatable=03
#f90 = /usr/local/lf6481/bin/lf95 -g --chk --sse3

# Options for speed
#f90 = gfortran 	-O3 -msse3 -march=native -g -pg
#f90 = pgfortran -fast -O4
#f90 = ifort -O -march=core2

FC = $(f90)

# packages directory
pd = /numerical/num/hsl2013/packages/

# blas files
DBLAS =
BLAS =

# Link to libraries

# This works on curtis
LIBS =  -L /opt/intel/Compiler/11.1/072/mkl/lib/em64t -lmkl_gf_lp64 -lmkl_sequential -lmkl_core \
        -L/usr/local/metis-4.0   -lmetis	



LIBS = -L/usr/local/GotoBLAS -lgoto_prof -lgoto  \
        -L/usr/local/metis-4.0   -lmetis
        
LIBS = -L/usr/local/GotoBLAS -lgoto  \
        -L/usr/local/metis-4.0   -lmetis


# This works on irons
LIBS = -L/opt/intel/mkl/lib/intel64/ -lmkl_gf_lp64 -lmkl_sequential -lmkl_core   \
        -L/usr/local/metis-4.0.3   -lmetis


# This works on fox
LIBS = -L/opt/intel/mkl/10.1.0.015/lib/em64t -lmkl_gf_lp64 -lmkl_sequential -lmkl_core  \
        -L/usr/local/metis-4.0   -lmetis


# This works on mitchell
LIBS = -L/opt/intel/mkl/lib/intel64/ -lmkl_gf_lp64 -lmkl_sequential -lmkl_core   \
        -L/usr/local/metis-4.0   -lmetis


# Dependencies


dtdeps=  $(fd15)fd15d.f $(fa14)fa14d.f  $(ym11)ym11d.f  $(mc54)mc54d.f \
         $(mc59)mc59d.f $(kb07)kb07i.f

ddeps90=$(hsl_mc79)hsl_mc79i.f90 $(hsl_mc78)hsl_mc78i.f90 \
          $(hsl_fa14)hsl_fa14d.f90 maxflow.f90

dtdeps90 = $(hsl_mc69)hsl_mc69d.f90


# What to run 
all: hsl_mc70dt

# Clean up directory
clean:
	rm a.out *.o  temp *deps* *.mod *deps90.f90 fort.8 coco.log *deps901.f90 *gcda *gcno

# Individual progs that can be run

hsl_mc70ds: ddeps90.o  hsl_mc70d.o hsl_mc70ds.f90  
	 $(f90) ddeps90.o  hsl_mc70d.o hsl_mc70ds.f90  $(LIBS)
	 ./a.out < hsl_mc70ds.data

hsl_mc70dt: ddeps90.o dtdeps90.o dtdeps.o 	hsl_mc70d.o hsl_mc70dt.o  
	 $(f90)  ddeps90.o dtdeps90.o dtdeps.o 	hsl_mc70d.o hsl_mc70dt.o  $(LIBS)
	 ./a.out  < hsl_mc70dt.data 


# Compile components needed for tests
ddeps90.f90:
	cat  $(ddeps90) >ddeps90.f90
	echo $(ddeps90) >ddeps90

dtdeps90.f90:
	cat  $(dtdeps90) >dtdeps90.f90
	echo $(dtdeps90) >dtdeps90
dtdeps.f:
	cat  $(dtdeps) >dtdeps.f
	echo $(dtdeps) >dtdeps
 
hsl_mc70d.o: ddeps90.o 	hsl_mc70d.f90
	$(f90) -c  ddeps90.o hsl_mc70d.f90
hsl_mc70ds.o:	hsl_mc70ds.f90
	$(f90) -c hsl_mc70ds.f90
hsl_mc70ds1.o:	hsl_mc70ds1.f90
	$(f90) -c hsl_mc70ds1.f90
hsl_mc70dt.o: hsl_mc70dt.f90
	$(f90) -c hsl_mc70dt.f90
ddeps90.o:	ddeps90.f90
	$(f90) -c ddeps90.f90
dtdeps90.o:	dtdeps90.f90
	$(f90) -c dtdeps90.f90
dtdeps.o:	dtdeps.f
	$(f90) -c dtdeps.f

# Directories containing different codes
hsl_fa14 = $(pd)hsl_fa14/
hsl_mc78 = $(pd)hsl_mc78/
hsl_mc79 = $(pd)hsl_mc79/
hsl_mc69 = $(pd)hsl_mc69/
fa14 = $(pd)fa14/
fd15 = $(pd)fd15/
kb07 = $(pd)kb07/
mc54 = $(pd)mc54/
mc59 = $(pd)mc59/
ym11 = $(pd)ym11/
