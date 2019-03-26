#fc=gfortran
fc=mpif90
myflgs = -fno-automatic -fcheck=bounds -W -O3
# before compiling, please export LAPACK_LIB=path_of_lapack and export BLAS_LIB=path_of_blas
#lib = -L$(LAPACK_LIB) -llapack -L$(BLAS_LIB) -lblas  
lib = 

obj =   \
	constants_mod.o\
	io_units_def_mod.o\
	input_defs_mod.o\
        manage_cubes.o\
	hermite.o\
	calc_wfn.o\
	nude.o

all: compila

compila: $(obj)  
	$(fc) $(myflgs) $(obj) -o nude.exe $(lib)

%.o : %.f
	$(fc) $(myflgs) -c $< 

%.o : %.f90
	$(fc) $(myflgs) -c $< 

clean: 
	rm -f *.o *.exe *.mod 





