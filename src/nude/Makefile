fc=gfortran
myflgs = -fno-automatic -fcheck=bounds -W -O3
lib = -L/home/marcom/LAPACK/lapack-3.6.0 -llapack -L/home/marcom/BLAS/BLAS-3.6.0 -lblas -llapack -lblas
#lib = -lblas -llapack 

obj =   \
	constants_mod.o\
	io_units_def_mod.o\
	input_defs_mod.o\
        manage_cubes.o\
	hermite.o\
	calc_wfn.o\
	nude.o

compila: $(obj)  
	$(fc) $(myflgs) $(obj) -o nude.exe $(lib)

%.o : %.f
	$(fc) $(myflgs) -c $< 

%.o : %.f90
	$(fc) $(myflgs) -c $< 

clean: 
	rm -f *.o *.exe *.mod 





