fc=gfortran
#fc=mpif90
myflgs = -O3
#myflgs = -fcheck=bounds -O0
#myflgs = -fno-automatic -fcheck=bounds -O3
# before compiling, please export LAPACK_LIB=path_of_lapack and export BLAS_LIB=path_of_blas
#lib = -L$(LAPACK_LIB) -llapack -L$(BLAS_LIB) -lblas  
lib = 

obj =   \
	prunesymm.o

all: compila

compila: $(obj)  
	$(fc) $(myflgs) $(obj) -o prunesymm.exe $(lib)

%.o : %.f
	$(fc) $(myflgs) -c $< 

%.o : %.f90
	$(fc) $(myflgs) -c $< 

clean: 
	rm -f *.o *.exe *.mod 





