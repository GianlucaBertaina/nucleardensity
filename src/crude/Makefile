fc=gfortran
#fc=mpif90
myflgs = -O3
#myflgs = -fcheck=bounds -O0
#myflgs = -fno-automatic -fcheck=bounds -O3
lib = 

obj =   \
	crude.F90

all: compile

compile: $(obj)  
	$(fc) $(myflgs) $(obj) -o crude $(lib)

%.o : %.f
	$(fc) $(myflgs) -c $< 

%.o : %.f90
	$(fc) $(myflgs) -c $< 

clean: 
	rm -f *.o crude *.mod 





