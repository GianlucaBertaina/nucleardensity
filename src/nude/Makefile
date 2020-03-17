#fc=gfortran
fc=mpif90
myflgs = -O3
#myflgs = -fno-automatic -fcheck=bounds -O3 
lib = 

obj =   \
	constants_mod.o\
	io_units_def_mod.o\
	input_defs_mod.o\
	manage_cubes.o\
	hermite.o\
	calc_wfn.o\
	nude.o

all: compile

compile: $(obj)  
	$(fc) $(myflgs) $(obj) -o nucleardensity $(lib)

%.o : %.f
	$(fc) $(myflgs) -c $< 

%.o : %.f90
	$(fc) $(myflgs) -c $< 

clean: 
	rm -f *.o nucleardensity *.mod 





