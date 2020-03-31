# Semiclassical Nuclear Densities

This code evaluates the nuclear density distributions, together with the bond, angle, and dihedral distributions, given the wavefunctions determined with a semiclassical method [https://doi.org/10.1063/1.5041911](https://doi.org/10.1063/1.5041911), using a Monte Carlo approach. 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for testing purposes and reproduction of the distributions for protonated Glycine.

### Prerequisites

These instructions are given for a Linux system. Both gfortran, MPICH, and VMD are needed.

### Installing and running the tests

Enter the directory 

```
run_Glyp
```
which contains the script for launching the Monte Carlo program for the evaluation of the nuclear densities of protonated Glycine in its ground state (ZPE) and excited OH stretch (OH) states at the harmonic and semiclassical anharmonic levels.

Edit the script file

```
evaluate.sh
```

and modify the initial PARAMETERS section as needed. In particular, please provide the number of available CPU cores to be dedicated to the parallel Monte Carlo calculation, the number of Monte Carlo samples to be drawn, and the number of voxels for each Cartesian direction in the output cube files.

Default values are provded in the script, which allow for fast testing on a laptop, while the values used in the manuscript are indicated (which imply a calculation duration of a few hours depending on the computational resources).


Assign execution permission and launch the script

```
chmod +x evaluate.sh
./evaluate.sh
```

Open the resulting density diffrences with VMD using the VMD scripts in the 

```
run_Glyp/visualize
```
directory, with the command (as an example):

```
cd run_Glyp/visualize
vmd -e compare_anharmonic-harmonic_ZPE.vmd
```
In VMD, modify the Isovalue parameter opening the Graphics/Representations menu.

## Authors

* **Marco Micciarelli** - *Wavefunction and density evaluation*
* **Chiara Donatella Aieta** - *Density evaluation and parallelization*
* **Gianluca Bertaina** - *Other observables and errorbars*
* **Michele Ceotto** - *Supervision*

## License

This project is copyrighted by the above Authors - see the [LICENSE](LICENSE) file for details

## Acknowledgments

* Fabio Gabas for running semiclassical simulations
* Riccardo Conte for useful discussions on the method 

