# Semiclassical Nuclear Densities

This code evaluates the nuclear density distributions, together with the bond, angle, and dihedral distributions, given the wavefunctions determined with a semiclassical method [https://doi.org/10.1063/1.5041911](https://doi.org/10.1063/1.5041911), using a Monte Carlo approach. 

If you use this code and/or the provided datasets for scientific publication, please cite the following articles:

C. Aieta, G. Bertaina, M. Micciarelli, and M. Ceotto "Representing Molecular Ground and Excited Vibrational Eigenstates with Nuclear Densities obtained from Semiclassical Initial Value Representation Molecular Dynamics"

C. Aieta, M. Micciarelli, G. Bertaina and M. Ceotto "Anharmonic quantum nuclear densities from full dimensional vibrational eigenfunctions with application to protonated glycine", Nat. Commun. (2020) [https://doi.org/10.1038/s41467-020-18211-](https://doi.org/10.1038/s41467-020-18211-).

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for testing purposes and reproduction of the distributions for water and protonated Glycine.

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

## IMPORTANT NOTICE
Indexes of normal modes are in their energetic order, so for water:

(bending, symm. stretching, asymm. stretching) -> BSA order, not SBA order, which is used in the article.

## Authors

* **Marco Micciarelli** - *Wavefunction and density evaluation*
* **Chiara Donatella Aieta** - *Density evaluation and parallelization*
* **Gianluca Bertaina** - *Other observables and errorbars*
* **Michele Ceotto** - *Supervision*

## License

This project is copyrighted by the above Authors - see the [LICENSE](LICENSE) file for details

## Acknowledgments

* Fabio Gabas for running semiclassical simulations
* Jaime Suarez for providing the results of the DVR calculations
* Riccardo Conte for useful discussions on the method 

