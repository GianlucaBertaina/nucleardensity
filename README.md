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

The "run_Glyp" and "run_H2O" directories contain contain the scripts for launching the Monte Carlo program for the evaluation of the nuclear densities of protonated Glycine and of water in their ground state (ZPE) and excited states at the harmonic and semiclassical anharmonic levels.

In particular, both directories contain subdirectories referring to the level of theory: "harmonic", "anharmonic" (semiclassical) and "DVR" (only for water. The deepest directory level refers to the state of interest: ground state (ZPE) and normal-mode related states (see article). 

For water, we evaluate densities for the following excited states:

```
100: bending
200: overtone of bending
010: symmetric stretch
001: asymmetric stretch
```
*IMPORTANT NOTICE*: Indexes of normal modes are in their energetic order, so for water:

(bending, symm. stretching, asymm. stretching) -> BSA order, not standard SBA order, which is used in the article.

For protonated Glycine, we evaluate densities for the following excited states:

```
23: N-H3 stretch + C2-H4/5 symmetric stretch out of phase
25: N-H1/2 symmetric stretch
26: N-H1/2 asymmetric stretch
27: O1-H6 stretch
```

Choose the molecule, enter the corresponding "run_..." directory and edit the script file

```
evaluate.sh
```

and modify the initial PARAMETERS section as needed. In particular, please provide the number of available CPU cores to be dedicated to the parallel Monte Carlo calculation, the number of Monte Carlo samples to be drawn, and the number of voxels for each Cartesian direction in the output cube files.

Default values are provided in the script, which allow for fast testing on a laptop, while the values used in the manuscript are indicated (which imply a calculation duration of a few hours depending on the computational resources).

Assign execution permission and launch the script

```
chmod +x evaluate.sh
./evaluate.sh
```

Open the resulting density diffrences with VMD using the VMD scripts in the "run_Glyp/visualize" or "run_H2O/visualize" directories, with the command (as an example):

```
cd run_Glyp/visualize
vmd -e compare_anharmonic-harmonic_ZPE.vmd
```
to plot the difference between the ZPE densities of protonated Glycine at the harmonic and anharmonic levels, or (as an example):

```
cd run_H2O/visualize
vmd -e plot_001_DVR.vmd
```
to plot the full density of state 001 for water, calculated at the DVR level.

In VMD, modify the Isovalue parameter opening the Graphics/Representations menu.

In the simulation directories, you also find, besides the cube nuclear density files, the bond, angle, and (for Glyp) dihedral distributions, which can be plotted with standard tools such as gnuplot.

## Authors

* **Marco Micciarelli** - *Wavefunction and density evaluation*
* **Chiara Donatella Aieta** - *Density evaluation and parallelization*
* **Gianluca Bertaina** - *Other observables and errorbars*
* **Michele Ceotto** - *Supervision*

## License

This project is copyrighted by the above Authors - see the [LICENSE](LICENSE) file for details

## Acknowledgments

* Fabio Gabas for running semiclassical simulations
* Jaime Suarez for providing the results of the DVR calculations, from [https://doi.org/10.1063/1.5041911](https://doi.org/10.1063/1.5041911)
* Riccardo Conte for useful discussions on the method 

