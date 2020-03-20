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
and edit the script file

```
evaluate.sh
```

as needed, then follow the instructions inside the file. Typical parameters to be changed are the number of available CPU cores on the computer, the number of Monte Carlo steps and the number of bins for each direction.

Open the resulting density diffrences with VMD using the VMD scripts in the 

```
run_Glyp/visualize
```
directory, as explained in the script file. Modify the Isovalue parameter in VMD opening the Graphics/Representations menu.

## Authors

* **Marco Micciarelli** - *Wavefunction evaluation*
* **Chiara Donatella Aieta** - *Density evaluation and parallelization*
* **Gianluca Bertaina** - *Other observables and errorbars*
* **Michele Ceotto** - *Supervision*

## License

This project is copyrighted by the above Authors - see the [LICENSE](LICENSE) file for details

## Acknowledgments

* Fabio Gabas for running semiclassical simulations
* Riccardo Conte for useful discussions on the method 

