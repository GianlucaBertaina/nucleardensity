#!/bin/bash
####################################################################################################
#    Semiclassical Nuclear Densities. Full notice in LICENSE file. See also README.md
#    Copyright (C) 2020 Chiara Donatella Aieta, Marco Micciarelli, Gianluca Bertaina, Michele Ceotto
####################################################################################################

# Script for launching the Monte Carlo program for the evaluation of the nuclear densities
# of protonated Glycine in its ground state (ZPE) and excited OH stretch (OH) states
# at the harmonic and semiclassical anharmonic levels
# 
# INSTRUCTIONS:
# 1) Modify the PARAMETERS section below
# 2) Assign execution permission: chmod +x evaluate.sh
# 3) Launch script: ./evaluate.sh
# 4) Open VMD scripts in visualize folder with, for example: 
#    cd visualize; vmd -e compare_anharmonic-harmonic_ZPE.vmd


###  PARAMETERS  ################################################################################

# Insert number of available CPU cores
export NPROC=2

# Insert number of Monte Carlo samples (100000000 in the manuscript, 10000 for testing)
export NSTEPS=10000

# Insert (odd) number of intervals per direction (161 in the manuscript, 61 for testing)
export NBIN=61

#################################################################################################


export BASE=$(pwd)
export natoms=11
export DIRDENSITY=${BASE}"/../src/density"
export DIRCRUDE=${BASE}"/../src/crude"

pushd $DIRDENSITY
make all
popd

pushd $DIRCRUDE
make all
popd

# Monte Carlo evaluation of nuclear densities
for level in anharmonic harmonic; do
  pushd $level
  for state in ZPE OH; do
    pushd $state
    rm *out *cube *gz expect*
    # Sets parameters
    sed -e "7s/^.*$/${NBIN} ${NBIN} ${NBIN}/" -e "9s/^.*$/${NSTEPS}/" -i input_nucleardensity.dat
    # Evaluates ddensities
    mpirun -n $NPROC $DIRDENSITY/nucleardensity |tee log.out
    popd
  done
  popd
done

# Evaluate total density as sum of marginal densities. Archive marginal densities
for level in anharmonic harmonic; do
  pushd $level
  for state in ZPE OH; do
    pushd $state
    # Sum
    sumfile=density_sum_${level}_${state}.cube
    cp nucl_1.cube ${sumfile}
    for i in $(seq 2 ${natoms}); do
      $DIRCRUDE/crude comb ${sumfile} nucl_${i}.cube 1.0 1.0
      mv cuberes.cube  ${sumfile}
    done
    # Archive
    tar -czvf nucl.tar.gz nucl_* ${sumfile}
    rm nucl_*.cube
    popd
  done
  popd
done

# Evaluate differences of total densities for visualization
pushd visualize    

# Harmonic OH minus Harmonic ZPE
level=harmonic
$DIRCRUDE/crude comb ${BASE}/$level/OH/density_sum_${level}_OH.cube ${BASE}/$level/ZPE/density_sum_${level}_ZPE.cube 1.0 -1.0; mv cuberes.cube density_diff_OH-ZPE_${level}.cube

# Anharmonic OH minus Anharmonic ZPE
level=anharmonic
$DIRCRUDE/crude comb ${BASE}/$level/OH/density_sum_${level}_OH.cube ${BASE}/$level/ZPE/density_sum_${level}_ZPE.cube 1.0 -1.0; mv cuberes.cube density_diff_OH-ZPE_${level}.cube

# Anharmonic ZPE minus Harmonic ZPE
state=ZPE
$DIRCRUDE/crude comb ${BASE}/anharmonic/$state/density_sum_anharmonic_${state}.cube ${BASE}/harmonic/$state/density_sum_harmonic_${state}.cube 1.0 -1.0; mv cuberes.cube density_diff_anharmonic-harmonic_${state}.cube

# Anharmonic OH minus Harmonic OH
state=OH
$DIRCRUDE/crude comb ${BASE}/anharmonic/$state/density_sum_anharmonic_${state}.cube ${BASE}/harmonic/$state/density_sum_harmonic_${state}.cube 1.0 -1.0; mv cuberes.cube density_diff_anharmonic-harmonic_${state}.cube
popd

exit 0

