#!/bin/bash

# Insert number of available CPU cores
export NPROC=20

export natoms=11
export DIRDENSITY=$(pwd)"/../src/density"
export DIRCRUDE=$(pwd)"/../src/crude"

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
    mpirun -n $NPROC $DIRDENSITY/nucleardensity |tee log.out
    popd
  done
  popd
done

# Evaluate totale density as sum of marginal densities. Archive marginal densities
for level in anharmonic harmonic; do
  pushd $level
  for state in ZPE OH; do
    pushd $state
    sumfile=density_sum_${level}_${state}.cube
    cp nucl_1.cube ${sumfile}
    for i in $(seq 2 ${natoms}); do
      $DIRCRUDE/crude comb ${sumfile} nucl_${i}.cube 1.0 1.0
      mv cuberes.cube  ${sumfile}
    done
    tar -czvf nucl.tar.gz nucl_* ${sumfile}
    rm nucl_*.cube
    popd
  done
  popd
done

exit 0

# archive cubes for wf
for j in ${states}; do pushd wf-${j}-symm-gs-pruned;tar -czvf nucl.tar.gz nucl_* ; popd; done
for j in ${states}; do pushd wf-${j}-symm-gs-pruned;rm nucl_{1..11}.cube ; popd; done
