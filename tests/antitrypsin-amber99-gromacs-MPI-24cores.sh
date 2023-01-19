#!/bin/bash

### ATTENTION: Hereafter, an example for executing MPI scripts on 48 cores in UniTn Cluster. Change it for your cluster/machine 

#PBS -l select=1:ncpus=24:mpiprocs=1:mem=4gb
#PBS -q VARIAMOLS_cpuQ

module load mpich-3.2

cd $PBS_O_WORKDIR
######################################################################

# Executing "block.py" and then "CANVAS.py" for Antitrypsin 
# NOTE: The output of block.py is the list of survived atoms, i.e. list-atoms-opt-3.txt, that will be the input of CANVAS.py script  


PYTHONDIR=../PYTHON-scripts

inputDIR=../input-files/ANTITRYPSIN/sim-w-amber99

python3 $PYTHONDIR/block-MPI.py choice3 -g $inputDIR/372_probs_TEMP.gro -l $inputDIR/list-ATres-3rd-choice.txt -n 24 

python3 $PYTHONDIR/CANVAS-MPI4.py -g $inputDIR/372_probs_TEMP.gro -l list-atoms-opt-3.txt -t $inputDIR/topol.top -n 24 

rm -rf test-antitrypsin-amber99-gromacs
mkdir test-antitrypsin-amber99-gromacs

mv run_simulation/ other-files/ analysis/ test-antitrypsin-amber99-gromacs
