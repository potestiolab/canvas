#!/bin/bash

### ATTENTION: Hereafter, an example for executing MPI scripts on 48 cores in UniTn Cluster. Change it for your cluster/machine 

#PBS -l select=1:ncpus=48:mpiprocs=1:mem=4gb
#PBS -q VARIAMOLS_cpuQ

module load mpich-3.2

cd $PBS_O_WORKDIR
######################################################################

# Executing "block.py" and then "CANVAS.py" for Antitrypsin 
# NOTE: The output of block.py is the list of survived atoms, i.e. list-atoms-opt-3.txt, that will be the input of CANVAS.py script  

FF_charmm_path=/home/raffaele.fiorentini/canvas-NEW-MPI/input-files/ANTITRYPSIN/sim-w-charmm36m/charmm36-jul2021.ff/

PYTHONDIR=../PYTHON-scripts

inputDIR=../input-files/ANTITRYPSIN/sim-w-charmm36m-FullyAT

echo $FF_charmm_path | python3 $PYTHONDIR/CANVAS-MPI4.py -g $inputDIR/Probs-charmm.gro -l $inputDIR/list-atoms.txt -t $inputDIR/topol.top -c lammps

rm -rf test-antitrypsin-charmm36m-lammps-FullyAT
mkdir test-antitrypsin-charmm36m-lammps-FullyAT

mv run_simulation/ other-files/ analysis/ test-antitrypsin-charmm36m-lammps-FullyAT
