#!/bin/bash

#################################################################################################################################
### ATTENTION: Hereafter, an example for executing MPI scripts on 48 cores on UniTn Cluster. Change it for your cluster/machine 

#PBS -l select=1:ncpus=48:mpiprocs=1:mem=4gb
#PBS -q VARIAMOLS_cpuQ

echo "loading modules"
module load mpich-3.2

echo "modules loaded"

echo "running simulations"

cd $PBS_O_WORKDIR
#################################################################################################################################

# Executing MPI version of "block.py" and then "CANVAS.py" for the APO-0 conformation of Pembrolizumab Antibody 
# NOTE: The output of block.py is the list of survived atoms, i.e. list-atoms-opt-2.txt, that will be the input of CANVAS.py script  

# In case you would like to execute the code in serial, please cancel the previuos code lines, and change block-MPI.py and CANVAS-MPI4.py
# with 'block.py' and 'CANVAS.py', respectively

CHARMM_PATH=../input-files/PEMBROLIZUMAB-ANTIBODY/sim-0A-charmm36m/charmm36-jul2021.ff

PYTHONDIR=../PYTHON-scripts

inputDIR=../input-files/PEMBROLIZUMAB-ANTIBODY/sim-0A-charmm36m

python3 $PYTHONDIR/block-MPI.py choice2 -g $inputDIR/apo-0A.gro -l $inputDIR/list-ATres-2nd-choice.txt 

echo $CHARMM_PATH | python3 $PYTHONDIR/CANVAS-MPI4.py -g $inputDIR/apo-0A.gro -l list-atoms-opt-2.txt -t $inputDIR/topol.top -d $inputDIR/dom.txt -c lammps

rm -rf test-apo-0A-charmm36m-lammps
mkdir test-apo-0A-charmm36m-lammps

mv run_simulation/ other-files/ analysis/ test-apo-0A-charmm36m-lammps
