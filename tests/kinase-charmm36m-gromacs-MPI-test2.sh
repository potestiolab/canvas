#!/bin/bash

# Executing MPI version of "block.py" and then "CANVAS.py" for the APO-0 conformation of Pembrolizumab Antibody 
# NOTE: The output of block.py is the list of survived atoms, i.e. list-atoms-opt-2.txt, that will be the input of CANVAS.py script  

# In case you would like to execute the code in serial, please cancel the previuos code lines, and change block-MPI.py and CANVAS-MPI4.py
# with 'block.py' and 'CANVAS.py', respectively 

CHARMM_PATH="../input-files/ADENYLATE-KINASE/sim-w-charmm36m/charmm36-jul2021.ff"
GMX_PATH="/opt/homebrew/Cellar/gromacs/2023.1/bin/gmx"
RES_PATH="/opt/homebrew/share/gromacs/top/residuetypes.dat"

PYTHONDIR="../PYTHON-scripts"
inputDIR="../input-files/ADENYLATE-KINASE/sim-w-charmm36m"


python3 $PYTHONDIR/block-MPI.py choice2 -g $inputDIR/frame_4ake_125ns.gro -l $inputDIR/list-ATres-2nd-choice.txt  

echo -e "$CHARMM_PATH\n$GMX_PATH\n$RES_PATH" | python3 $PYTHONDIR/CANVAS.py -g $inputDIR/frame_4ake_125ns.gro -l list-atoms-opt-2.txt -t $inputDIR/topol.top -d $inputDIR/dom.txt -c lammps 


#python3 $PYTHONDIR/CANVAS.py -g $inputDIR/frame_4ake_125ns.gro -l list-atoms-opt-2.txt -t $inputDIR/topol.top -d $inputDIR/dom.txt 

rm -rf test-kinase-charmm36m-gromacs-nosolv
mkdir test-kinase-charmm36m-gromacs-nosolv

mv run_simulation/ other-files/ analysis/ test-kinase-charmm36m-gromacs-nosolv
