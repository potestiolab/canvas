#!/bin/bash

# Executing "block.py" and then "CANVAS.py" for Antitrypsin 
# NOTE: The output of block.py is the list of survived atoms, i.e. list-atoms-opt-3.txt, that will be the input of CANVAS.py script  


PYTHONDIR=../PYTHON-scripts

GMX_PATH="/opt/homebrew/Cellar/gromacs/2023.1/bin/gmx"
RES_PATH="/opt/homebrew/share/gromacs/top/residuetypes.dat"


inputDIR=../input-files/TAMAPIN/sim-w-charmm27/

python3 $PYTHONDIR/block-MPI.py choice2 -g $inputDIR/6d93.gro -l $inputDIR/list-ATres-2nd-choice.txt 

echo $GMX_PATH | python3 $PYTHONDIR/CANVAS-MPI4-prova.py -g $inputDIR/6d93.gro -l list-atoms-opt-2.txt -t $inputDIR/topol.top

rm -rf test-tamapin-charmm27-gromacs
mkdir test-tamapin-charmm27-gromacs

mv run_simulation/ other-files/ analysis/ test-tamapin-charmm27-gromacs
