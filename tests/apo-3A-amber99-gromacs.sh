#!/bin/bash

# Executing serail version of "block.py" and then "CANVAS.py" for the APO-5 conformation of Pembrolizumab Antibody 
# NOTE: The output of block.py is the list of survived atoms, i.e. list-atoms-opt-2.txt, that will be the input of CANVAS.py script  

PYTHONDIR=../PYTHON-scripts

inputDIR=../input-files/PEMBROLIZUMAB-ANTIBODY/sim-3A-amber99

python3 $PYTHONDIR/block.py choice2 -g $inputDIR/apo-3A.gro -l $inputDIR/list-ATres-2nd-choice.txt 

python3 $PYTHONDIR/CANVAS-MPI4.py -g $inputDIR/apo-3A.gro -l list-atoms-opt-2.txt -t $inputDIR/topol.top -d $inputDIR/dom.txt

rm -rf test-apo-3A-amber99-gromacs
mkdir test-apo-3A-amber99-gromacs

mv run_simulation/ other-files/ analysis/ test-apo-3A-amber99-gromacs 
