#!/bin/bash

# Executing serial version of "block.py" and then "CANVAS.py" for the APO-4 conformation of Pembrolizumab Antibody 
# NOTE: The output of block.py is the list of survived atoms, i.e. list-atoms-opt-2.txt, that will be the input of CANVAS.py script  

PYTHONDIR=../PYTHON-scripts

inputDIR=../input-files/PEMBROLIZUMAB-ANTIBODY/sim-1A-amber99

python3 $PYTHONDIR/block.py choice2 -g $inputDIR/apo-1A.gro -l $inputDIR/list-ATres-2nd-choice.txt 

python3 $PYTHONDIR/CANVAS.py -g $inputDIR/apo-1A.gro -l list-atoms-opt-2.txt -t $inputDIR/topol.top -d $inputDIR/dom.txt

rm -rf test-apo-1A-amber99-gromacs
mkdir test-apo-1A-amber99-gromacs

mv run_simulation/ other-files/ analysis/ test-apo-1A-amber99-gromacs 
