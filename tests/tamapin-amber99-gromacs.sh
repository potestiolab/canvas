#!/bin/bash

# Executing "block.py" and then "CANVAS.py" for Antitrypsin 
# NOTE: The output of block.py is the list of survived atoms, i.e. list-atoms-opt-3.txt, that will be the input of CANVAS.py script  


PYTHONDIR=../PYTHON-scripts

inputDIR=../input-files/TAMAPIN/sim-w-amber99/

python3 $PYTHONDIR/block.py choice2 -g $inputDIR/6d93.gro -l $inputDIR/list-ATres-2nd-choice.txt 

python3 $PYTHONDIR/CANVAS.py -g $inputDIR/6d93.gro -l list-atoms-opt-2.txt -t $inputDIR/topol.top

rm -rf test-tamapin-amber99-gromacs
mkdir test-tamapin-amber99-gromacs

mv run_simulation/ other-files/ analysis/ test-tamapin-amber99-gromacs
