#!/bin/bash

# Executing "block.py" and then "CANVAS.py" for Adenylate Kinase 
# NOTE: The output of block.py is the list of survived atoms, i.e. list-atoms-opt-2.txt, that will be the input of CANVAS.py script  

PYTHONDIR=../PYTHON-scripts

inputDIR=../input-files/ADENYLATE-KINASE/sim-w-amber99

python3 $PYTHONDIR/block.py choice2 -g $inputDIR/frame_4ake_125ns.gro -l $inputDIR/list-ATres-2nd-choice.txt  

python3 $PYTHONDIR/CANVAS.py -g $inputDIR/frame_4ake_125ns.gro -l list-atoms-opt-2.txt -t $inputDIR/topol.top -d $inputDIR/dom.txt -c lammps

rm -rf test-kinase-amber99-lammps
mkdir test-kinase-amber99-lammps

mv run_simulation/ other-files/ analysis/ test-kinase-amber99-lammps 
