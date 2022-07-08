<div align="center">

<img src="images/canvas-logo.jpg" alt="Scheme" width="950">
</div>

# 0 - Introduction 

The **CANVAS** model is a novel Multiple-Resolution approach which allows one to model at an atomistic resolution only the precise subset
of degrees really  necesssary for the study of a given phenomenon, even when this leads to a boundary between resolutions which falls 
within a bio-molecule. 

CANVAS is the acronym of **C**oarse-grained **A**nisotropic **N**etwork model for **VA**riable resolution **S**imulation. 

This model enables one to set the level of resolution of the coarse-grained subdomain(s) in a quasi-continuous range, spanning 
from the all-atom level to a degree o coaresening higher than one bead per amino acid. We emphasize that this is the novelty of the method: 
in paricular, one has an extremely high freedom in the choice of the level of coarse-graining. Moreover, the parametrization of interactions 
on the basis of the model chosen has an automatic construction: they do not require reference simulations, but rather than only the all-atom stucture of bio-molecule. 

The system under examination is the present work for validating the CANVAS model is the **pembrolizumab** antibody and **Adenylate Kinase** 
protein in acqueous solution.   

All the details of the work can be found in [REF: Fiorentini, Tarenzi, Potestio]

<br />

# 1 - Requirements  

The only requirements are the following: 

* **`GROMACS`**: it is a versatile package to perform molecular dynamics. The installation guide is reported 
  on [this link](https://manual.gromacs.org/documentation/2018-current/install-guide/index.html). 
  Please, note that we tested our code with GROMACS-2018, but nothing should change using another (older or newer) 
  version of this simulating package.

* **`Python3`**: it is an interpreted, object-oriented, high-level programming language with dynamic semantics. 
  The installation guide is provided [Here](https://docs.python-guide.org/starting/installation/). 
  If you are working on _Linux_ or _MacOs_ system, Python3 should be already installed. 
  On the other hand, if you are using Windows operating system, it is not certain for its presence.
  Please, be care of working with Python3 (3.7 or 3.9 is the best choice) as the code would return an error if using Python2.
    
* **`VMD`**: it is a molecular visualization program for displaying, animating, and analyzing large biomolecular systems
  using 3-D graphics and built-in scripting. 
  [Here](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD) it is possible to find the last 
  version of VMD, while the installation guide can be found on [this link](https://www.ks.uiuc.edu/Research/vmd/current/ig/node6.html). The biomolecule visualization and some analyses are based on its use, therefore **we strongly 
  recommend its installation**, even tough it is not mandatory. 

* **`LAMMPS`**:  it is a versatile package to perform molecular dynamics. The installation guide is reported 
  on [this link](https://docs.lammps.org/Install.html). Its installation is mandatory only in case of canvas 
  simulation of BioMolecule in Lammps (instead of GROMACS).  
<br />

# 2 - Tree Diagram 

Before looking the usage of this code, a tree diagram, representing the sequence of files and directories, is present for easier global reading.

```mermaid
flowchart LR
    id1((CANVAS))-->id2[/PYTHON-scripts/]
    style id1 fill:#f96
    style id2 fill:#9accdc
    style id3 fill:#9accdc
    style id4 fill:#9accdc
    style id5 fill:#9accdc
    style id6 fill:#9accdc
    style idT1 fill:#9accdc
    style idT2 fill:#9accdc
    style idT3 fill:#9accdc
    style id7 fill:#fef5cc
    style id8 fill:#fef5cc
    style id9 fill:#fef5cc

    id1-->id3[/images/];
    id1-->id4[/input-files/];
    id1-->id5[/lib/]; 
    id1-->id6[/tests/];
    id1-->id7([README.md]);

    id2-->id8([CANVAS.py]); 
    id2-->id9([block.py])

    id3-.->idI[images for README.md]
    
    id4-->idT1[/ADENYLATE-KINASE/]
    id4-->idT2[/ANTITRYPSIN/]
    id4-->idT3[/PEMBRLIZUMAB-ANTIBODY/]

    id5-.->idA[python libraries] 

    id6-.->idTest[.sh scripts for executing CANVAS.py and block.py]  
```

<br />

Directory folders are shown in light blue; files are shown in light yellow, while grey boxes explain only what is inside the folders. In this latter case, files and/or directory  are not indicated in this diagram for sake of clarity and compactness.  

# 3 - Usage 

The typical usage of the program consists in a call to _`block.py`_ and _`CANVAS.py`_ in succession 
by using Python3. Afterwards, it is possible to simulate the BioMolecule through Gromacs or Lammps, as proposed in **Sec. 7**. 

* **`block.py`**: has the purpose to write a file containing the list of survived atoms, that must be used 
                in _CANVAS.py_ as mandatory argument as explained in **Sec. 4.1**. 
                
* **`CANVAS.py`** has the purpose to write the input files needed for simulating a solvated system in 
                Multiple Resolution in GROMACS or LAMMPS and analyzing it. 

                

Before running the python scripts, read carefully the next section that provides a detailed explaination of each task 
and argument. Moreover, take care to not moving them outside the main folder (`canvas/`) otherwise a fatal error is 
printed on screen.

<br />

# 4 - _block.py_

## 4.1 - Tasks 

_`block.py`_ is inside the `PYTHON-SCRIPT/` directory and it allows the user to select one of three possible options, depending on the type of 
_`atomistic/medium-grained/coarse grained`_ subdivision that would like to obtain: 

* **{-choice1-}**: One or more central atomistic residues that require an atomistic description is/are known. 
               Since the high-resolution region is not completely defined, an atomistic sphere with radius _R_, 
               defined by the user, is traced (around the central residue(s)). Then, the latter, is sourrounded by 
               a 3D-annulus of width _D_ that defines a transition/hybrid region where only the backbone atoms ($N$\, 
               $C_\alpha$, $C$, $O$) are retained. The remainder, is modelled coarse-grained, where only 
               the $C_\alpha$ atoms are kept. A schematic represention is shown hereafter: 

<div align="center">

<img src="images/choice1.jpg" alt="Scheme" width="385">
</div>

<br /><br />

* **{-choice2-}**: All residues that require an atomistic description are known; therefore the higher resolution region
                   is completely defined. Around, an hybrid region of width _D_ will be traced where only the backbone 
                   atoms ($N$, $C_\alpha$, $C$, $O$) are retained. Note that in *choice1* the atomistic region 
                   is not completely defined a priori, but it will be described  by a radius _R_ starting from the 
                   knowledge of one or more central residues.  

<div align="center">

<img src="images/choice2.jpg" alt="Scheme" width="400">
</div>

<br /><br />

* **{-choice3-}**: All residues that require an atomistic and hybrid description(where only the backbone atoms 
                   are retained i.e. $N$, $C_\alpha$, $C$, $O$) are known. Consequently, the residues that 
                   will be modelled as corse grained (where only $`C_\alpha`$ atoms are kept) are also automatically 
                   enstablished. 

<div align="center">

<img src="images/choice3.jpg" alt="Scheme" width="400">
</div>

<br /><br />

Each task can require different input files, provided to the program in the form of command-line options. A short explaination of tasks and arguments is given by launching the command `python3 block.py -h` or `python3 block.py --help`. Alternatively, for printing a short usage message, please type: `python3 block.py` or `python3 block.py -u`

<br />

### 4.1.1 - choice1 
---------

**{-Choice1-}** option requires two mandatory files, i.e. the _{+coordinate FILE+}_, and the _{+list AT-bb-CG FILE+}_ 
and one optional argument that is the _{+diameter of hybrid region+}_ in the CANVAS model. The arguments are described in **Sec. 4.2**

In order to launch the **choice1** task the command-line is the following:

```bash
python3 block.py choice1 -g <Coordinate FILE> -l <list AT-bb-CG FILE> [-D <diameter hybrid region>] 

or 

python3 block.py choice1 --gro  <Coordinate FILE> --list <list AT-bb-CG FILE> [--diameter <diameter hybrid region>] 
```

The output of the program is the list of survived atoms. For further information, please type on terminal `python3 choice1`.

<br />

### 4.1.2 - choice2 
---------

**{-Choice2-}** option requires two mandatory files, i.e. the _{+coordinate FILE+}_, and the _{+list AT-bb-CG FILE+}_ 
and one optional argument that is the _{+diameter of hybrid region+}_ in the CANVAS model. The arguments are described in **Sec. 4.2**. 

In order to launch the **choice2** task the command-line is the following:

```bash
python3 block.py choice2 -g <Coordinate FILE> -l <list AT-bb-CG FILE> [-D <diameter hybrid region>] 

or 

python3 block.py choice2 --gro  <Coordinate FILE> --list <list AT-bb-CG FILE> [--diameter <diameter hybrid region>] 
```

The output of the program is the list of survived atoms. For further information, please type on terminal `python3 choice2`.

<br />

### 4.1.3 - choice3 
---------

**{-Choice3-}** option requires two mandatory files, i.e. the _{+coordinate FILE+}_, and the _{+list AT-bb-CG FILE+}_. The arguments are described in **Sec. 4.2**

In order to lunch the **choice3** task the command-line is the following:

```bash
python3 block.py choice2 -g <Coordinate FILE> -l <list AT-bb-CG FILE> 

or 

python3 block.py choice2 --gro  <Coordinate FILE> --list <list AT-bb-CG FILE> 
```

The output of the program is the list of survived atoms. For further information, please type on terminal `python3 choice3`.

<br />

## 4.2 - Arguments 

As shown in Section 4.1 the _coordinate FILE_ and list _AT-bb-CG FILE_ are always mandatory whatever the task chosen. 
On the other hand, if the option taken is **choice1** or **choice2**, the user has to opportunity to change the default value 
of diameter or the hydrid region (1.0 nm is the default one). 

Please note that the _list AT-bb-CG FILE_ is organized in different way according the task chosen.

A short explaination of files and of the diameter value is the following:

* **{-Coordinate FILE-}**: File of atom Coordinates in .gro format 

* **{-List At-bb-CG FILE (for choice1)-}**: File with a list of central atomistic(s) residue(s)
                                    and corresponding atomistic(s) radius(ii) organized in two columns:
```                                                           
|-----------------|-----------------------|                                        
| residue1 (int)  |   radius1 (int/float) |  
| residue2        |   radius2             | 
| ........        |   ........            |
| residueN        |   radiusN             |  
|-----------------|-----------------------|                                     
```                                    
                                    
* **{-List AT-bb-CG FILE (for choice2)-}**: File organized in only one column that contains 
                                    the list of atomistic residues:
```
|----------------|  
| residue1 (int) |  
| residue2       |  
| ........       |  
| residueM       |  
|----------------| 
```                                        
                                        
* **{-List At-bb-CG FILE (for choice3)-}**: File organized in two columns that contains the list of atomistic residues
                                    and the list of residues that require an hybrid resolution:
```
|-----------------|-----------------|  
| res1_AT (int)   |  res1_hy (int)  |  
| res2_AT         |  res2_hy        |  
| res3_AT         |  ........       |  
| res4_AT         |  resM_hy        |  
| ........        |                 |  
| resN_AT         |                 |  
|-----------------|-----------------|  
```

* **{-Diameter hybrid region-}**: Value (in nm) of diameter of the hybrid region. Default value: 1.0 nm
                                         

In **Appendix**  we focus on each argument discussed breafly before.

<br />

# 5 - _CANVAS.py_

_CANVAS.py_ does not have _tasks_ in the sense explained in **Sec. 4.1**. Indeed, this script
requires three mandatory files: _{+coordinate FILE+}_, _{+list survived atoms FILE+}_, and the _{+topology FILE+}_, 
and three optional argument, namely the _{+Domains FILE+}_, the flag _{+-c/--code lammps+}_ and the flag _{+-r/--resc14 Y+}_ . The arguments are described in **Sec. 5.1**

In order to launch the **CANVAS.py** script, the command-line is the following:

```bash 
python3 CANVAS.py [-h] -i <Coordinate FILE> -l <List survived atoms FILE> -t <Topology FILE> [-d <Domains FILE>] [-r Y] [-c lammps]

   or: 
   
python3 CANVAS.py [--help] --in <Coordinate FILE> --list <List survived atoms FILE> --top <Topology FILE> [--dom <Domains FILE>] [--resc14 Y] [--code lammps]
```
Please, take in account that _list survived atoms FILE_ is the output file obtained after launching _block.py_. 

A short explaination of arguments is provided by launching the command `python3 CANVAS.py -h` or `python3 CANVAS.py --help`. Alternatively, for printing a short usage message, please type: `python3 CANVAS.py` or `python3 CANVAS.py -u`

The output of the program consists of three directories, detailled described in **Sec. 7**, with the purpose of simulating with Gromacs or Lammps a Biomolecule with the CANVAS model, and analyze the resulting data: 
 
* **other-files/**
* **run_simulation/**
* **analysis/**

<br />

## 5.1 - Arguments 

As shown in **Sec. 5** the coordinate file, the file containing the list of survived atoms, and the topology one are always, mandatory. On the other hand, the domain division file is strongly recommended, even though is optional, the rescaled non-bonded 1-4 interactions flag 
(`-r/--resc14 Y`) and the string that indicates if LAMMPS or GROMACS input-files are created (`-c/--code lammps`) are optional. In the latter case, if not specified, gromacs input-files are returned.  

A short explaination of the above mentioned files is the following:

* **{-Coordinate FILE-}**: File of atom Coordinates in .gro format

* **{-List survived atoms FILE-}**: File containing the list of survived atoms organized in two columns:
                                It is the output file of block.py if correctly executed. 
```                                
|----------------|--------------------------| 
| at_num 1 (int) | caption 1 ('at' or 'cg') |
| at_num 2       | caption 2                |
| at_num 3       | caption 3                |
| at_num 4       | caption 4                |
| ....           | ....                     |
| at_num N       | caption N                |
|----------------|--------------------------|  
```

* **{-Topology FILE-}**: File containing the topology of fullyAT representation

* **{-Domain Division FILE-}**: Optional - but strongly recommended - File with the list of atoms divided in domains.
                            The number of columns is equal to the number of domains
                            Each row contains the at_numbers of the atoms that belong
                            to the same block separated by spaces:
```
|---------------------------------------------------|
| AT_Num-1   AT_Num-2   AT_Num-3   .....  AT_Num-N  |
| AT_Num-10  AT_Num-11  AT_Num-12  .....  AT_Num-M  |
| .....      .....      .....      .....  .....     |
|---------------------------------------------------| 
```
* **{-rescaled14-}**: String containing the word `Y` or `y`. Other strings
                  are not allowed. If `-r/--resc14 Y` is set, than all pairs
                  where ONLY one CG bead is kept if in the corresponding
                  all-atom representation that pair is present (by default 
                  only fully-atomistic rescaled non-bonded 1-4 interaction are mantained).
                  Keep attention, as it might create artifacts in MD simulation. 
                  If using charmm.ff this command is ignored. 

* **{-codestring-}**: String containing `lammps` word. Other strings are not
                      allowed. If `-c/--code lammps` is set, then this program
                      produces the input files needed for simulating the CANVAS
                      model in lammps. If `-c/--code gromacs` is set, or in case
                      this flag is ignored, this program returns the input files
                      for simulating the CANVAS model in gromacs.


In the **Appendix** we focus on each argument discussed breafly before.

<br />

# 6 - Examples 

Inside the `tests/` directory there is the complete list of example files for the Pembrolizumab, Adelynate Kinase, and 
Antytripsin biomolecules, allowing the user to try **block.py** and **CANVAS.py**. 

Hereafter, for the sake of clarity, only three examples are reported. 

 
```bash 
PYTHONDIR=../PYTHON-scripts

inputDIR=../input-files/ANTITRYPSIN

python3 $PYTHONDIR/block.py choice3 -g $inputDIR/372_probs_TEMP.gro -l $inputDIR/list-ATres-3rd-choice.txt

python3 $PYTHONDIR/CANVAS.py -g $inputDIR/372_probs_TEMP.gro -l list-atoms-opt-3.txt -t $inputDIR/topol.top
```

```bash
PYTHONDIR=../PYTHON-scripts

inputDIR=../input-files/ADENYLATE-KINASE/kinase-w-charmm36m

python3 $PYTHONDIR/block.py choice2 -g $inputDIR/frame_4ake_125ns.gro -l $inputDIR/list-ATres-2nd-choice.txt

python3 $PYTHONDIR/CANVAS.py -g $inputDIR/frame_4ake_125ns.gro -l list-atoms-opt-2.txt -t $inputDIR/topol.top -d $inputDIR/dom.txt
```

```bash 
PYTHONDIR=../PYTHON-scripts

inputDIR=../input-files/PEMBROLIZUMAB-ANTIBODY/APO-0

python3 $PYTHONDIR/block.py choice2 -g $inputDIR/apo0.gro -l $inputDIR/list-ATres-2nd-choice.txt

python3 $PYTHONDIR/CANVAS.py -g $inputDIR/apo0.gro -l list-atoms-opt-2.txt -t $inputDIR/topol.top -d $inputDIR/dom.txt
```

The output files of each test can be also found in `canvas/output-files/` directory.  

<br />

# 7 - Simulating CANVAS model with GROMACS, analyze data, and Visualizing Proteins   

After executing in sequence _block.py_ and _CANVAS.py_, three subfolders are created:

1. **{-other-files-}**: it contains a lot of intermediate files created by CANVAS.py: they are not useful per running simulation and for analyzing data. Thus, they can be ignored. 

2. **{-run-simulation-}**: it contains all the ingredients needed for simulating a BioMolecule in Multiple Resolution through the CANVAS model. 

3. **{-analysis-}**: it contains .tcl script and text files for analyzing and visualizing the biomolecule in VMD

<br />

## 7.1 - Launching the Simulation in Gromacs 

In  order to launch the simulation of a biomolecule trhough the CANVAS model, go inside the **{+run-simulation/+}** directory and follow the standard simulation protocol: minimization, equilibration in nvt (50 ps), equilibration in npt (50 ps), and finally the run production (500 ns). 

```bash 
gmx_mpi grompp -f em.mdp -c solvated_ions.gro -p topol_new.top -o em.tpr
gmx_mpi mdrun -v -deffnm em
```
```bash 
gmx_mpi grompp -f nvt.mdp -c em.gro -r em.gro -p topol_new.top -o nvt.tpr
gmx_mpi mdrun -v -deffnm nvt
```
```bash 
gmx_mpi grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol_new.top -o npt.tpr
gmx_mpi mdrun -v -deffnm npt
```
```bash 
gmx_mpi grompp -f run.mdp -c npt.gro -t npt.cpt -p topol_new.top -o run.tpr
gmx_mpi mdrun -v -deffnm run
```

Please, take in account that when using **charmm36-jul2021.ff**, a fatal error when launching the minimization
is printed on terminal. The latter is due to the fact that sodium and chloride ions in our files are labelled with **NA** and **CL** respectively, as for any other forcefield. Indeed, this version requires that the latter must be marked with **SOD** and **CLA**. 

Therefore, until this issue is solved in the next version of Charmm36m.ff, it is possible to fix the error, substituting manually **NA** and/or **CL** ions present in `solvated_ions.gro` and in `topol_new.itp` files, with **SOD** and **CLA** respectively. 

After this change, the standard simulating protocol can be followed without errors. 

<br />

## 7.2 - Analysing Data & Visualizing Trajectory 

VMD installation is not mandatory; however it is useful for visualing trajectory and for making preliminary analysis, such as the Root Mean Square Fluctuation (RMSF) of only $`C_\alpha`$ carbon (atomistic, medium-grained and 
coarse-grained). This folder contains the following files: 

* **{-CA_list.txt-}**: list of the $`C_\alpha`$ atoms indexes separated by one white space, in one line. 

* **{-radius_charge-allatom.txt-}**: this file is divided in four columns containing the index, the atom-name, 
  the radius (in Angstrom units) obtained by the value of $`\sigma/2`$, and the corresponding charge. This file 
  is then read by the tkl script, described hereafter. 
```        
|-------|----------------|-----------------|-----------------|  
| index |    at-type     | VdW_radius (AA) |    charge       |                            
|-------|----------------|-----------------|-----------------| 
|   1   | at-type1 (str) | radius1 (float) | charge1 (float) |
|   2   | at-type2       | radius2         | charge2         |    
|   3   | at-type3       | radius3         | charge3         |
|   4   | at-type4       | radius4         | charge4         |
|  .... |  ......        | .......         | ......          |
|   N   | at-typeN       | radiusN         | chargeN         |
|-------|----------------|-----------------|-----------------|  
```

* **{-radius_charge-allatom.tcl-}**: this script has the purpose of reading the _radius_charge-allatom.txt_ file. 
  In order to execute the script by means of VMD, please launch the command `source radius_charge-allatom.tcl` in
  TkConsole. Then, in order to visualize the radius size of each CG bead with a different coloration given by 
  the charge value, go on `Graphics/Representation/Create Rep` and create the following representation:

  ```bash 
  Select Atoms = resname MUL, Drawing method = VDW, Coloring Method = beta
  ```
  In this way, each bead will have a different radius according the value of $`\sigma/2`$, and a color given by the 
  Beta Coloring Method: in particular, blue shades are indicative of positive charges, red shades are indicative of
  negative charges, while the white color corresponds to a neutral charge. The figure below displays the apo-0 (4A) form 
  of Pembrolizumab in terms of this new representation just mentioned.  

<div align="center">
 <img src="images/apo0-4A.jpg" alt="Scheme" width="350">
</div>

<br />

* **rmsf.tcl**: this script can be executed on the VMD TkConsole by launching in sequence the following commands:

  ```bash
  source rmsf.tcl
  rmsf_all_ca 
  ```
  If everything went fine, `rmsf.dat` file is returned in output. The latter contains the fluctuation values of each 
  $`C_\alpha`$ atom (atomistic, medium-grained, and coarse-grained)

<br />

# 8 - Fixing Errors when using Charmm36m forcefield 

This section reports two errors that it is possible to encounter if using charmm36m forcefields, and how to solve them. The first one regards the creation of all-atom topology given the coordinate file(.pdb or .gro). The second problem concernes the addition of ions for neutralizing the charge.  Thus, it is clear that the latter are forcefield dependent, and they are not related to the CANVAS model.

<br />

## 8.1 - Creating topology file 

When the coordinate (.gro) and the topology file (.top) is created using the Charmm36m forcefield using the command provided in Appendix-A below, a fatal error like this can be encountered: 

``` 
Fatal error:
Residue 194 named GLY of a molecule in the input file was mapped
to an entry in the topology database, but the atom CB used in
that entry is not found in the input file. Perhaps your atom
and/or residue naming needs to be fixed.

For more information and tips for troubleshooting, please check the GROMACS
website at http://www.gromacs.org/Documentation/Errors
-------------------------------------------------------
```

As proposed in this [link](https://gromacs.bioexcel.eu/t/newest-charmm36-port-for-gromacs/868/8), the problem is that the COO- terminus entry in `aminoacids.c.tdb` uses atom _CB_ as the third control atom for adding _OT1_ and _OT2_. Until it is fixed in the next release of the forcefield, please edit the line `2 8 OT C CA CB` in **aminoacids.c.tdb** to `2 8 OT C CA N` that fix the issue. 

Note that the charmm36m.ff present in `input-files` directory is modified according what above explained.

<br />

## 8.2 - Simulating CANVAS Model 

As already seen in **Sec. 7.1**, another issue can be encountered when using the last version (july 2021) of Charmm36m is due to a different denomination of ions. In every force-field sodium and chloride ions are indicated with NA and CL, respectively. However, in charmm36m 
they are named SOD and CLA. A solution for fixing the issue is modifying manually the `topol_new.itp` and `solvated_ions.gro` files, substituing: 
  *  **{-NA-}** with **{+SOD+}**  
  *  **{-CL-}** with **{+CLA+}** 

<br />

# 9 - Contacts 

Raffaele Fiorentini: raffaele.fiorentini@unitn.it


<br /><br />


# Appendix 

Hereafter, we focus on the arguments discussed breafly before in **Sec. 4.2**  and **Sec. 5.1**. 

## A - Coordinate input FILE 

The coordinate file is mandatory for both _block.py_ and _CANVAS.py_.
The same file is recommended when launching these two python scripts.  

It must be provided in _gro_ format. It contains a molecular structure, that is the coordinates 
of each atom in the reference structure. As an example a the coordinate gro file of a biomolecule named Adenylate 
Kinase with 3241 atoms is the following: 


```latex
Adenylate Kinase 
 3341
    1MET      N    1   4.508   3.059   5.730
    1MET     H1    2   4.415   3.098   5.714
    1MET     H2    3   4.488   2.964   5.756
    1MET     H3    4   4.556   3.093   5.813
    1MET     CA    5   4.587   3.078   5.609
    ....
    ....  
  214GLY    HA2 3338   4.388   2.523   4.183
  214GLY      C 3339   4.445   2.373   4.325
  214GLY    OC1 3340   4.324   2.345   4.340
  214GLY    OC2 3341   4.543   2.319   4.386
   9.97074   9.97074   9.97074
```
> For simplicity we cut the file. The integral version can be found in `input-files/ADENYLATE-KINASE/kinase-w-amber99`
   
Detailed information can be found in [Gromacs-GroFile](https://manual.gromacs.org/archive/5.0.3/online/gro.html)

It is, also, possible to find coordinate files with pdb extension. Indeed, the crystallographic structure 
of a biomolecule is, in general, in the latter format. The most of proteins/biomolecules can be dowloaded 
from [_Protein Data Bank_](https://www.rcsb.org) in pdb format.

However, it is possible to transform pdb format in the gro one by means of a [pdb2gmx](https://manual.gromacs.org/documentation/current/onlinehelp/gmx-pdb2gmx.html#gmx-pdb2gmx)
GROMACS tool.

```bash
pdb2gmx -f <pdb FILE> -o <gro FILE> -water tip3p -ignh 
```

**Amber99sb-ildn** and **charmm36m** force-fields have been tested for performing a Multiple Resolution Simulation using CANVAS model; however, any other previous version of **Amber** and **Charmm** should not be a big deal. 

> Note: in order that a forcefield is found, it is necessary that the main folder containing the latter is inside the gromacs installing package, or in current directory where the pdb file is present. 

Please, take in account that the _coordinate FILE_ is the reference structure from which the CANVAS model is constructed; 
therefore choose it carefully. Indeed, you can select the crystallographic version of the biomolecule, or other equilibrated frames.  

At the end, if everything went fine, the coordinate file in .gro format, and the topology file (topol.top) will be created.
The latter will be discussed in **Appendix B**.

<br />

## B - Topology FILE 

The topology file is mandatory for _CANVAS.py_. It is not used when launching _block.py_. 

This file contains the parameters and the bonded interactions of the reference structure. 
The topology file usually specify _`[bonds]`_ (2 atoms connected), _`[angles]`_ (3 atoms connected), and _`[dihedrals]`_ 
(4 atoms connected linearly).

*Amber99sb-ildn* and *Charmm36m* forcefields has been successfully tested and employed. 

As example, we report in [Urea-top](https://manual.gromacs.org/documentation/current/reference-manual/file-formats.html#top) the topology file for **Urea** in water. The topology files of Adelynate Kinase and Pembrolizumab employed in our work can be found in **`input-files/`** folder.

Further and detailed information can be found in [Gromacs-Topology](https://manual.gromacs.org/documentation/current/reference-manual/topologies/topology-file-formats.html)

A top file can be generated by [pdb2gmx](https://manual.gromacs.org/documentation/current/onlinehelp/gmx-pdb2gmx.html#gmx-pdb2gmx) after choosing a pdb or gro file: 

```bash
pdb2gmx -f <pdb/gro FILE> -water tip3p
```

At the end, if everything went fine, the topology file (topol.top) will be created.

<br />

## C - List AT-bb-CG FILE 

_list AT-bb-CG_ is a mandatory file when launching _block.py_. As already explained in **Sec. 4.2**, this file  
is organized in different way according the task chosen: 

* **{-List At-bb-CG FILE (for choice1)-}**: File with a list of central atomistic(s) residue(s)
                                        and corresponding atomistic(s) radius(ii) organized in two columns. 
                                        The left column must present only integer numbers; the right one 
                                        must show only integer or float number. No strings, or special characters 
                                        are permitted. An error message is printed on screen if some condition is not fulfilled. 
```                                                                       
|-----------------|-----------------------|                                        
| residue1 (int)  |   radius1 (int/float) |  
| residue2        |   radius2             | 
| ........        |   ........            |
| residueN        |   radiusN             |  
|-----------------|-----------------------|                                     
```                                   
                                    
* **{-List AT-bb-CG FILE (for choice2)-}**: File organized in only one column that contains 
                                        the list of atomistic residues. Only integer numbers are, 
                                        thus, permitted. No strings, or special characters are feasible. 
                                        An error message is printed on screen if some condition is not fulfilled. 
```
|----------------|  
| residue1 (int) |  
| residue2       |  
| ........       |  
| residueM       |  
|----------------| 
```                                        
                                        
* **{-List At-bb-CG FILE (for choice3)-}**: File organized in two columns that contains the list of atomistic residues
                                        and the list of residues that require an hybrid resolution. Thus, both columns 
                                        must present only integer numbers. No strings, or special characters are 
                                        feasible. 
                                        An error message is printed on screen if some condition is not fulfilled. 
                                        
```
|-----------------|-----------------|  
| res1_AT (int)   |  res1_hy (int)  |  
| res2_AT         |  res2_hy        |  
| res3_AT         |  ........       |  
| res4_AT         |  resM_hy        |  
| ........        |                 |  
| resN_AT         |                 |  
|-----------------|-----------------|  
``` 

<br />

## D - List survived atoms FILE

The file containing the list of survived atoms is a mandatory argument when launching _CANVAS.py_.
The user is not require to write this file, as the latter is the output of _block.py_.

Hereafter, we report an example of list: 
   
```latex
   ...   ..
   449   at
   450   at
   451   at
   452   cg
   454   cg
   461   cg
   475   cg
   ...   ..
```

This file is organized in two columns: 

* The {+1st column+} corresponds to the atom number of the atoms that survive from the fullyAT 
  reference 
  
* The {+2nd column+} corresponds to the label of each atom that survives: 
     * **`at`** stands for **atomistic**: it means that the the atom token in account conserves it own properties. 
  
     * **`cg`** stands for **coarse-grained or hybrid**: it means that the atom token in account will be treated 
       as CG bead and it will have average properties of the atoms (decimated) that it represents. This label should be 
       used for both the coarse-grained part (only $`C_\alpha`$ atoms are retained) and hybrid part (backbone atoms 
       are kept: $`N`$, $`C_\alpha`$, $`C`$, $`O`$)

For the sake of completeness, please note that it is possible to write by yourself a text file 
containing the list of the atoms, each one labelled with the caption (`at` or `cg`) if you know 
in advance all the atoms that survived.  However, we do not encurage this precedure to continue
especially in case of biomolecules with tens of thousand atoms, since the probability to make 
mistakes becomes higher and higher.

<br />

## E - Domains FILE 

The _Domains FILE_ is an optional argument when executing _CANVAS.py_, but we strongly encurage to make use of it. 
It contains the list of atoms divided in different domains. 

Given the fully atomistic representation, we _transform_ the latter in multiple resolution ready to be simulated in GROMACS: according with the standard procedure, a bond is placed between atoms where at least a CG bead is involved if their distance lies below a cutoff. However, sometimes we would like to avoid bonds between atoms that belong 
to different domains. If we do not, the CANVAS model might result very rigid. 

For each domain, we report in a row the list of all the atoms that belong to it as reported hereafter: 
 
```latex 
1 2 3 4 5 6 7 8 9 10 11 12 13
14 15 16 17 18 19 20 21 22 23
24 25 26 27 28 29 30 31 32 33
```

Therefore: 

* The number of columns of this file is equals to the number of domains. 

* Each number corresponds to the atomistic number (atnum) of the atom in all-atom representation.  

* Each row contains all the atoms that belong to a certain domain, each one separated by spaces. 

There is no automatic procedure to obtain automatically this file. Thus, you must write your 
own txt file. The following script could be useful to construct the file _dom.txt_ for N domains: 

```bash 
#!/bin/bash

rm dom.txt
touch dom.txt

dom1=$(for i in {A1..An}; do echo -n $i ""; done)
echo $dom1 >> dom.txt

dom2=$(for i in {B1..Bn}; do echo -n $i ""; done)
echo $dom2 >> dom.txt

domX=$(for i in {M1..Mn}; do echo -n $i ""; done)
echo $domX >> dom.txt
```

where in place of $`A_1, A_2, B_1, B_2, M_1, M_2`$ you should write the atomistic numbers. Please, note  that the script works good if each domain 
is constituted by consecutive at_numbers.

<br />

## F - Diameter hybrid region 

The value of the diameter of the hybrid region is set to the default value _1.0 nm_. If you want to change it use the  
flag `-D <value (nm)>` when launching _block.py_. Integer or decimal numbers are equally accepted. be sure that _`D`_ 
is higher or equals to 0. _`D = 0`_, means that no hybrid region is constructed and the system assumes only two resolutions: **fullyAT** and **Coarse-grained** (keeping only C-alpha atoms).

As already explained in **Sec. 4.2**, the value of diameter for hybrid region can be set only when selecting **choice1** and **choice2** options. Indeed, **choice3** does not require the setting of its value since is unnecessary. 

<br />

## G - Rescaled14

String containing the word `Y` or `y`. Other strings are not allowed. If `-r/--resc14 Y` is set, than all pairs where ONLY one CG bead is kept if in the corresponding all-atom representation that pair is present (by default only fully-atomistic rescaled non-bonded 1-4 interaction are mantained). Keep attention, as it might create artifacts in MD simulation. If using charmm.ff this command is ignored. 

<br />

## H - Codestring 

String containing `lammps` word. Other strings are not allowed. If `-c/--code lammps` is set, then this program produces the input files needed for simulating the CANVAS model in lammps. If `-c/--code gromacs` is set, or in case this flag is ignored, this program returns the input files
for simulating the CANVAS model in gromacs.


