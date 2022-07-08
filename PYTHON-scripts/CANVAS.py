# -*- coding: utf-8 -*-

import os
import sys
import math
from math import sqrt, exp, log 
from operator import itemgetter

import argparse

import itertools
import numpy as np 

import glob
from os.path import expanduser
#from mpmath import mpf

from datetime import datetime 


PYTHONPATH = os.path.abspath(os.getcwd())

spl_word = "canvas-NEW-MPI" # We find CANVAS, cut the entire path until CANVAS and add /lib in order to find our libraries. 

python_modules_path = PYTHONPATH.split(spl_word)[0] + spl_word + "/lib"

sys.path.append(python_modules_path)

from read_gromacs import * 
from check_gromacs_files import *
from write_gromacs_files import * 
from write_VMD_files import * 
from inp_out import *   
from write_lammps_files import * 

REAL_start = datetime.now() 

print("####################################################\n")
print("'{}' running...".format(os.path.basename(sys.argv[0]))) 
print("---------------------------------\n")

#n1. Input Arguments 

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, add_help=False )

group_in=parser.add_argument_group("Required Arguments")

group_in.add_argument('-g', '--gro',     dest='grofile', action='store', metavar = 'FILE', help = argparse.SUPPRESS)   # Mandatory 
group_in.add_argument('-l', '--list',    dest='listatomsfile', metavar = 'FILE',help = argparse.SUPPRESS)              # Mandatory
group_in.add_argument('-t', '--top',     dest='topfile',    metavar = 'FILE', help = argparse.SUPPRESS)                # Mandatory

group_in.add_argument('-h', '--help',    action='help', help = argparse.SUPPRESS)			               # Optional
group_in.add_argument('-d', '--dom',     dest='domfile',  metavar = 'FILE', help = argparse.SUPPRESS)                  # Optional 
group_in.add_argument('-r', '--resc14',  dest='rescaled14', metavar = 'STR',  help = argparse.SUPPRESS)                # Optional
group_in.add_argument('-c', '--code',    dest='codestring', metavar = 'STR', help = argparse.SUPPRESS) 		       # Optional 
group_in.add_argument('-s', '--solvate', dest='solvatestring', metavar = 'STR', help = argparse.SUPPRESS)              # Optional

#n2. Printing help message if the script does not have any arguments 
if len(sys.argv)==1:
    print_usage_main_CANVAS()
    sys.exit(1)

if(sys.argv[1].strip() == "--usage") or (sys.argv[1].strip() == "-u"):
    print_usage_main_CANVAS()
    quit()

if(sys.argv[1].strip() == "--help" or sys.argv[1].strip() == "-h"):
    print_help_main_CANVAS()
    quit()

#n3. Printing help message if the script does not present valid arguments
check_argv_errors_CANVAS()

#n4. Parsing arguments, printing error and help message if code presents not allowed arguments, and checking if mandatory files are present  
try:
    args  = parser.parse_args()
except SystemExit:
    print("\nArguments with no flag are not allowed. Check that each flag (e.g. -g) is followed by its specific argument")
    print("Look below for more information.\n")
    print_help_main_CANVAS()
    sys.exit()


grofile         = args.grofile
topfile         = args.topfile
listatomsfile   = args.listatomsfile
domfile         = args.domfile 
rescaled14      = args.rescaled14 
codestring      = args.codestring 
solvatestring   = args.solvatestring 

mandatory_files_present_CANVAS(grofile, listatomsfile, topfile)

#n5. Code Flag: LAMMPS or GROMACS. This program will return the input files for simulating the CANVS model in GROMACS (default) or LAMMPS. 
FlagLammps = False 

if(codestring is None):
    print("'-c/--code' is not set. This program will return the input files for simulating the CANVAS model in GROMACS\n")
    FlagLammps  = False 

if(codestring is not None):   # If codestring exists... 
    if(codestring == "lammps" or codestring == "Lammps" or codestring == "LAMMPS"): 
        print("'-c/--code {}' set. This script will return the input files for simulating the CANVAS model in LAMMPS... 1% completed\n".format(codestring)) 
        FlagLammps  = True
 
    elif(codestring == "gromacs" or codestring == "Gromacs" or codestring == "GROMACS"):
        print("'-c/--code {}' set. This script will return the input files for simulating the CANVAS model in GROMACS... 1% completed\n".format(codestring))
        FlagLammps  = False
    
    else:
        print("ERROR. '-c/--code {}' is set, but not recognized. Only 'lammps' or 'gromacs' strings are allowed.\n".format(codestring)) 
        print("o If you would like to return the input files for simulating the CANVAS model in GROMACS, please type:") 
        print("'-c/--code gromacs' or simply ignore this optional flag.\n")
        print("o If you would like to return the input files for simulating the CANVAS model in LAMMPS, please type:" ) 
        print("-'c/--code lammps'\n")
        print("Look below for further help\n")
        print_help_main_CANVAS()
        quit()

#n6. Solvate Flag. This program solvates or not the system for simulating the CANVAS model. 
if(solvatestring is None):
    if(codestring is not None):  
        print("'-s/--solvate' is not set. This program will solvate the system before simulating the CANVAS model in {}\n".format(codestring))
    else:
        print("'-s/--solvate' is not set. This program will solvate the system before simulating the CANVAS model in GROMACS\n")
    solvate_Flag = True 

if(solvatestring is not None):  
    if(solvatestring == "y" or solvatestring == "Y"):
        print("'-s/--solvate y' set. This script will solvate the system before simulating the CANVAS model in {}\n".format(codestring))
        solvate_Flag = True 

    elif(solvatestring == "n" or solvatestring == "N"): 
        print("'-s/--solvate n' set. This script will not solvate the system before simulating the CANVAS model in {}".format(codestring))
        print("This flag might be useful for implicit solvent simulations\n")
        solvate_Flag = False 

    else: 
        print("ERROR. '-s/--solvate {}' set, but not recognized. Only 'n/N' or y/Y' string are allowed.\n".format(solvatestring))
        print("o If you would like to solvate (and then neutralize) the system  before simulating the CANVAS model, please type:")
        print("'-s/--solvate y' or '-s/--solvate Y' or simply ignore this optional flag.\n")
        print("o If you would like to NOT solvate the system  before simulating the CANVAS model, please type:")
        print("'-s/--solvate n' or '-s/--solvate N'\n")
        print("Look below for further help\n")
        print_help_main_CANVAS()
        quit() 

                             #####################################
############################## (1) READING COORDINATE.gro FILE  ############################################################################
                             #####################################
''' 
    Reading the gromacs coordinate file .gro using the module "readgro_all" of "read_gromacs.py" library. The latter returns 5 outputs:

    # column[0] = list of lists containing: res_number, res_name, at_name, at_number, pos_x, pos_y, pos_x (and vel_x, vel_y and vel_z if they are present)
    # column[1] = number of total atoms of gro file 
    # column[2] = Box size, x-direction: Lx
    # column[3] = Box size, y-direction: Ly 
    # column[4] = Box size, z-direction: Lz
    # column[5] = Title of .gro file
'''

#1.1- Checking if grofile is present  
if not os.path.isfile(grofile):     
    print("Error while opening the file. '{}' does not exist\n".format(grofile))
    quit()


#1.2- Opening and reading grofile 
f = open(grofile, "r")
check_empty_file(grofile)

groo = readgro_all(f)

atom    = groo[0]
n_atoms = groo[1]
Lx      = groo[2]
Ly      = groo[3]
Lz      = groo[4]
title   = groo[5]

print("\nOriginal coordinate file correctly read... 2% completed\n")

                      #########################################
######################## (2) READING LIST OF SURVIVED ATOMS  ###############################################################
                      #########################################

''' 
   Reading the "listatomsfile" for the knowledge of the number and which atoms survive. It is organized in two columns: 
 
     -1st-  at_number of the atom that survives 
     -2nd-  type (AT or CG) of the atom that survives.   

    This part of program returns an error and exit if the following condition are not fullfilled: 

     o  The column file must not contain empty rows 
     o  The list of survived atoms is an integer, therefore we want to avoid strings inside the file 
     o  The indexes of the survived atoms has to be between 1 and n_atoms. The code can not recognize indexes less than 1 
        and higher than the total number of atoms: actually the at_number goes between 1 and n_atoms 
     o  A specific index of survived atoms can not be repetead. Each element must be exclusive. 
'''

#2.1- Checking if listatomsfile is present 
if not os.path.isfile(listatomsfile):
    print("Error while opening the file. '{}' does not exist\n".format(listatomsfile))
    quit()


#2.2- Opening and reading the file, and printing errors 
f_list_survived  = open(listatomsfile,"r")     
check_empty_file(listatomsfile)

count_line = 1

for line in f_list_survived: 
    check_empty_rows_CANVAS(line, listatomsfile)
    line = line.split() 

    if(len(line) != 2):             #N_lines = 2 
        print("The input file '{}' cannot be read correctly.".format(listatomsfile))
        print("\nError. Each line must contain 2 columns: the {}th row has {} columns".format(count_line, len(line)))
        print("\nBe sure that '{}' is organized in two different columns\n".format(listatomsfile))
        print(" o 1st column       INT                List of atomistic numbers of each survived atom.")
        print(" o 2nd column       STR                Type of each survived atom. The strings accepted are 'at' and 'cg'")
        print("\nLook below for further information\n")
        print_help_main_CANVAS()
        quit()

    else: 
        if not line[0].isdigit():  #1st column = INT
            print("The input file '{}' cannot be read correctly.".format(listatomsfile))
            print("\nError. The 1st column must be an integer number. '{}' at {}th row is not". format(line[0], count_line))
            print("\nBe sure that '{}' is organized in two different columns\n".format(listatomsfile))
            print(" o 1st column       INT                List of atomistic numbers of each survived atom.")
            print(" o 2nd column       STR                Type of each survived atom. The strings accepted are 'at' and 'cg'")
            print("\nLook below for further information\n")
            print_help_main_CANVAS()
            quit()

        if(line[1]!="CG" and line[1]!="cg" and line[1]!="Cg" and line[1]!="AT" and line[1]!="at" and line[1]!="At"):  #2nd column = 'AT' or 'CG' string 
            print("The input file '{}' cannot be read correctly.".format(listatomsfile))
            print("\nError. The 2nd column must contain the string AT (in case of atomstic) or CG (in case of Coarse Grained). '{}' at {}th row is not"\
                  .format(line[1], count_line))
            print("\nBe sure that '{}' is organized in two different columns\n".format(listatomsfile))
            print(" o 1st column          INT                List of atomistic numbers of each survived atom.")
            print(" o 2nd column          STR                Type of each survived atom. The strings accepted are 'at' and 'cg'")
            print("\nLook below for further information\n")
            print_help_main_CANVAS()
            quit()

    count_line = count_line + 1 


f_list_survived.seek(0)                        # read again file if everything went fine. 

dict_survived = {int(line.split()[0]) : line.split()[1] for line in f_list_survived} 

for k,v in dict_survived.items():
    if(k<0 or k>n_atoms):
        print("Error. The list of survived atoms '{}' must have index between 1 and {}".format(listatomsfile, n_atoms))
        print("\nLook below for further information\n")
        print_help_main_CANVAS()
        quit()

f_list_survived.seek(0)                        # read again file if everything went fine.
n_lines_f = sum(1 for line in f_list_survived) # compute n lines of listatomsfile

if(n_lines_f != len(dict_survived)):
    print("Error. The list of survived atoms '{}' contains, at least, two identical indexes at 1st column".format(listatomsfile))   # Avoid repeated elements
    print("\nLook below for further information\n")
    print_help_main_CANVAS()  
    quit()


#2.3 Number of survived atoms  
n_survived  = len(dict_survived)     

print("\nlist of survived atoms correctly read... 4% completed\n")

                      ###############################
######################## (3) READING TOPOLOGY FILE  ##############################################################
                      ###############################

"""
   Reading the topology file using 'read_gromacs.readtop' module. Topology file is divided in seven sub-parts: 

   1- Reading the atoms with their own properties (specifically the charge, the mass, and the at_type of each element). 
   2- Reading the atoms involved in bonds (couple of atoms) in the fullyAT system. 
   3- Reading the atoms involved in angles (triplets of atoms) in the fullyAT system. 
   4- Reading the atoms involved in dihedrals (quartet of atoms) in the fullyAT system. 
   5- Reading the atoms involved in pairs (couple of atoms) in the fullyAT system. 
   6- Reading the path where the forcefield is present. 
   7- Reading the forcefield used for simulating systems. 
"""

#3.1- Checking if topfile is present  
if not os.path.isfile(topfile):
    print("Error while opening the file. '{}' does not exist\n".format(topfile))
    quit()

#3.2- Opening and reading topology file
ft = open(topfile, "r")
check_empty_file(topfile)

toop = readtop(ft)

atoms_topol                          = toop[0]
original_bonds_protein_fully_at      = toop[1]
original_angles_protein_fully_at     = toop[2]
original_dihedrals_protein_fully_at  = toop[3]
original_pairs_protein_fully_at      = toop[4]
path_ff                              = toop[5]
forcefield                           = toop[6]
original_cmaps_protein_fully_at      = toop[7]  # In case of amber-forcefield, this list is empty


Flag_cmaps = True
if not original_cmaps_protein_fully_at: # If CMAP list is empty:
   Flag_cmaps = False

print("\nFully-atomistic topology file topol.top correctly read... 5% completed\n")

                      ##################################
######################## (4) READING FORCEFIELD PATH  ##############################################################
                      ##################################

"""
   Checking if the forcefield path found after reading the topology file is correct. Indeed, the 'path_ff' found before 
   is correct only if the forcefield directory is inside the main gromacs folder. If not, the forcefield path shoulb be 
   inside our current simulating directory.  
"""

path_ff_new    = check_forcefield_dir(forcefield, path_ff) 

                      ##################################
######################## (5) READING FFNONBONDED FILE ##############################################################
                      ##################################
"""
   Reading ffnonbondefile through 'read_gromacs.read_ffnonbonded'. It is not an external file, but it can be read inside the forcefield path 
   in GROMACS. If not, the program asks the user to insert the correct path. 'ffnonbondedfile' is divided in four subparts: 

   1- Reading a dictionary in which, for each atom type, is associated the value of sigma (dict_sigma) 
   2- Reading a dictionary in which, for each atom type, is associated the value of epsilon (dict_epsilon) 
   3- Reading a dictionary in which, for each atom type, is associated the value of mass (dict_mass)
   4- Reading a dictionary in which, for each atom type, is associated the value of atom number (dict_atnum) 
"""

#5.1- Finding the correct path for GROMACS ffnonbonded.itp file

ffnonbondedfile = check_ffnonbonded(path_ff_new) 

#5.2- Opening and reading ffnonbonded.itp file 
f_nb=open(ffnonbondedfile, "r")

ffnonbon = read_ffnonbonded(f_nb)

dict_sigma   = ffnonbon[0]
dict_epsilon = ffnonbon[1]
dict_mass    = ffnonbon[2]
dict_atnum   = ffnonbon[3]
pairtypes    = ffnonbon[4]   # In case of amber-forcefield, this list is empty. 

Flag_pairtypes = True 

if not pairtypes: # If pairtypes list is empty:
   Flag_pairtypes = False  

print("\nffnonbonded.itp correctly read... 8% completed\n")

                      ###############################
######################## (6) READING FFBONDED FILE ##############################################################
                      ###############################
"""
   Reading ffbondefile through 'read_gromacs.read_ffbonded'. It is not an external file, but it can be read inside the forcefield path 
   in GROMACS. The ffbondedfile is divided in five subparts: 

   1- Reading a list containing the bonded parameters for each kind of bondtype: [a1, a2, bon_lenght, const_lenght] 
   2- Reading a list containing the bonded parameters for each kind of angletype: [a1, a2, a3, theta, theta_ener] 
   3- Reading a list containing the bonded parameters for each kind of dihedraltype with funct=4: [a1, a2, a3, a4, 4, phi_0, k_phi, mult] 
   4- Reading a list containing the bonded parameters for each kind of dihedraltype with funct=9: [a1, a2, a3, a4, 9, phi_0, k_phi, mult]
   5- Reading a list containing the bonded parameters for each kind of dihedraltype called torsions: [string_tors, phi_0, k_phi, mult] 
"""

#6.1- Finding correct path of ffbonded.itp file 

ffbondedfile = check_ffbonded(path_ff_new) 

#6.2- Opening and reading ffbonded.itp file
f_b=open(ffbondedfile, "r")

ffbon = read_ffbonded(f_b)

bondtypes          = ffbon[0]
angletypes         = ffbon[1]
dihedraltypes_dih4 = ffbon[2]     # Found in Amber ff.
dihedraltypes_dih9 = ffbon[3]     # Found in Amber and Charmm ff.
dihedraltypes_tors = ffbon[4]     # Found in Amber99sb-ildn ff. 
dihedraltypes_dih2 = ffbon[5]     # Found in Charmm ff. 

print("\nffbonded.itp correctly read... 9% completed\n")


                     ###############################
######################## (7) READING CMAP FILE     ##############################################################
                      ###############################
"""
   Reading cmap file through 'check_cmap'. It is not an external file, but it can be read inside the forcefield path 
   in GROMACS.
"""

if Flag_cmaps == True: 
   cmapfile = check_cmap(path_ff_new)

   f_cmap = open(cmapfile, "r")
   cmaptypes = readcmap(f_cmap)

   print("\ncmap.itp correctly read... 10% completed\n")

                    ##############################################################
####################### (8) SPLIT ATOM LIST in SURVIVED and NOT-SURVIVED ATOMS  #############################################
                    ##############################################################
"""
   Updating "atom" with type AT-CG-NS (where NS stands for not-survived) 
"""

list_survived_atoms = [ i for k,v in dict_survived.items() for i in atom if k==i[3]]
list_other_atoms    = [ i for i in atom if i not in list_survived_atoms]  


list_survived_atoms = [ i + [v] for k,v in dict_survived.items() for i in list_survived_atoms if k==i[3]]
list_other_atoms    = [ i + ["ns"] for i in list_other_atoms]


                    #################################################################
####################### (9) UPDATE ATOM LIST with AT_TYPE of each survived atom    #############################################
                    #################################################################

"""
   Updating "atom" with the at_type of each survived atom
"""

# NOTE: j[0] of "atoms_topol" corresponds to i[3]-1 of "list_survived_atoms"
# Example: the element of index i[3] = 2000 can be found in atoms_topol[1999] 
#          Since we would like to append the at_type of each survived atom (column j[1] of atoms_topol), if i[3] = 2000, we append atoms_topol[1999][1]

# Generalizing, we append atoms_topol[i[3]-1][1]

list_survived_atoms = [ i + [atoms_topol[i[3]-1][1]] for i in list_survived_atoms]
list_other_atoms    = [ i + [atoms_topol[i[3]-1][1]] for i in list_other_atoms]


                    ######################################################################################
####################### (10)  UPDATE LIST of SURVIVIVED ATOMS with its own REF_BLOCK                    #############################################
                    ######################################################################################

"""
   Updating "list_survived_atoms" with the ref_block of each survived atoms from 1 to n_survived 
"""
list_survived_atoms = [list_survived_atoms[count-1] + [count]  for count in range(1, n_survived+1)]

                    ######################################################################################
####################### (11) UPDATE LIST ATOMS with the REF_BLOCK of each survived atom                  #############################################
                    ######################################################################################

"""
   Updating "atom" with ref_block of survived atoms based on the minimum distance. For each NS atom, 
   the ref_block is the same of its closest CG survived atom 
"""

#11.1- Compute "min_distances" list of lists that contains

index=0
min_distances=[]
for i in list_other_atoms:
    dist=[]
    min_dist=1000000
    for j in list_survived_atoms:

        if(j[7]=="cg" or j[7]=="CG" or j[7]=="Cg"):  # We require that the survived atom of reference has to be CG.
            x2=(j[4]-i[4])**2
            y2=(j[5]-i[5])**2
            z2=(j[6]-i[6])**2
            distance=sqrt(x2+y2+z2)
            dist.append(distance)
            if(distance<min_dist):
                min_dist=distance
                index=j[-1]


    list_dist=[min_dist,index]
    min_distances.append(list_dist)



""" MUCH FASTER, BUT IN CASE OF BIG SYSTEMS, A MEMORY ERROR OCCURS; THEREFORE THE PREVIOUS METHOD IS PREFERABLE, DESPITE IT IS SLOWER
#11.1- Compute "min_distances" list of lists that contains

# NOTE: "Else 10000" is due to the fact that if the survived atom is AT, the distance between the other_atoms and the AT-survived-atom
# cannot never be the minimum one. In this way the AT-survived-atoms cannot include other atoms by definition. The other atoms can only 
# belong to CG-survived-atoms, and never to the AT ones. Indeed, some decimated atoms, could have as a minimum distance an atomistic 
# survived atom. But we cannot allow it (through ELSE 10000)  

x_surv  = np.array([i[4] if (i[7]=="cg" or i[7]=="CG" or i[7]=="Cg") else 10000 for i in list_survived_atoms])
y_surv  = np.array([i[5] for i in list_survived_atoms])
z_surv  = np.array([i[6] for i in list_survived_atoms])

x_other = np.array([i[4] for i in list_other_atoms])
y_other = np.array([i[5] for i in list_other_atoms])
z_other = np.array([i[6] for i in list_other_atoms])

dx = x_surv - x_other[:, np.newaxis]
dy = y_surv - y_other[:, np.newaxis]
dz = z_surv - z_other[:, np.newaxis]

x2 = np.square(dx)
y2 = np.square(dy)
z2 = np.square(dz)

distance = np.sqrt(x2+y2+z2) 

argmin  = np.argmin(distance, axis = 1)  # position of each minimum of the subarray of each distance with respect the survived atom of reference  
minimum = np.amin(distance, axis = 1)    # minimum for the 2D-array for each row

ref_block = []
for i in argmin: 
    ref_block.append(list_survived_atoms[i][-1])

refblock = np.array(ref_block)

combined = np.transpose((minimum, refblock))

min_distances = combined.tolist()

for i in min_distances:
    i[1] = int(i[1])
"""

#11.2- Append the ref_block of each NS atom according with its closest CG survived atom 

for i in range(1, len(list_other_atoms)+1):        # from 1 to the number of not-survived atoms 
    ref_block = min_distances[i-1][1]
    list_other_atoms[i-1].append(ref_block)

                    ###############################################################
####################### (12) UPDATE ATOM LIST with the MASS of each atom         #############################################
                    ###############################################################
"""
   Updating "atom" with the mass of each atom. The latter can be found in the column[7] of the list "atoms_topol"
"""

list_survived_atoms = [ i + [atoms_topol[i[3]-1][7]] for i in list_survived_atoms]
list_other_atoms    = [ i + [atoms_topol[i[3]-1][7]] for i in list_other_atoms]


                    ###############################################################
####################### (13) UPDATE ATOM LIST with the CHARGE of each atom       #############################################
                    ###############################################################

"""
   Updating "atom" with the charge of each atom:
      o The charge of each NS survived atom (list_other_atoms) remains the same 
      o The charge of each AT survived atom remains the same because the atom is treated atomistically 
      o The charge of each CG survived atom is the algebraic sum of the atoms charge that belong at the same ref_block
"""

#13.1- Update "list_survived_atoms" and "list_other_atoms" with the original value of charge (column[6] of "atoms_topol) 

list_survived_atoms =  [i + [atoms_topol[i[3]-1][6]] for i in list_survived_atoms]
list_other_atoms    =  [i + [atoms_topol[i[3]-1][6]] for i in list_other_atoms]


#13.2- We substitute the value of charge ONLY for the CG survived atoms, as the new value is the sum of the atoms charge 
#      that belong at the same ref_block. 


for block in range(1,n_survived+1):

    q = 0 

    for i in list_other_atoms: 
        if(block==i[9]):                				# i[9] is the ref_block of list_other_atoms 
            q = q + i[-1]

    for j in list_survived_atoms: 	                                # j[9] is the ref_block of list_survuved_atoms 
        if(block==j[9]):
            q = q + j[-1]						# Q is now the sum of charges of the decimated atoms 
 									# and the survived atom that belong to a particular block 

    for k in list_survived_atoms:                                       # Change the value of charge k[-1] for each survived atom 
        if(block==k[9] and (k[7]=="cg" or k[7]=="CG" or k[7]=="Cg")):   # with the value of the charge total sum of each block   
            k[-1] = float("{:3.6f}".format(q))    


                    ##############################################################
####################### (14) UPDATE ATOM LIST with the SIGMA of each atom       #############################################
                    ##############################################################

"""
    Updating "atom" with the sigma of each atom.
    **Sigma** is a measure of how much close two atoms are, therefore it corresponds at the Van der Waals radius of a given atom 

    Two cases are possible:
      o In case the block consists of only one atom (as in case of AT-survived atoms), sigma conserves its original value 
        (look into ffnonbonded.itp file) 
      o In case the block consists of two or more atoms, sigma is equals to the double of Radius of Gyration (Rg) as sigma is 
        a diameter. 

    In the latter case, it turns out that the Radius of Gyration squared is given by: 
                         --------------------------------
                        |  Rg^2 = (1/N)*sum(r_i-r_cm)^2  |
                         --------------------------------     
    where:
      o  N    = total atoms of considered block 
      o  r_i  = position of each atom of the block (x,y,z)
      o  r_cm = center of mass of the spefic block: r_cm = (1/N)*sum(r_i)
"""

#14.1- Calculating, the center of mass of each block (r_cm) 
for block in range(1,n_survived+1):

    x_cm = 0
    y_cm = 0
    z_cm = 0 

    n_atoms_block = 0 						# number atoms in each specific block 

    for i in list_other_atoms:
        if(block==i[9]):                                        # i[9] is the ref_block of list_other_atoms 
           
            x_cm = x_cm + i[4]					# i[4] is the value of x of each atom 
            y_cm = y_cm + i[5]					# i[5] is the value of y of each atom 
            z_cm = z_cm + i[6]  				# i[6] is the value of z of each atom 

            n_atoms_block = n_atoms_block + 1 

    for j in list_survived_atoms:                               # j[9] is the ref_block of list_survived_atoms 
        if(block==j[9]):

            x_cm = x_cm + j[4]                                  # j[4] is the value of x of each atom 
            y_cm = y_cm + j[5]                                  # j[5] is the value of y of each atom 
            z_cm = z_cm + j[6]                                  # j[6] is the value of z of each atom 

            n_atoms_block = n_atoms_block + 1

    x_cm = x_cm / n_atoms_block
    y_cm = y_cm / n_atoms_block
    z_cm = z_cm / n_atoms_block

    
#14.2- ...and computing Rg^2 for each block
    Rg_squared = 0 
  
    for i in list_other_atoms: 
        if(block==i[9]):

            x_squared  = (i[4] - x_cm)**2
            y_squared  = (i[5] - y_cm)**2
            z_squared  = (i[6] - z_cm)**2

            Rg_squared = Rg_squared + (x_squared + y_squared + z_squared)


    for i in list_survived_atoms:
        if(block==i[9]):

            x_squared  = (i[4] - x_cm)**2
            y_squared  = (i[5] - y_cm)**2
            z_squared  = (i[6] - z_cm)**2
   
            Rg_squared = Rg_squared + (x_squared + y_squared + z_squared)


    Rg_squared = Rg_squared / n_atoms_block

    Rg = sqrt(Rg_squared)

    Rg = Rg * 2     							     # because sigma is a diameter and not a radius (nm)  

    Rg = float("{:12.9f}".format(Rg))


#14.3- Appending the value of sigma in list_survived_atoms accordings with the rules of one or more atoms per block  
    for i in list_survived_atoms: 
        
        if(n_atoms_block!=1 and block==i[9]):
            i.append(Rg)

        if(n_atoms_block==1 and block==i[9]):
            for k,v in dict_sigma.items(): 
                if(k==i[8]):                               # if the at_type is the same in dict_sigma and list_survived_atoms
                    i.append(v)  			   # append the original value of sigma of the specific at_type   


#14.4- Calculating the maximum value of sigma, since when creating .mdp files, the coulomb and VdW cutoff should be equals to 2.5 * max_sigma  
max_sigma = 0 
for i in list_survived_atoms:
    if(i[12] > max_sigma):
        max_sigma = i[12]


                    ################################################################
####################### (15) UPDATE ATOM LIST with the EPSILON of each atom       #############################################
                    ################################################################

"""
   Updating "atom" with the epsilon of each atom.

   Two cases are possible: 

      o The block consists of one or more atoms whose epsilons are all equals to zero. In this case: eps_block = 0 

      o The block consists of one or more atoms, in which at least, the epsilon of one of them, is not equals to zero.
        In this case, the epsilon of the survived atom, representing the specific block, is given by: 

             ---------------------------------                 --------------------------------- 
        	Eps_block = [Prod(Eps_i)]^(1/N) or equivalently  Eps_block = e^[(1/N)Sum(Eps_i)]
             ---------------------------------                 ---------------------------------
    where: 

      o Eps_i is the epsilon !=0 of an atom (survived or not) that belongs to a specific block.
      o N     is the number if atoms, whose epsilon is != 0, of a specific block  

    The condition that epsilon has to be not-equal to zero, stems from the fact that, if Eps_i = 0 => Eps_block = 0,
    even though the block consists of more atoms with Eps_i not-equal to 0    

    We can also notice that, at difference with 'Sigma', if the block consists of only one atom,
    Eps_block = Eps_atom, therefore the formula used works good in both explained cases. 
"""

#15.1- Calculating the original value of epsilon for each atom comparing with dict_epsilon, and adding it in "list_survived_atoms"  
#      and "list_other_atoms". 

list_other_atoms    = [i + [v] for i in list_other_atoms for k,v in dict_epsilon.items() if(k==i[8])]
list_survived_atoms = [i + [v] for i in list_survived_atoms for k,v in dict_epsilon.items() if(k==i[8])]

#15.2- Calculating la value of epsilon for each block. For simplicity it is better to compute the logarithm of epsilon.
#      The reason is because the N-times product of very small very quantities (like epsilon) is equal to a very small 
#      result. A solution is the usage of "mpf" module of library mppmath. The best solution is the log calculation without 
#      needing to add any library or working with too small values. 

for block in range(1,n_survived+1):

    ln_eps    = 0 
    count_eps = 0

    for i in list_other_atoms:
        if(block==i[9] and i[-1]!=0):                         # i[9] is the ref_block of list_other_atoms 

            ln_eps = ln_eps + log(i[-1])  
            count_eps = count_eps + 1                         # count_eps increases ONLY if the epsilon value is not-equal to Zero.

  
    for i in list_survived_atoms: 
        if(block==i[9] and i[-1]!=0):

            ln_eps = ln_eps + log(i[-1]) 
            count_eps = count_eps + 1

    if(count_eps !=0):                                           
        eps=float(exp((1/count_eps)*ln_eps))                  # if count_eps != 0 means that it exists at least one atom, whose eps is not-equal to zero. 
							      # Therefore, Eps_block = e^[(1/N)Sum(Eps_i)]    (2nd case)
    else:
        eps=float(0)                                          # if count_eps = 0 means that all the atoms included in the block have eps = 0 
						              # Therefore, Eps_block = 0   (1st case) 
    for i in list_survived_atoms: 
        if(block==i[9]): 

            i[-1] = eps                                       # Substitute the orignal value of epsilon for each block.


                    #####################################################################################################
####################### (16) UPDATE ATOM LIST with the number of survived atoms for each residue (four, one, or atom)  #########################
                    #####################################################################################################
"""
   Updating "atom" with the number of survived atoms for each atomistic residue (four, one, or atom). Three cases are possible: 

     o We retain the backbone atoms of the specific residue, that are FOUR atoms (N, C, CA, and O). In this case we append the string "four"
     o We retain only the C_alpha atom of the specific residue, that is ONE atom (CA). In this case we append the string "one" 
     o The residue is completely atomistic. In this case we append the string "atom". 
"""

list_coun = [i[0] for i in list_survived_atoms]              # Creation of auxiliary list containing only the res_number of each survived atom 

dict_coun = {i:list_coun.count(i) for i in list_coun}        # Creation of a dictionary containing  the number of survived atoms for each residue. 
							     # key = res_number, value = number of surv atoms for each res_number 
							     # In case of at residues : value > 4 
						             # In case of CA cg       : value = 1 
							     # In case of backbone cg : value = 4 

list_survived_atoms = [i + [v] for i in list_survived_atoms for  k,v in dict_coun.items() if (k==i[0])]

for i in list_survived_atoms: 
    if(i[-1]==1): 
        i[-1] = "one"

    elif(i[-1]==4): 
        i[-1] = "four"

    else: 
        i[-1] = "atom" 

                    ####################################################################################################
####################### (17) UPDATE ATOM LIST with the atomic number of each survived atom                            #########################
                    ####################################################################################################

"""
   Updating "atom" with the atomic number (H=1, C=6, and so on...) for each survived atom. It can be found it in dict_atnum. 
"""

list_survived_atoms = [ i + [v] for i in list_survived_atoms for k,v in dict_atnum.items() if (i[8]==k)]

                   #####################################################################################################
####################### (18) UPDATE ATOM LIST with the new atom type of each survived atom                            #########################
                    ####################################################################################################

"""
   Updating "atom" with the new atom_type for each survived atom. Two cases are possible: 

      o AT-survived atom: The original atom type does not change, since the atom is completely atomistic. 

      o CG-survived atom: the original atom type is replaced by the new one. Specifically, the atom type of 
                          this bead given by a letter (starting from A) followed by a number (from 101 to 9999):
 			   i.e. A101, A2,...,A9999,  B101,...,B9999,  ,....,  Z101,...,Z9999  

      Please Note that the enumeration starts from 101 and not 1, because could create some artifacts. 
      For instance, C3 is a Carbon already defined in the atomistic force. 
      Starting from 101, avoid these kind of issues.  
"""

ch = "A"   # initilialized character 

j = 101

for i in list_survived_atoms:
 
    if(i[7]=="cg" or i[7]=="CG" or i[7]=="Cg"): 

        new_attype = ch + str(j)

        if(len(new_attype)<6):                      # i.e the value of J overcomes 9999 (letter + number lenght higher than 5)
            i.append(new_attype)
            j = j + 1
 
        else: 
            j = 101
            ch = chr(ord(ch) + 1)
            new_attype = ch + str(j)      
            i.append(new_attype)
            j = j + 1  

 
    if(i[7]=="at" or i[7]=="AT" or i[7]=="At"):
       i.append(i[8])                                # i.e. in this case the atom type does not change, we append the value of original at_type  
						     # already written in the column[8] 

                   #####################################################################################################
####################### (19) UPDATE ATOM LIST with the new residue number of each survived atom                       #########################
                    ####################################################################################################

"""
   Updating "atom" with the new residue number (i.e. the 1st column for the .gro coordinate file). Specifically:

      o For each CG survived bead, the counter "i" increases.
      o For each AT survived bead, the counter "i" changes, increases by 1, only when the residue changes. 
"""

count = 0 

for i in list_survived_atoms: 

    if(i[7]=="cg"): 
        count = count + 1 
        i.append(count)     

    if(i[7]=="at"):
        if(i[2]=="N"):         # if at_name is "N" it means that a new amino acid is started. Evert residue begins with "N"
            count = count + 1
            i.append(count) 

        else:                  # Otherwise, the atomistic amino acid has not changed, therefore count does not change. 
            i.append(count)   

print("\nlist of survived atoms updated with all properties... 27% completed\n")

                     #######################################################
######################## (20) CHECKING WHERE GROMACS IS INSTALLED         ###############################################################
                     ####################################################### 

"""
   So far, we have already found where the forcefield is present, therefore it means that Gromacs is correctly installed in 
   our laptop/cluster. Now, we just need to find the module GMXRC. 
   Usually, it is present is the folder /PREFIX_FOLDER/bin/

   We already know: path_ff    = PREFIX_FOLDER/share/gromacs/top/FF_FOLDER.ff 
                  : forcefield = FF_FOLDER.ff

   Therefore PREFIX_FOLDER is equal to the difference between path_ff, forcefield, and share/gromacs/top  
   Finally, path_GMXRC = PREFIX_FOLDER + /bin

   Check also for gmx command that could be 'gmx' or 'gmx_mpi'. We need it for solvating system and adding ions.   
"""

path_GMXRC = check_for_GRXRC(path_ff, forcefield)
        
os.system("source {}".format(path_GMXRC))

gmx_command = check_for_gmx_command(path_ff, forcefield) 


          ###############################################################
############# (21) WRITING NEW .GRO FILE ACCORDING WITH MODIFICATIONS #############################################################################
          ###############################################################

'''
   Writing two files: 

   o **outsurv.gro** that contains ONLY the coordinates of survived atoms ready to be solvated.

   o **initial-frame.gro**, an initial coordinate file, that contains ONLY the coordinates of survived atoms, but with the new 
       at_types and res_name for CG-survived atoms. The remeinder does not change. This file is useful for a comparison between fully-AT 
       and CANVAS simulation because it contains ONLY the survived atoms of protein, with the original value of coordinates. 
       The point is that, after solvating the system, the protein coordinates are shifted, not allowing a perfect comparison, 
       for instance in terms of RMSD or RMSF between fullyAT and CANVAS model.  
'''

#21.1- Writing outsurv.gro 

fw2=open("outsurv.gro", "w+")          # w+ means that the file can be written and read, as well. If the file exists, will be removed and created again. 
                                       # https://www.geeksforgeeks.org/reading-writing-text-files-python/

write_outsurv_file(fw2, n_survived, list_survived_atoms, Lx, Ly, Lz)
fw2.close()

#21.2- Writing initial-frame.gro

f_init = open("initial-frame.gro", "w") 
write_initial_frame_file(f_init, n_survived, list_survived_atoms, Lx, Ly, Lz)
f_init.close()


          ###############################################################
############# (22) SOLVATING PROTEIN WITH GROMACS                     ############################################################
          ###############################################################

#22.1- Creation of box.gro and solvated_survived_1st.gro. Moreover:  > dev/null 2>&1 hides the gromacs output on screen. 
os.system("{} editconf -noversion -quiet -f outsurv.gro -o box.gro -c -d 1.4 -bt cubic > /dev/null 2>&1".format(gmx_command))

if(solvate_Flag == True):
    os.system("{} solvate -noversion -quiet -cp box.gro -cs spc216.gro -o solvated_survived.gro > /dev/null 2>&1".format(gmx_command))

    #22.2- Opening and reading the solvated file just created
    f_s = open("solvated_survived.gro", "r")

    gro_solv = readgro_all(f_s)

    atoms            = gro_solv[0]    #atoms = [res_number, res_name, at_name, at_number, pos_x, pos_y, pos_z]
    n_solvated_atoms = gro_solv[1]
    Lx               = gro_solv[2]
    Ly               = gro_solv[3]
    Lz               = gro_solv[4]

    CA_list = [x[3] for x in atoms if x[2]=="CA"] 

    f_solv_new = open("newsolvated.gro", "w+")
    write_newsolvated_file(f_solv_new, atoms, list_survived_atoms, n_solvated_atoms, Lx, Ly, Lz)
    f_solv_new.close() 


if(solvate_Flag == False):
    f_box = open("box.gro", "r")
 
    gro_box = readgro_all(f_box) 
    
    atoms_box   = gro_box[0]    #atoms = [res_number, res_name, at_name, at_number, pos_x, pos_y, pos_z] 
    n_box_atoms = gro_box[1]
    Lx          = gro_box[2]
    Ly          = gro_box[3]
    Lz          = gro_box[4]

    CA_list = [x[3] for x in atoms_box if x[2]=="CA"]

    f_box_new = open("newbox.gro", "w+")
    write_newsolvated_file(f_box_new, atoms_box, list_survived_atoms, n_box_atoms, Lx, Ly, Lz)
    f_box_new.close() 


#22.3- Counting the number of SOL atoms 
count_SOL = 0

if(solvate_Flag == True): 
    for i in atoms:
        if(i[2]=="OW"):
            count_SOL = count_SOL + 1
    print("\nSystem solvated... 33% completed\n")

if(solvate_Flag == False): 
    count_SOL = 0 
    print("\nSystem not solvated as required... 33% completed\n")


          ########################
############  (A) BONDS         ###################################################
          ########################

'''
   Knowledge of the number of bonds and which bonds are present in the CANVAS model. 
   Bonds between two survived atoms could be atomistic or not. Two cases are possible: 

   1. Two survived atoms, in the fully atomistic representation are covalently bonded, therefore the same atomistic bond 
      must be reproduced in the Mul-Res representation. 
 
   2. Two survived atoms, are not bonded in the fully-at case, but their distance il less than  the cutoff equals to 1.4 nm. 
      In such case, two atoms will be connected by a spring constant if they belong at two consecutive residues; otherwise they will be 
      linked by a weak variable elastic constant. An exhaustive discussion of the interactions in the CANVAS model (in particular in case the model 
      presents also the backbone beads, and not only C-alphas) is presented in the reference paper. 
'''

#A.1- After reading the topology file, and which couples of atoms have covalent bond in the fullyAT representation
#     (original_bonds_protein_fully_at = toop[1]), we check if two survived atoms keep this property. 

f_list_survived.seek(0) 
list_survived_only_num = [int(line.split()[0]) for line in f_list_survived]     # list 1st column survived atoms 

bond_prot_fully_at = []


for i in original_bonds_protein_fully_at:
    result=all(x in list_survived_only_num for x in i)
    if(result==True):
        bond_prot_fully_at.append(i)         # list survived atoms that have a covalent bond (original index of fullyAT representation) 



#A.2- Appending, for each couple of survived atoms the "at_type" and "at/cg type" (i[8], and i[7] of list_survived_atoms, respectively)   
#     [a1, a2, at_type-1, at_type-2, at/cg-1, at/cg-2]

bond_prot_fully_at = [i + [k[8]] for i in bond_prot_fully_at for k in list_survived_atoms if (k[3] == i[0])]
bond_prot_fully_at = [i + [k[8]] for i in bond_prot_fully_at for k in list_survived_atoms if (k[3] == i[1])]
bond_prot_fully_at = [i + [k[7]] for i in bond_prot_fully_at for k in list_survived_atoms if (k[3] == i[0])]
bond_prot_fully_at = [i + [k[7]] for i in bond_prot_fully_at for k in list_survived_atoms if (k[3] == i[1])]


#A.3- After reading the bondtype in ffbonded file (ffbon[0]), we append the value if b0 and kb for each couple of 
#     survived atoms bonded covalently. 'bondtypes' and 'bond_prot_fully_at' look like as follow: 
#
#               bond_prot_fully_at : [[1,2,X,Y,at,cg],[2,3,A,B,cg,cg],[...]]
#               bondtypes          : [[X,Y,b0,kb],[...]] 
 
for i in bond_prot_fully_at:
    for j in bondtypes:
        if((j[0]==i[2] and j[1]==i[3]) or (j[0]==i[3] and j[1]==i[2])):

            i.append(float(j[3]*836.8))     # Append kb value (in Groamacs units)
            i.append(float(j[2]/10))        # Append b0 value (in Gromacs units)


#A.4- Reconverting the original index (of fullyAT representation) in the new index after decimated the most of atoms in CANVAS model
#     and substitute the original value (1st and 2nd column) 

for i in list_survived_atoms: 
    for j in bond_prot_fully_at:
        if(i[3] == j[0]): 
            j[0] = i[9]

for i in list_survived_atoms:
    for j in bond_prot_fully_at:
        if(i[3] == j[1]):
            j[1] = i[9]

#A.5- After checking for atoms covalently bonded we write all the other couple  of atoms (or beads) that will be connected 
#     by stiff or weak (variable) spring. This result will be stored in BOND_SURV_PROTEIN = []

	#A.5.1- Defining a function depending on the value of elastic constant (k_el) and a second function where the value 
        #       of the elastic constant (k_el) is dependent on the distance between the couple of atoms involved in the bond. 

def add_variable_bonds_ENM():

    A        =  81.2273
    B        =  1322.06
    C        = -1.25948
    boltz    = 0.008314
    temperat = 300

    x = j[4] - i[4]
    y = j[5] - i[5]
    z = j[6] - i[6]

    dist = sqrt(x**2 + y**2 + z**2)

    if(dist <= 1.4):    # We insert a spring only IF the distance between the couple of atoms is less than 1.4 nm. 

        orig = []

        if(i[9]<j[9]):                  # because we want the 1st element starting from left has the lowest index. 
            k_el = (A+B*(dist**C))/(boltz*temperat)

            orig.append(i[9])           # a1
            orig.append(j[9])		# a2 
            orig.append(i[8])           # type-1
            orig.append(j[8])           # type-2
            orig.append(i[7])           # AT/CG-1 
            orig.append(j[7])           # AT/CG-2 
            orig.append(k_el)           # k_el 
            orig.append(dist)           # b0

        else:
            k_el = (A+B*(dist**C))/(boltz*temperat)

            orig.append(j[9])
            orig.append(i[9])
            orig.append(j[8])
            orig.append(i[8])
            orig.append(j[7])
            orig.append(i[7])
            orig.append(k_el)
            orig.append(dist)

        bonds_surv_proteina.append(orig)


def add_bonds_ENM(k_el):

    x = j[4] - i[4]
    y = j[5] - i[5]
    z = j[6] - i[6]

    dist = sqrt(x**2 + y**2 + z**2)

    if(dist <= 1.4):   

        orig = []

        if(i[9]<j[9]):       
            orig.append(i[9])
            orig.append(j[9])
            orig.append(i[8])
            orig.append(j[8])
            orig.append(i[7])
            orig.append(j[7])
            orig.append(k_el)
            orig.append(dist)

        else:
            orig.append(j[9])
            orig.append(i[9])
            orig.append(j[8])
            orig.append(i[8])
            orig.append(j[7])
            orig.append(i[7])
            orig.append(k_el)
            orig.append(dist)

        bonds_surv_proteina.append(orig)


        #A.5.2- After defining the values of weak and stiff spring, we take care to put spring in which is involved at least of CG bead. 

#k_weak    = 500
k_strong  = 50000

k_interface = 50    # elastic constant on the interface AT-CG where the beads and atoms involved are not consecutive

bonds_surv_proteina = []

for i in list_survived_atoms: 
    for j in list_survived_atoms:
      
        if(j[9]>i[9]):                                                              # We avoid repeated elements (at_number of J > at_number of I)

            """ 1.1 Coarse Grained C-alpha and atomistic residue """

            if(i[7] =="cg" and i[14]=="one" and j[14]=="atom"):                     # i = cg-CA, j = at-residue 

                # If consecutive
                if(i[0] == j[0] - 1):                                               # i (cg-CA) precedes j (at-residue)
                
                    if(j[2]=="N"):                                                  # i = cg-CA, j = atomistic N 
                        add_bonds_ENM(k_strong)           
      

                # If NOT-consecutive --> ONLY interface (very weak) bond between cg-CA e at-CA 
                if(i[0] != j[0] - 1 and i[0] != j[0]):  

                    if(j[2] =="CA"):				                    # i = cg-CA, j = atomistic CA 
                        add_bonds_ENM(k_interface)
 
            """ 1.2 Atomistic residue and Coarse Grained C-alpha """

            if(i[14] == "atom" and j[7]=="cg" and j[14]=="one"):                    # i = at-residue; j = cg-CA  

                # If consecutive
                if(i[0] == j[0] - 1):                                               # i (at-residue) preceedes j (cg-CA)

                    if(i[2] == "C"):                                                # i = atomistic C, j = cg-CA 
                        add_bonds_ENM(k_strong)


                # If NOT-consecutive --> ONLY interface (very weak) bond between at-CA e cg-CA 
                if(i[0] != j[0] - 1 and i[0] != j[0]):

                    if(i[2] =="CA"):                                                # i = atomistic CA; j = cg-CA 
                        add_bonds_ENM(k_interface)


            """ 2. Coarse Grained C-alpha & Coarse-grained C-alpha """

            if(i[7] =="cg" and i[14]=="one" and j[7] == "cg" and j[14]=="one"):     # i = cg-CA ; j = cg-CA 

                # If consecutive
                if(i[0] == j[0] - 1):
                    add_bonds_ENM(k_strong)

                # If NOT-consecutive
                if(i[0] != j[0] - 1 and i[0] != j[0]):
                    add_variable_bonds_ENM()


            """ 3.1 Coarse-Grained C-alpha & Backbone-CG """ 

            if(i[7] =="cg" and i[14]=="one" and j[7]=="cg" and j[14] == "four"):    # i = cg-CA ; j = cg-backbone

                # If consecutive
                if(i[0] == j[0] - 1):   
               
                    if(j[2]=="N"):                                                  # i = cg-CA, j = backbone N 
                        add_bonds_ENM(k_strong)


                # If NOT-consecutive
                if(i[0] != j[0] - 1 and i[0] != j[0]):
       
                     if(j[2] == "CA"):
                        add_variable_bonds_ENM()

 
            """ 3.2 Backbone-CG & Coarse-Grained C-alpha """

            if(i[7] == "cg" and i[14]=="four" and j[7] =="cg" and j[14]=="one"):    # i = cg-backbone ; j = cg-CA

                # If consecutive
                if(i[0] == j[0] - 1):                                       

                    if(i[2] == "C"):                                                # i = backbone C, j = cg-CA 
                        add_bonds_ENM(k_strong)



                # If NOT-consecutive
                if(i[0] != j[0] - 1 and i[0] != j[0]):

                    if(i[2] == "CA"):
                        add_variable_bonds_ENM()


            """ 4. Backbone-CG & Backbone-CG """

            if(i[7] =="cg" and i[14]=="four" and j[7]=="cg" and j[14] == "four"):   # i = cg-backbone ; j = cg-backbone

            # If consecutive
                if(i[0] == j[0] - 1):  

                    if(i[2] == "C" and j[2] == "N"):                                # avoid stiff connection between the 1st N and the 2nd C, 
                        continue       	                                            # as they are covalently connected  


            # If NOT-consecutive 
                if(i[0] != j[0] - 1 and i[0] != j[0]): 

                    if(i[2] == "CA" and j[2] == "CA"):                              # ONLY weak bonds between the two CAs (guarantee less rigidity)
                        add_variable_bonds_ENM()


            """ 5.1 Backbone-CG & Atomistic residue """

            if(i[7] == "cg" and i[14] == "four" and j[14] == "atom"):               # i = backbone-cg ; j = at-residue 

                # If consecutive 
                if(i[0] == j[0] - 1):

                    if(i[2] == "C" and j[2] == "N"):              		    # avoid stiff connection between the 1st C and the 2nd N, 
                        continue                                 		    # as they are covalently connected 


                # If NOT-consecutive --> ONLY weak bond between backbone-CG e at-CA 
                if(i[0] != j[0] - 1 and i[0] != j[0]):

                    if(i[2] =="CA" and j[2] == "CA"):                               # i = backbone-CA-cg; j = atomistic CA
                           add_bonds_ENM(k_interface)  
       

            """ 5.2 Atomistic residue & Backbone-CG """

            if(i[14] == "atom" and j[7] == "cg" and j[14] == "four"):               # i = at-residue ; j = backbone-cg

            # If consecutive 
                if(i[0] == j[0] - 1):

                    if(i[2] == "C" and j[2] == "N"):                                # avoid stiff connection between the 1st C and the 2nd N, 
                        continue                                                    # as they are covalently connected 


            # If NOT-consecutive --> ONLY weak bond between atomistic CA e backbone-cg 
                if(i[0] != j[0] - 1 and i[0] != j[0]):

                    if(i[2] =="CA" and j[2] == "CA"):                                   # i = atomistic CA ; j = backbone-CA-cg
                        add_bonds_ENM(k_interface) 


bonds_surv_proteina = sorted(bonds_surv_proteina, key=itemgetter(0,1)) 						    # sorting 1st, and then 2nd column
bonds_surv_proteina =list(bonds_surv_proteina for bonds_surv_proteina,_ in itertools.groupby(bonds_surv_proteina))  # remove duplicates if present


        #A.5.3- Creating a new list called "bonds_mult_res_prot" which is the sum of two lists:
        #
        #       o bond_prot_fully_at : list of lists with all atomistic bonds (kb and b0 taken from ffbonded.itp file)
        #       o bonds_surv_proteina: list of lists with couple of atoms (and beads) not covalently bonded. At least one of two elements is a CG-bead.
    
bonds_mult_res_prot = bond_prot_fully_at + bonds_surv_proteina
bonds_mult_res_prot = sorted(bonds_mult_res_prot, key=itemgetter(6), reverse = True)
bonds_mult_res_prot = sorted(bonds_mult_res_prot, key=itemgetter(0,1))                    # sort 1st, 2nd column and k constant 
								                          # in ascending order (7th column).

        #A.5.4- Removing couple of atoms (and beads) if repetead.  

bonds_mult_res_prot = sorted(bonds_mult_res_prot, key = itemgetter(6), reverse=True)      # descending order sorting of elastic constant value
bonds_mult_res_prot = sorted(bonds_mult_res_prot, key = itemgetter(0,1))                  # "a1" and "a2" are ascending sorted 


        #A.5.5- Some couple of atoms in "bonds_mult_res_prot" could have a covalent bond, strong, and weak spring at the same time 
        #       therefore, we impose that covalent bond wins vs. strong bond, that wins vs. weak bond.
        #       We use the following code, defining that two or more internal list are duplicated if the first two element are equals.
        #       For example, [1,2,50] is equal to [1,2,70]. The code maintains the list founded as first. Indeed, the value of kb 
        #       has been sorted in descing order (x[6]), so that if two internal lists are "identical" (in the sense described before), 
        #       only the list with higher kb is kept.
        #       covalent bond > k_strong > k_weak (variable)

bonds_mult_res_prot=[list(my_iterator)[0] for g, my_iterator in itertools.groupby(bonds_mult_res_prot, lambda x: [x[0],x[1]])]   

print("\nlist of new bonds updated... 48% completed\n")

# A.6- Reading domfile: if Domain Division FILE is present, we exclude bonds between atoms of different domanain, where both CG bead is involved.
#      We cannot exclude CG-AT bonds, otherwise each domain is completely unconnected the the rest of system 

# A.6.1- Checking for empty file 
if domfile is not None:
    f_dom = open(domfile,"r") 
    check_empty_file(domfile)

    for line in f_dom:
        check_empty_rows_CANVAS(line, domfile)
    
    f_dom.seek(0)

    
    dom = [line.split() for line in f_dom]    # create list of lists from file 


    #A.6.2- Checking that only integer numbers are present 
    for i in dom: 
        for j in i:
            if(not j.isdigit()):
                print("Error. '{}' is not integer. Each domain must contain only the atomistic number of the atoms, not strings, not float.".format(j))
                print("\nLook below for further help\n")
                print_help_main_CANVAS()
                quit() 


    dom = [[int(j) for j in i] for i in dom]  # transform list of lists, from string in INT 


    #A.6.3- Creating dom$i_fullyAT and dom$i where i belongs to [1, n_domains] 
    for i, value in enumerate(dom):
        exec('dom%d_fullyAT = %s' % (i+1,value))

    n_domains = i + 1
    for i in range(1, n_domains+1):
        exec('dom%d = []' % i)

    
    for count in range(1, n_domains+1):                
        for i in list_survived_atoms:
            if eval('i[3] in dom%d_fullyAT and i[7] == "cg"' % count): 
                eval('dom%d.append(i[9])' % count)


    #A.6.4- Appending in 'bonds_mult_res_prot' the number of domain (1,2,...,N) in case of CG beads 
    for i in bonds_mult_res_prot:
        for count in range(1, n_domains+1):
            if eval('i[0] in dom%d' % count):
                i.append(count)

    for i in bonds_mult_res_prot:
        for count in range(1, n_domains+1):
            if eval('i[1] in dom%d' % count):
                i.append(count)

    #A.6.5- We keep: 
    #          o Bonds between CG beads having same domanin
    #          o Bonds where ONLY ONE CG bead is involved (even if belonging to different domains) 
    #          o Bonds between ATOMS (even if they belonging to different domains)
    # 
    #       We discard Bonds between CG beads that belong to different domains.

    bonds_CG_SameDomain = [ i for i in bonds_mult_res_prot if len(i) == 10 if i[8] == i[9]] 

    bonds_AT1= [ i for i in bonds_mult_res_prot if len(i) == 8]
    bonds_AT2= [ i for i in bonds_mult_res_prot if len(i) == 9]    # If ONLY one CG bead is present 

    bonds_mult_res_prot = bonds_CG_SameDomain + bonds_AT1 + bonds_AT2 

    bonds_mult_res_prot = sorted(bonds_mult_res_prot, key=itemgetter(0,1))

    bonds_mult_res_prot = list(bonds_mult_res_prot for bonds_mult_res_prot,_ in itertools.groupby(bonds_mult_res_prot))

    # A.6.6- Final part: 
    #            If len(list) = 8  => Nothing to do 
    # 	         If len(list) = 9  => Removing i[-1]
    # 		 If len(list) = 10 => Removing i[-1] and i[-2]

    for i in bonds_mult_res_prot:
        if(len(i) == 10):      # Two "cg" beads: removing the last two elements 
            del(i[-1])
            del(i[-1])

    for i in bonds_mult_res_prot:
        if(len(i) ==9):        # One "cg" beads: removing the last element  
            del(i[-1])

else: 
    print('Please, note that No Domain Division File is present. Take care of it.')

          ############################################################
############ (B) PAIRS                                             ####################################
          ############################################################

"""
   Knowledge of the number of pairs and which pairs are present in CANVAS model.
   Pairs involve only 2 atoms with an extra LJ interaction 1-4, therefore we can find it both for atomistic and mult-res protein. 

   We have that: 
  
   For two survived atoms, if in the fully atomistic representation, a pair exists, the same atomistic pair must be reproduced 
   in the Mul-Res representation. However, we exclude all the pairs in which, at least one CG beads are involved. Indeed, we do not want 
   an extra rescaled Coulomb and Lennard Jones in the CG part (only AT-AT pairs are kept).

   However, if the flag "-r/--resc14 y" is set, then we maintein the pair where a CG bead is involved,
   discarding only pairs where both CG beads are involved (AT-AT and CG-AT are kept). 
   In case of charmm forcefield the last operation is not possible, since it not possible to know in advance the values of sigma14 and epsilon14
   for CG beads, because the latter are parametrized, only for the fullyAT representation in ffnonbonded.itp file in "[ pairtypes ]" section. 
   At the contrary, in Amber ff. the [ pairtypes ] section does not exists, as the coulomb for each couple of atoms is rescaled by a factor 0.833 
   (fudgeQQ = 0.833), while the VdW is rescaled by a factor 0.5 (fudgeLJ = 0.5)
""" 


#B.1- After reading the topology file, and which couples of atoms have a pair in the fullyAT representation (original_pairs_protein_fully_at = toop[4]) 
#     we check if two survived atoms keep this property.

 
pairs_mult_res_prot = []

for i in original_pairs_protein_fully_at:
    result=all(x in list_survived_only_num for x in i)
    if(result==True):
        pairs_mult_res_prot.append(i)           # list survived atoms that have a pair (original index of fullyAT representation)


#B.2- Appending, for each couple of survived atoms the "at_type" and "at/cg type" (i[8], and i[7] of list_survived_atoms, respectively)   
#     [a1, a2, at_type-1, at_type-2, at/cg-1, at/cg-2]


pairs_mult_res_prot = [i + [k[8]] for i in pairs_mult_res_prot for k in list_survived_atoms if(k[3] == i[0])]
pairs_mult_res_prot = [i + [k[8]] for i in pairs_mult_res_prot for k in list_survived_atoms if(k[3] == i[1])] 

pairs_mult_res_prot = [i + [k[7]] for i in pairs_mult_res_prot for k in list_survived_atoms if(k[3] == i[0])]
pairs_mult_res_prot = [i + [k[7]] for i in pairs_mult_res_prot for k in list_survived_atoms if(k[3] == i[1])]


#B.3- After reading the pairtype in ffnonbonded file (ffnonbon[4]), we append the value if "sigma14" and "epsilon14" for each couple of 
#     atoms that require a rescaled LJ and Coul. 'pairtypes' and 'pairs_mult_res_prot' look like as follow: 
#
#               pairs_mult_res_prot : [[1,2,X,Y,at,at],[2,3,A,B,at,at],[...]]
#               pairtypes           : [[X,Y,sigma14,eps14],[...]] 


for i in pairs_mult_res_prot:
    for j in pairtypes:
        if((j[0]==i[2] and j[1]==i[3]) or (j[0]==i[3] and j[1]==i[2])):

            i.append(float(j[2]))        # Append sigma14 value ####(in Groamacs units)
            i.append(float(j[3]))        # Append Eps14 value ###(in Gromacs units)

#B.4- Reconverting the original index (of fullyAT representation) in the new index after decimated the most of atoms in CANVAS model
#     and substitute the original value (1st and 2nd column) 

for i in list_survived_atoms:
    for j in pairs_mult_res_prot:
        if(i[3] == j[0]):
            j[0] = i[9]

for i in list_survived_atoms:
    for j in pairs_mult_res_prot:
        if(i[3] == j[1]):
            j[1] = i[9]

#B.5- Discarding all the couples of PAIR in which a bond is present between them. Indeed, if two atoms or beads are bonded, cannot 
#     also have a rescaled Coulomb and LJ term (pairs). 

bonds_only_atnum = set([(elem[0], elem[1]) for elem in bonds_mult_res_prot])   #--> {(a1,a2), (an, am), ...(xxx,xxx)} 

pairs_mult_res_prot = [elem for elem in pairs_mult_res_prot if not (elem[0], elem[1]) in bonds_only_atnum]

#B.6- o  If '-r/--resc14' optional flag is set and the associate string is "Y/y", then 
#        the code manteins from PAIRS all the couples where, at least a CG beads is involved comparing with the all-atom representation
#        i.e.  CG-AT and AT-AT pairs are mantained
#
#     o  Otherwise, if '-r/--resc14' flag is set and the associate string is not Y or y, then the program returns an error. 
#
#     o  Finally, if '-r/--resc14' flag is not set, then the code keeps in PAIRS only AT-AT pairs. All the coouples where, where at least 
#        a CG bead is involved comparing with the all-atom representation is discarded.  

if(Flag_pairtypes == False):   # i.e. the pairtypes list is empty, then the rescaled14 flag is allowed as optional  
 
    if(rescaled14 is not None):
        if(rescaled14 == "y" or rescaled14 == "Y"):
            pairs_mult_res_prot = [x for x in pairs_mult_res_prot if x[4]=="at" or x[5]=="at"]
            print("\nWarning! All pairs where a CG bead are kept if the corresponding all-atom representation that pair is present!")
            print("It could create some artifacts... 49% completed.\n") 
        else:
            print("Error! The accepted arguments of '-r/--resc14'  are 'y' or 'Y' for indicating that all pairs where a CG bead is present in the all-atom representation are kept! Other strings are not allowed.\n")
            print("Look below for further help\n") 
            print_help_main_CANVAS()
            quit() 
    else:
        pairs_mult_res_prot = [x for x in pairs_mult_res_prot if x[4]=="at" and x[5]=="at"]



if(Flag_pairtypes == True):   # i.e. the pairtypes list is NOT-empty, then the rescaled14 flag is not permitted, and only at-at pairs are allowed  
    if(rescaled14 is not None): 
        print("\nWarning! The forcefield used ({}) does not allow 'rescaled14' flag. Only pairs between atoms (AT-AT) are permitted.".format(forcefield))
        print("The code will ignore this flag... 49% completed.\n") 
    else:
       pairs_mult_res_prot = [x for x in pairs_mult_res_prot if x[4]=="at" and x[5]=="at"]

print("\nlist of pairs updated... 50% completed\n")


         ############################################################
############ (C) ANGLES  				          ##########################################
          ###########################################################

'''
   In this part we need to know the number of angles and which angles are between in CANVAS model 
   We have to take in account that the angle among three survived atoms could be only atomistic.

   Basically, only one scenario is possible: 
   
   1. Three survived atoms, in the fully atomistic representation have an angle, 
      therefore the same atomistic angle must be reproduced in the Mul-Res representation 

   The strategy applied is the same of BOND part (A). 
'''

#C.1- After reading the topology file, and which couples of atoms have an angle in the fullyAT representation
#     (original_angles_protein_fully_at = toop[2]), we check if three survived atoms keep this property. 

angles_prot_fully_at = []

for i in original_angles_protein_fully_at:
    result=all(x in list_survived_only_num for x in [i[0],i[1],i[2]])
    if(result==True):
        angles_prot_fully_at.append(i)         # list survived atoms that have an angle (original index of fullyAT representation) 


#C.2- Appending, for each triplet of survived atoms the "at_type" (i[8] of list_survived_atoms)   
#           [a1, a2, a3, at_type-1, at_type-2, at_type-3]

angles_prot_fully_at = [i + [k[8]] for i in angles_prot_fully_at for k in list_survived_atoms if(k[3] == i[0])]
angles_prot_fully_at = [i + [k[8]] for i in angles_prot_fully_at for k in list_survived_atoms if(k[3] == i[1])]
angles_prot_fully_at = [i + [k[8]] for i in angles_prot_fully_at for k in list_survived_atoms if(k[3] == i[2])]


#C.3- After reading the angletypes in ffbonded file (angletypes = ffbon[1]), we append the value of angle elastic constant (cth) 
#     and the angle (th0) for each triplet of survived atoms
#
#     Basically angletypes and angles_prot_fully_at look like as follow: 
#
#               angles_prot_fully_at : [[1,2,3, func, X,Y,Z],[2,3,5, A,B,C],[...]]
#               angletypes           : [[X,Y,Z, func, th0,cth],[...]]                      if func = 1 
#               angletypes           : [[X,Y,Z, func, th0,cth,ub0,kub],[...]]             if func = 5
# 
# 							       where 'func' could be 1 or 5 (usually). 

for i in angles_prot_fully_at:
    for j in angletypes:
        if((j[0]==i[4] and j[1]==i[5] and j[2]==i[6]) or (j[0]==i[6] and j[1]==i[5] and j[2]==i[4])):

            i.append(float(j[5]*8.368))       # Append cth value (recoverted in Gromacs units)
            i.append(float(j[4]))             # Append th0 value (in Gromacs units equals also to Lammps units)
 
            if(i[3]==5):                      # If func=5, then we add also "ub0" and "kub"
                i.append(float(j[6]/10))      # Append ub0 value (recoverted in Gromacs units)
                i.append(float(j[7]*836.8))   # Append kub value (recoverted in Gromacs units)


#C.4- Reconverting the original index (of fullyAT representation) in the new index after decimated the most of atoms in CANVAS model
#     and substitute the original value (1st, 2nd, and 3rd columns) 

for i in list_survived_atoms:
    for j in angles_prot_fully_at:
        if(i[3] == j[0]):
            j[0] = i[9]

for i in list_survived_atoms:
    for j in angles_prot_fully_at:
        if(i[3] == j[1]):
            j[1] = i[9]

for i in list_survived_atoms:
    for j in angles_prot_fully_at:
        if(i[3] == j[2]):
            j[2] = i[9]

#C.5- Adding the CG/AT-type for each bead and atoms. 


angles_prot_fully_at = [j + [i[7]] for j in angles_prot_fully_at for i in list_survived_atoms if(i[9] == j[0])]
angles_prot_fully_at = [j + [i[7]] for j in angles_prot_fully_at for i in list_survived_atoms if(i[9] == j[1])]
angles_prot_fully_at = [j + [i[7]] for j in angles_prot_fully_at for i in list_survived_atoms if(i[9] == j[2])]


#C.6- Adding the new at_type (A{$i}) in case of CG beads. The at_type does not change for AT atoms. 

for i in list_survived_atoms:
    for j in angles_prot_fully_at: 
        if(i[9] == j[0]):
            j[4] = i[16]

for i in list_survived_atoms:
    for j in angles_prot_fully_at:
        if(i[9] == j[1]):
            j[5] = i[16]

for i in list_survived_atoms:
    for j in angles_prot_fully_at:
        if(i[9] == j[2]):
            j[6] = i[16] 

#C.7- Creating a new list called "angles_mult_res_prot" (in analogy with bonds and pairs)  which is simply angles_prot_fully_at, 
#     and sorting the list. 


angles_mult_res_prot = angles_prot_fully_at
angles_mult_res_prot = sorted(angles_mult_res_prot, key=itemgetter(0,1,2))  # sorting 1st, then 2nd and finally the 3rd column.

print("\nlist of angles updated... 55% completed\n")

         ############################################################
############ (D) DIHEDRALS  					  ########################################
          ###########################################################

"""
   In this part we need to know the number of dihedrals and which dihedrals are between the Multiple resolution protein. 
   We have to take in account that the dihedral between four survived atoms could be ONLY atomistic, like for angles. 


   As first, in Gromacs there is the distintion between propers and impropers dihedrals. Indeed, the dihedrals, 
   according with the functional form could be treated in: 

   1) dihedrals with functional form 4 --> dih4 
   2) dihedrals with functional form 0 --> dih9 
   3) special dihedrals (torsion)      --> tors    (they can be found only by using "Amber99sb-ildn.ff" forcefield) 
   4) dihedrals with functional form 2 --> dih2 

   Therefore, because the aforementioned subdivision, this part requires more steps with respect the three previous parts. 
"""


#D.1- After reading the topology file, and which quartet of atoms have an dihedral in the fullyAT representation
#     (original_dihedrals_protein_fully_at = toop[3]), we first split dihedrals according with he functional form, 
#     i.e. 4, 9, torsion, or 2. 


original_dihedrals_protein_fully_at_dih4 = []
original_dihedrals_protein_fully_at_dih9 = []
original_dihedrals_protein_fully_at_tors = []
original_dihedrals_protein_fully_at_dih2 = []

for i in original_dihedrals_protein_fully_at:

    if(len(i) > 5):
        original_dihedrals_protein_fully_at_tors.append(i)

    if(len(i) == 5):
        if(i[4]==4):                                               # i.e. lenght = 5 and funct = 4:
            original_dihedrals_protein_fully_at_dih4.append(i)

        if(i[4]==9):                                               # i.e. lenght = 5 and funct = 9.
            original_dihedrals_protein_fully_at_dih9.append(i)

        if(i[4]==2):						   # i.e. lenght = 5 and funct = 2 
            original_dihedrals_protein_fully_at_dih2.append(i)     


#D.2- Checking if four survived atoms keep one the three previous properties and create three new lists: 
#     dihedrals_prot_fully_at_dih4, dihedrals_prot_fully_at_dih9, dihedrals_prot_fully_at_tors,  dihedrals_prot_fully_at_dih2

dihedrals_prot_fully_at_dih4 = []

for i in original_dihedrals_protein_fully_at_dih4:
    result=all(x in list_survived_only_num for x in [i[0],i[1],i[2],i[3]])
    if(result==True):
        dihedrals_prot_fully_at_dih4.append(i)



dihedrals_prot_fully_at_dih9 = []

for i in original_dihedrals_protein_fully_at_dih9:
    result=all(x in list_survived_only_num for x in [i[0],i[1],i[2],i[3]])
    if(result==True):
        dihedrals_prot_fully_at_dih9.append(i)



dihedrals_prot_fully_at_tors = []

for i in original_dihedrals_protein_fully_at_tors:
    result=all(x in list_survived_only_num for x in [i[0],i[1],i[2],i[3]])
    if(result==True):
        dihedrals_prot_fully_at_tors.append(i)


dihedrals_prot_fully_at_dih2 = []

for i in original_dihedrals_protein_fully_at_dih2:
    result=all(x in list_survived_only_num for x in [i[0],i[1],i[2],i[3]])
    if(result==True):
        dihedrals_prot_fully_at_dih2.append(i)


#D.3- Appending, for each quartet of survived atoms the "at_type" (i[8] of list_survived_atoms)   
#           [a1, a2, a3, a4, at_type-1, at_type-2, at_type-3, at_type-4]
#
#     This operation is done for dih4, dih9, tors, and dih2.  

dihedrals_prot_fully_at_dih4 = [ i + [k[8]] for i in dihedrals_prot_fully_at_dih4 for k in list_survived_atoms if(k[3] == i[0])]
dihedrals_prot_fully_at_dih4 = [ i + [k[8]] for i in dihedrals_prot_fully_at_dih4 for k in list_survived_atoms if(k[3] == i[1])]
dihedrals_prot_fully_at_dih4 = [ i + [k[8]] for i in dihedrals_prot_fully_at_dih4 for k in list_survived_atoms if(k[3] == i[2])]
dihedrals_prot_fully_at_dih4 = [ i + [k[8]] for i in dihedrals_prot_fully_at_dih4 for k in list_survived_atoms if(k[3] == i[3])]

dihedrals_prot_fully_at_dih9 = [ i + [k[8]] for i in dihedrals_prot_fully_at_dih9 for k in list_survived_atoms if(k[3] == i[0])]
dihedrals_prot_fully_at_dih9 = [ i + [k[8]] for i in dihedrals_prot_fully_at_dih9 for k in list_survived_atoms if(k[3] == i[1])]
dihedrals_prot_fully_at_dih9 = [ i + [k[8]] for i in dihedrals_prot_fully_at_dih9 for k in list_survived_atoms if(k[3] == i[2])]
dihedrals_prot_fully_at_dih9 = [ i + [k[8]] for i in dihedrals_prot_fully_at_dih9 for k in list_survived_atoms if(k[3] == i[3])]

dihedrals_prot_fully_at_tors = [ i + [k[8]] for i in dihedrals_prot_fully_at_tors for k in list_survived_atoms if(k[3] == i[0])]
dihedrals_prot_fully_at_tors = [ i + [k[8]] for i in dihedrals_prot_fully_at_tors for k in list_survived_atoms if(k[3] == i[1])]
dihedrals_prot_fully_at_tors = [ i + [k[8]] for i in dihedrals_prot_fully_at_tors for k in list_survived_atoms if(k[3] == i[2])]
dihedrals_prot_fully_at_tors = [ i + [k[8]] for i in dihedrals_prot_fully_at_tors for k in list_survived_atoms if(k[3] == i[3])]

dihedrals_prot_fully_at_dih2 = [ i + [k[8]] for i in dihedrals_prot_fully_at_dih2 for k in list_survived_atoms if(k[3] == i[0])]
dihedrals_prot_fully_at_dih2 = [ i + [k[8]] for i in dihedrals_prot_fully_at_dih2 for k in list_survived_atoms if(k[3] == i[1])]
dihedrals_prot_fully_at_dih2 = [ i + [k[8]] for i in dihedrals_prot_fully_at_dih2 for k in list_survived_atoms if(k[3] == i[2])]
dihedrals_prot_fully_at_dih2 = [ i + [k[8]] for i in dihedrals_prot_fully_at_dih2 for k in list_survived_atoms if(k[3] == i[3])]

#D.4- After reading the dihedraltypes in ffbonded file,
#     (dihedraltypes_dih4 = ffbon[2], dihedraltypes_dih9 = ffbon[3], dihedraltypes_tors = ffbon[4],  dihedraltypes_dih2 = ffbon[5])
#     we append the value of the dihedral elastic constant (k_phi), the phase (phi0) and the multiplicity (n)
#     for each quartet of survived atoms whose funtional form is 4 (dih4) [later we treat dih9, tors, and dih2]
#
#     'dihedraltypes_dih4' and 'dihedral_prot_fully_at_dih4' look like as follow: 
#
#               dihedrals_prot_fully_at_dih4  : [[1,2,3,4, funct, A,N,M,V],[2,3,5,10, funct, A,B,C,R],[...]]
#               dihedraltypes_dih4            : [[A,B,C,D, funct, phi_0,k_phi,n],[...]]  where funct = 4. 
#
#      Attention: At difference with bonds, or angles, some diehdraltypes_dih4 (from ffbonded.itp), in the first 4 columns 
#      are indicated with the letter "X". "X" is not a specific atomtype, but it means that it may represent 
#      whatever atomtype. Therefore, as first, we compare the two list of lists in case they find precisely those 4 atomtypes (avoiding X). 
#      Contrarily (it is not so rare), we need to take in account that those dihedraltypes are signed with X. 

#      For each quartet we can find once or twice X. Therefore we may encounter the following possibilities: 

#       1 --> [X, b, c, d] or [d, c, b, X]               "or" indicates its symmetric     
#       2 --> [a, X, c, d] or [d, c, X, a]     
#       3 --> [a, b, X, d] or [d, X, b, a]    
#       4 --> [a, b, c, X] or [X, c, b, a]

#       1'--> [X X c d] or [d c X X]
#       2'--> [X b X d] or [d X b X]
#       3'--> [X b c X] or [X c b X]
#       4'--> [a X X d] or [d X X a]
#       5'--> [a X c X] or [X c X a]
#       6'--> [a b X X] or [X X b a]

for i in dihedrals_prot_fully_at_dih4:
    flag_dih4 = True
    for j in dihedraltypes_dih4:

        if((j[0]==i[5] and j[1]==i[6] and j[2]==i[7] and j[3]==i[8]) or (j[0]==i[8] and j[1]==i[7] and j[2]==i[6] and j[3]==i[5])):  # consider symmetric 
            i.append(j[5])           # Append phase value  phi_0 (the same both for Gromacs and Lammps) in deg units.
            i.append(j[6]*4.184)     # Append the dihedral constant energy k_phi (re-converted in Gromacs units)
            i.append(j[7])           # Append the value of the multiplicity n (The same for Gromacs and Lammps).
            flag_dih4 = False

        if(flag_dih4 == True):       # It means that no quartet of atomtypes was found in dihedraltype 4 and give a look to the case in which X is present

            # 4 cases (plus symmetric ones): only one X... 
            if((j[0]=="X" and j[1]==i[6] and j[2]==i[7] and j[3]==i[8]) or (j[0]==i[8] and j[1]==i[7] and j[2]==i[6] and j[3]=="X")): # (1)
                i.append(j[5])
                i.append(j[6]*4.184)
                i.append(j[7])

            if((j[0]==i[5] and j[1]=="X" and j[2]==i[7] and j[3]==i[8]) or (j[0]==i[8] and j[1]==i[7] and j[2]=="X" and j[3]==i[5])): # (2)   
                i.append(j[5])
                i.append(j[6]*4.184)
                i.append(j[7])

            if((j[0]==i[5] and j[1]==i[6] and j[2]=="X" and j[3]==i[8]) or (j[0]==i[8] and j[1]=="X" and j[2]==i[6] and j[3]==i[5])): # (3)
                i.append(j[5])
                i.append(j[6]*4.184)
                i.append(j[7])

            if((j[0]==i[5] and j[1]==i[6] and j[2]==i[7] and j[3]=="X") or (j[0]=="X" and j[1]==i[7] and j[2]==i[6] and j[3]==i[5])): # (4)           
                i.append(j[5])
                i.append(j[6]*4.184)
                i.append(j[7])

            # 6 cases (plus symmetric ones): twice X.
            if((j[0]=="X" and j[1]=="X" and j[2]==i[7] and j[3]==i[8]) or (j[0]==i[8] and j[1]==i[7] and j[2]=="X" and j[3]=="X")):   # (1')
                i.append(j[5])
                i.append(j[6]*4.184)
                i.append(j[7])

            if((j[0]=="X" and j[1]==i[6] and j[2]=="X" and j[3]==i[8]) or (j[0]==i[8] and j[1]=="X" and j[2]==i[6] and j[3]=="X")):   # (2')
                i.append(j[5])
                i.append(j[6]*4.184)
                i.append(j[7])

            if((j[0]=="X" and j[1]==i[6] and j[2]==i[7] and j[3]=="X") or (j[0]=="X" and j[1]==i[7] and j[2]==i[6] and j[3]=="X")):   # (3')
                i.append(j[5])
                i.append(j[6]*4.184)
                i.append(j[7])

            if((j[0]==i[5] and j[1]=="X" and j[2]=="X" and j[3]==i[8]) or (j[0]==i[8] and j[1]=="X" and j[2]=="X" and j[3]==i[5])):   # (4')
                i.append(j[5])
                i.append(j[6]*4.184)
                i.append(j[7])

            if((j[0]==i[5] and j[1]=="X" and j[2]==i[7] and j[3]=="X") or (j[0]=="X" and j[1]==i[7] and j[2]=="X" and j[3]==i[5])):   # (5')
                i.append(j[5])
                i.append(j[6]*4.184)
                i.append(j[7])

            if((j[0]==i[5] and j[1]==i[6] and j[2]=="X" and j[3]=="X") or (j[0]=="X" and j[1]=="X" and j[2]==i[6] and j[3]==i[5])):   # (6')
                i.append(j[5])


#D.5- Sometimes (and it is not so rare) (it should not happens for dih4, but only for dih9 happens), some quartets of atoms are repeated more than once 
#      according with the multiplicity (The energy and the phase change, as well). For example in dih9 the quartet of atomtypes C N CT CT in ffbonded.itp 
#      is written 3 times (n = 1, 2, 3). Therefore according with the previous calculation the dihedral 1230 1232 1234 126 (for lysozyme, for example) 
#      has to be splitted in 3, based on its multiplicity.

#      dihedrals_prot_fully_at_dih4 = [a1,a2,a3,a4, funct, attype1,attype2,attype3,attype4, phi0_1, kphi_1,n1, phi0_2, kphi_2, n2, phi0_3, kphi_3, n3]

#      The latter, after splitting, becomes: 

#      [at1, at2, at3, at4, funct, attype1, attype2, attype3, attype4, phi0_1, kphi_1, n1]
#      [at1, at2, at3, at4, funct, attype1, attype2, attype3, attype4, phi0_2, kphi_2, n2]
#      [at1, at2, at3, at4, funct, attype1, attype2, attype3, attype4, phi0_3, kphi_3, n3]

#      It shows that:
#
#      if mult = 1  => lenght of internal list of dihedrals_prot_fully_at_dih4 = 12
#      if mult = 2  => lenght of internal list of dihedrals_prot_fully_at_dih4 = 15
#      if mult = 3  => lenght of internal list of dihedrals_prot_fully_at_dih4 = 18
#      if mult = 4  => lenght of internal list of dihedrals_prot_fully_at_dih4 = 21    (mult=4 should be the maximum possible value of multiplicity)
#      if mult = 5  => lenght of internal list of dihedrals_prot_fully_at_dih4 = 24    (but, for safety, we split lists until mult=6) 
#      if mult = 6  => lenght of internal list of dihedrals_prot_fully_at_dih4 = 27

for i in dihedrals_prot_fully_at_dih4:

    if(len(i) == 15):
        d1 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[9],i[10],i[11]]
        d2 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[12],i[13],i[14]]

        dihedrals_prot_fully_at_dih4.append(d1)
        dihedrals_prot_fully_at_dih4.append(d2)
 

    if(len(i) == 18):
        d1 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[9],i[10],i[11]]
        d2 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[12],i[13],i[14]]
        d3 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[15],i[16],i[17]]

        dihedrals_prot_fully_at_dih4.append(d1)
        dihedrals_prot_fully_at_dih4.append(d2)
        dihedrals_prot_fully_at_dih4.append(d3)


    if(len(i) == 21):
        d1 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[9],i[10],i[11]]
        d2 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[12],i[13],i[14]]
        d3 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[15],i[16],i[17]]
        d4 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[18],i[19],i[20]]

        dihedrals_prot_fully_at_dih4.append(d1)
        dihedrals_prot_fully_at_dih4.append(d2)
        dihedrals_prot_fully_at_dih4.append(d3)
        dihedrals_prot_fully_at_dih4.append(d4)


    if(len(i) == 24):
        d1 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[9],i[10],i[11]]
        d2 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[12],i[13],i[14]]
        d3 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[15],i[16],i[17]]
        d4 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[18],i[19],i[20]]
        d5 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[21],i[22],i[23]]

        dihedrals_prot_fully_at_dih4.append(d1)
        dihedrals_prot_fully_at_dih4.append(d2)
        dihedrals_prot_fully_at_dih4.append(d3)
        dihedrals_prot_fully_at_dih4.append(d4)
        dihedrals_prot_fully_at_dih4.append(d5)


    if(len(i) == 27):
        d1 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[9],i[10],i[11]]
        d2 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[12],i[13],i[14]]
        d3 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[15],i[16],i[17]]
        d4 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[18],i[19],i[20]]
        d5 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[21],i[22],i[23]]
        d6 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[24],i[25],i[26]]

        dihedrals_prot_fully_at_dih4.append(d1)
        dihedrals_prot_fully_at_dih4.append(d2)
        dihedrals_prot_fully_at_dih4.append(d3)
        dihedrals_prot_fully_at_dih4.append(d4)
        dihedrals_prot_fully_at_dih4.append(d5)
        dihedrals_prot_fully_at_dih4.append(d6)



dihedrals_prot_fully_at_dih4 = [l for l in dihedrals_prot_fully_at_dih4 if len(l)==12]             # keep lists with lenght = 12,  after splitting

dihedrals_prot_fully_at_dih4 = sorted(dihedrals_prot_fully_at_dih4, key=itemgetter(0,1,2,3,11))    # Sorting the 1st column (at1), then the 2nd (at2)
                                                                                                   #, the 3rd (at3), the 4th (at4) and 12th (n) 

#D.6- Repeating the previous four steps (D.4, and D.5) for dih9 

for i in dihedrals_prot_fully_at_dih9:
    flag_dih9 = True
    for j in dihedraltypes_dih9:

        if((j[0]==i[5] and j[1]==i[6] and j[2]==i[7] and j[3]==i[8]) or (j[0]==i[8] and j[1]==i[7] and j[2]==i[6] and j[3]==i[5])):  # consider symmetric 
            i.append(j[5])           # Append phase value  phi_0 (the same both for Gromacs and Lammps) in deg units.
            i.append(j[6]*4.184)     # Append the dihedral constant energy k_phi (re-converted in Gromacs units)
            i.append(j[7])           # Append the value of the multiplicity n (The same for Gromacs and Lammps).
            flag_dih9 = False

        if(flag_dih9 == True):       # It means that no quartet of atomtypes was found in dihedraltype 4 and give a look to the case in which X is present

            # 4 cases (plus symmetric ones): only one X... 
            if((j[0]=="X" and j[1]==i[6] and j[2]==i[7] and j[3]==i[8]) or (j[0]==i[8] and j[1]==i[7] and j[2]==i[6] and j[3]=="X")): # (1)
                i.append(j[5])
                i.append(j[6]*4.184)
                i.append(j[7])

            if((j[0]==i[5] and j[1]=="X" and j[2]==i[7] and j[3]==i[8]) or (j[0]==i[8] and j[1]==i[7] and j[2]=="X" and j[3]==i[5])): # (2)   
                i.append(j[5])
                i.append(j[6]*4.184)
                i.append(j[7])

            if((j[0]==i[5] and j[1]==i[6] and j[2]=="X" and j[3]==i[8]) or (j[0]==i[8] and j[1]=="X" and j[2]==i[6] and j[3]==i[5])): # (3)
                i.append(j[5])
                i.append(j[6]*4.184)
                i.append(j[7])

            if((j[0]==i[5] and j[1]==i[6] and j[2]==i[7] and j[3]=="X") or (j[0]=="X" and j[1]==i[7] and j[2]==i[6] and j[3]==i[5])): # (4)           
                i.append(j[5])
                i.append(j[6]*4.184)
                i.append(j[7])

            # 6 cases (plus symmetric ones): twice X.
            if((j[0]=="X" and j[1]=="X" and j[2]==i[7] and j[3]==i[8]) or (j[0]==i[8] and j[1]==i[7] and j[2]=="X" and j[3]=="X")):   # (1')
                i.append(j[5])
                i.append(j[6]*4.184)
                i.append(j[7])

            if((j[0]=="X" and j[1]==i[6] and j[2]=="X" and j[3]==i[8]) or (j[0]==i[8] and j[1]=="X" and j[2]==i[6] and j[3]=="X")):   # (2')
                i.append(j[5])
                i.append(j[6]*4.184)
                i.append(j[7])

            if((j[0]=="X" and j[1]==i[6] and j[2]==i[7] and j[3]=="X") or (j[0]=="X" and j[1]==i[7] and j[2]==i[6] and j[3]=="X")):   # (3')
                i.append(j[5])
                i.append(j[6]*4.184)
                i.append(j[7])

            if((j[0]==i[5] and j[1]=="X" and j[2]=="X" and j[3]==i[8]) or (j[0]==i[8] and j[1]=="X" and j[2]=="X" and j[3]==i[5])):   # (4')
                i.append(j[5])
                i.append(j[6]*4.184)
                i.append(j[7])

            if((j[0]==i[5] and j[1]=="X" and j[2]==i[7] and j[3]=="X") or (j[0]=="X" and j[1]==i[7] and j[2]=="X" and j[3]==i[5])):   # (5')
                i.append(j[5])
                i.append(j[6]*4.184)
                i.append(j[7])

            if((j[0]==i[5] and j[1]==i[6] and j[2]=="X" and j[3]=="X") or (j[0]=="X" and j[1]=="X" and j[2]==i[6] and j[3]==i[5])):   # (6')
                i.append(j[5])
                i.append(j[6]*4.184)
                i.append(j[7])



for i in dihedrals_prot_fully_at_dih9:

    if(len(i) == 15):
        d1 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[9],i[10],i[11]]
        d2 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[12],i[13],i[14]]

        dihedrals_prot_fully_at_dih9.append(d1)
        dihedrals_prot_fully_at_dih9.append(d2)
 

    if(len(i) == 18):
        d1 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[9],i[10],i[11]]
        d2 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[12],i[13],i[14]]
        d3 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[15],i[16],i[17]]

        dihedrals_prot_fully_at_dih9.append(d1)
        dihedrals_prot_fully_at_dih9.append(d2)
        dihedrals_prot_fully_at_dih9.append(d3)


    if(len(i) == 21):
        d1 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[9],i[10],i[11]]
        d2 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[12],i[13],i[14]]
        d3 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[15],i[16],i[17]]
        d4 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[18],i[19],i[20]]

        dihedrals_prot_fully_at_dih9.append(d1)
        dihedrals_prot_fully_at_dih9.append(d2)
        dihedrals_prot_fully_at_dih9.append(d3)
        dihedrals_prot_fully_at_dih9.append(d4)


    if(len(i) == 24):
        d1 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[9],i[10],i[11]]
        d2 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[12],i[13],i[14]]
        d3 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[15],i[16],i[17]]
        d4 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[18],i[19],i[20]]
        d5 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[21],i[22],i[23]]

        dihedrals_prot_fully_at_dih9.append(d1)
        dihedrals_prot_fully_at_dih9.append(d2)
        dihedrals_prot_fully_at_dih9.append(d3)
        dihedrals_prot_fully_at_dih9.append(d4)
        dihedrals_prot_fully_at_dih9.append(d5)


    if(len(i) == 27):
        d1 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[9],i[10],i[11]]
        d2 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[12],i[13],i[14]]
        d3 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[15],i[16],i[17]]
        d4 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[18],i[19],i[20]]
        d5 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[21],i[22],i[23]]
        d6 = [i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[24],i[25],i[26]]

        dihedrals_prot_fully_at_dih9.append(d1)
        dihedrals_prot_fully_at_dih9.append(d2)
        dihedrals_prot_fully_at_dih9.append(d3)
        dihedrals_prot_fully_at_dih9.append(d4)
        dihedrals_prot_fully_at_dih9.append(d5)
        dihedrals_prot_fully_at_dih9.append(d6)


dihedrals_prot_fully_at_dih9 = [l for l in dihedrals_prot_fully_at_dih9 if len(l)==12]             # keep lists with lenght = 12, after splitting

dihedrals_prot_fully_at_dih9 = sorted(dihedrals_prot_fully_at_dih9, key=itemgetter(0,1,2,3,11))    # Sorting the 1st column (at1), then the 2nd (at2)
                                                                                                   #, the 3rd (at3), the 4th (at4) and 12th (n) 

#D.7- After reading the special dihedraltypes defined as torsion (dihedraltypes_tors = ffbon[4]), 
#     we append the value of the dihedral elastic constant (k_phi), the phase (phi0) and the multiplicity (n)
#     for each quartet of "torsion" survived atoms.
#
#     Basically, "dihedraltypes_tors" and "dihedrals_prot_fully_at_tors" look like as follows: 
#
#        dihedrals_prot_fully_at_tors:  [[a1,a2,a3,a4, funct, string_tors, A,B,C,D], [b1,b2,b3,b4, funct, C,Y,Z,K], [...]]
#        dihedraltypes_tors          :  [[string_tors, phi_0, k_phi, n],[....],[....]]
#
#     This case of much easier than dih4 and dih9 because each string is unique, and the problem of "X" atomtype is not present.
#

for i in dihedrals_prot_fully_at_tors:
    for j in dihedraltypes_tors:
        if(j[0]==i[5]):
            i.append(j[1])           # Append the phase phi0
            i.append(j[2]*4.184)     # Append the energy k_phi (re-converted for Gromacs units) 
            i.append(j[3])           # Append the multiplicity n 

dihedrals_prot_fully_at_tors = sorted(dihedrals_prot_fully_at_tors, key=itemgetter(0,1,2,3,12))

#D.8- About dih2, the situation is analogue to dih4 and dih9, but dih2 do not have multiplicity (n) that 
#     is present in dih4 and dih9. Therefore we repeat only *D.4* section, keeping attention that j[7] of 
#     "dihedraltypes_dih2" (multiplicity) is absent here. The *D.5* section is unnecessary because no quartet of atoms 
#     can be repeated, thus no list-splitting according with the multiplicity must be performed. We have that:  
#
#     'dihedraltypes_dih2' and 'dihedral_prot_fully_at_dih2' look like as follow: 
#
#               dihedrals_prot_fully_at_dih2  : [[1,2,3,4, funct, A,N,M,V],[2,3,5,10, funct, A,B,C,R],[...]]
#               dihedraltypes_dih2            : [[A,B,C,D, funct, phi_0,k_phi],[...]]  where funct = 2. 

for i in dihedrals_prot_fully_at_dih2:
    flag_dih2 = True
    for j in dihedraltypes_dih2:

        if((j[0]==i[5] and j[1]==i[6] and j[2]==i[7] and j[3]==i[8]) or (j[0]==i[8] and j[1]==i[7] and j[2]==i[6] and j[3]==i[5])):  # consider symmetric 
            i.append(j[5])           # Append phase value  phi_0 (the same both for Gromacs and Lammps) in deg units.
            i.append(j[6]*8.368)     # Append the dihedral constant energy k_phi (re-converted in Gromacs units, incl. 1/2 factor)
            flag_dih2 = False

        if(flag_dih2 == True):       # It means that no quartet of atomtypes was found in dihedraltype 4 and give a look to the case in which X is present

            # 4 cases (plus symmetric ones): only one X... 
            if((j[0]=="X" and j[1]==i[6] and j[2]==i[7] and j[3]==i[8]) or (j[0]==i[8] and j[1]==i[7] and j[2]==i[6] and j[3]=="X")): # (1)
                i.append(j[5])
                i.append(j[6]*8.368)

            if((j[0]==i[5] and j[1]=="X" and j[2]==i[7] and j[3]==i[8]) or (j[0]==i[8] and j[1]==i[7] and j[2]=="X" and j[3]==i[5])): # (2)   
                i.append(j[5])
                i.append(j[6]*8.368)

            if((j[0]==i[5] and j[1]==i[6] and j[2]=="X" and j[3]==i[8]) or (j[0]==i[8] and j[1]=="X" and j[2]==i[6] and j[3]==i[5])): # (3)
                i.append(j[5])
                i.append(j[6]*8.368)

            if((j[0]==i[5] and j[1]==i[6] and j[2]==i[7] and j[3]=="X") or (j[0]=="X" and j[1]==i[7] and j[2]==i[6] and j[3]==i[5])): # (4)           
                i.append(j[5])
                i.append(j[6]*8.368)

            # 6 cases (plus symmetric ones): twice X.
            if((j[0]=="X" and j[1]=="X" and j[2]==i[7] and j[3]==i[8]) or (j[0]==i[8] and j[1]==i[7] and j[2]=="X" and j[3]=="X")):   # (1')
                i.append(j[5])
                i.append(j[6]*8.368)

            if((j[0]=="X" and j[1]==i[6] and j[2]=="X" and j[3]==i[8]) or (j[0]==i[8] and j[1]=="X" and j[2]==i[6] and j[3]=="X")):   # (2')
                i.append(j[5])
                i.append(j[6]*8.368)

            if((j[0]=="X" and j[1]==i[6] and j[2]==i[7] and j[3]=="X") or (j[0]=="X" and j[1]==i[7] and j[2]==i[6] and j[3]=="X")):   # (3')
                i.append(j[5])
                i.append(j[6]*8.368)

            if((j[0]==i[5] and j[1]=="X" and j[2]=="X" and j[3]==i[8]) or (j[0]==i[8] and j[1]=="X" and j[2]=="X" and j[3]==i[5])):   # (4')
                i.append(j[5])
                i.append(j[6]*8.368)

            if((j[0]==i[5] and j[1]=="X" and j[2]==i[7] and j[3]=="X") or (j[0]=="X" and j[1]==i[7] and j[2]=="X" and j[3]==i[5])):   # (5')
                i.append(j[5])
                i.append(j[6]*8.368)

            if((j[0]==i[5] and j[1]==i[6] and j[2]=="X" and j[3]=="X") or (j[0]=="X" and j[1]=="X" and j[2]==i[6] and j[3]==i[5])):   # (6')
                i.append(j[5])
                i.append(j[6]*8.368)


dihedrals_prot_fully_at_dih2 = sorted(dihedrals_prot_fully_at_dih2, key=itemgetter(0,1,2,3))       # Sorting the 1st column (at1), then the 2nd (at2)
                                                                                                   #, the 3rd (at3), the 4th (at4)

#D.9- Creating the list dihedrals_prot_fully_at which is the sum of the three previous lists (dih4, dih9, tors & dih2) 

dihedrals_prot_fully_at = dihedrals_prot_fully_at_dih4 + dihedrals_prot_fully_at_dih9 + dihedrals_prot_fully_at_tors + dihedrals_prot_fully_at_dih2

#D.10- Reconverting the original index (of fullyAT representation) in the new index after decimated the most of atoms in CANVAS model
#      and substitute the original value (1st, 2nd, and 3rd, 4th columns) and sort.

for i in list_survived_atoms:
    for j in dihedrals_prot_fully_at:
        if(i[3] == j[0]):
            j[0] = i[9]

for i in list_survived_atoms:
    for j in dihedrals_prot_fully_at:
        if(i[3] == j[1]):
            j[1] = i[9]

for i in list_survived_atoms:
    for j in dihedrals_prot_fully_at:
        if(i[3] == j[2]):
            j[2] = i[9]

for i in list_survived_atoms:
    for j in dihedrals_prot_fully_at:
        if(i[3] == j[3]):
            j[3] = i[9]


dihedrals_prot_fully_at = sorted(dihedrals_prot_fully_at, key=itemgetter(0,1,2,3,-1))   # Sorting a1, a2, a3, a4, and finally the multiplicity n.  
											# Note that in case of dih2, "-1" is kb, not multiplicity 
											# but it is not a problem since we sort first in terms of 
											# a1, a2, a3, a4; dih2 has no "n", therefore there is no 
										        # possibility to have identical a1,a2,a3,a4. Sorting based 
											# on "-1" is done only for identical a1,a2,a3,a4 for dih4 and 
											# dih9 having different multiplicity.  

#D.11- Adding the AT/CG type as for angletypes (C.5)

dihedrals_prot_fully_at = [j + [i[7]] for j in dihedrals_prot_fully_at for i in list_survived_atoms if(i[9] == j[0])]
dihedrals_prot_fully_at = [j + [i[7]] for j in dihedrals_prot_fully_at for i in list_survived_atoms if(i[9] == j[1])]
dihedrals_prot_fully_at = [j + [i[7]] for j in dihedrals_prot_fully_at for i in list_survived_atoms if(i[9] == j[2])]
dihedrals_prot_fully_at = [j + [i[7]] for j in dihedrals_prot_fully_at for i in list_survived_atoms if(i[9] == j[3])]


#D.12- Adding the new at_type (A{$i}) in case of CG beads. The at_type does not change for AT atoms. 

for i in list_survived_atoms:
    for j in dihedrals_prot_fully_at:
        if(i[9] == j[0]):

            if(not j[5].startswith("torsion")):    # in case of torsion    : the atomtypes are j[6], j[7], j[8], and j[9]
                j[5] = i[16]                       # in case of not-torsion: the atomtypes are j[5], j[6], j[7], and j[8]
            else:
                j[6] = i[16] 

for i in list_survived_atoms:
    for j in dihedrals_prot_fully_at:
        if(i[9] == j[1]):

            if(not j[5].startswith("torsion")):
                j[6] = i[16]
            else:
                j[7] = i[16]


for i in list_survived_atoms:
    for j in dihedrals_prot_fully_at:
        if(i[9] == j[2]):

            if(not j[5].startswith("torsion")):
                j[7] = i[16]
            else:
                j[8] = i[16]


for i in list_survived_atoms:
    for j in dihedrals_prot_fully_at:
        if(i[9] == j[3]):

            if(not j[5].startswith("torsion")):
                j[8] = i[16]
            else:
                j[9] = i[16]


#D.13- Creating a new list called "dihedrals_mult_res_prot" (in analogy with bonds, pairs, and angles) which is simply dihedrals_prot_fully_at, 

dihedrals_mult_res_prot = dihedrals_prot_fully_at

#D.14- When we need the input files for Lammps, dih4, dih9 and torsion require an harmonic series expansion of
#      k_phi in 5 coefficients. The latter expansion changes according with phi0 and the multiplicity n. 

#      In general, when n = 0,1,2,3,4         => we have to use the multi/harmonic style with 5 coefficients 
#                  when n = 5,6 and k_phi =0  => we have to use the multi/harmonic style with the 5 coeffcients set to 0.0 
#                  when n = 5,6 and k_phi!=0  => we have to use charmm style 

#      MULTI/HARMONIC STYLE: <count> multi/harmonic <coeff1> <coeff2> <coeff3> <coeff4> <coeff5>
#      CHARMM STYLE        : <count> charmm <k_phi> <n> <phi0> <0.0>     

#      In the charmm style the last coeffcient is 0.0 corresponding to the weight factor. It must be set to 0.0 always for amber ff. 

#      As follow, there is the summerizing table showing the 5 coefficients for the multi/harmonic style (for simplicity we call k_phi as K)

#       #################################################
#       #   n    phi0  |  c1    c2    c3    c4    c5    #
#       # --------------------------------------------  #
#       #   0     0    |  2k     0     0     0     0    #
#       #   0   180    |   0     0     0     0     0    #
#       #                                               # 
#       #   1     0    |   k     k     0     0     0    #
#       #   1   180    |   k    -k     0     0     0    #
#       #                                               # 
#       #   2     0    |   0     0    2k     0     0    #
#       #   2   180    |  2k     0   -2k     0     0    #
#       #                                               #
#       #   3     0    |   k    -3k    0    4k     0    #
#       #   3   180    |   k     3k    0   -4k     0    #
#       #                                               #
#       #   4     0    |  2k     0   -8k     0    8k    #
#       #   4   180    |   0     0    8k     0   -8k    #
#       #                                               #
#       #   5     0    |   0     0     0     0     0    #
#       #                                               #
#       #   6     0    |   0     0     0     0     0    #
#       #################################################  


#       On the other hand, in case of dih2, the functional form is different i.e. k(X-Xi)^2   where X is the greek letter CHI. 
#       Therefore, in Lammps the dihedrals dih2 are defined as IMPROPERS, with the HARMONIC STYLE i.e.: 

#       HARMONIC: <count> harmonic <k_phi> <phi0>    (it is like bonds or angles of type 1) 
#

for i in dihedrals_mult_res_prot:

    if(i[5].startswith("torsion")):                  # We are dealing with torsion => (phi0 --> i[10], kphi --> i[11], n --> i[12])

        if(i[12] ==0 and i[10] ==0.0):               # n=0 & phi0=0
            i.append("multi")
            i.append(2*i[11]/4.184)                  # i[11] is the energy column, Lammps units 
            i.append(0.0)
            i.append(0.0)
            i.append(0.0)                            # We add the string MULTI before and after the 5 coefficients for later easier use. 
            i.append(0.0)
            i.append("multi")

        if(i[12] ==0 and i[10] ==180.0):             # n=0 & phi0=180
            i.append("multi")
            i.append(0.0)
            i.append(0.0)
            i.append(0.0)
            i.append(0.0)
            i.append(0.0)
            i.append("multi")



        if(i[12] ==1 and i[10] ==0.0):             # n=1 & phi0 = 0
            i.append("multi")
            i.append(i[11]/4.184)
            i.append(i[11]/4.184)
            i.append(0.0)
            i.append(0.0)
            i.append(0.0)
            i.append("multi")

        if(i[12] ==1 and i[10] ==180.0):             # n=1 & phi0 = 180.0
            i.append("multi")
            i.append(i[11]/4.184)
            i.append(-1*i[11]/4.184)
            i.append(0.0)
            i.append(0.0)
            i.append(0.0)
            i.append("multi")


        if(i[12] ==2 and i[10] ==0.0):             # n=2 & phi0 = 0
            i.append("multi")
            i.append(0.0)
            i.append(0.0)
            i.append(2*i[11]/4.184)
            i.append(0.0)
            i.append(0.0)
            i.append("multi")

        if(i[12] ==2 and i[10] ==180.0):             # n=2 & phi0 = 180.0
            i.append("multi")
            i.append(2*i[11]/4.184)
            i.append(0.0)
            i.append(-2*i[11]/4.184)
            i.append(0.0)
            i.append(0.0)
            i.append("multi")



        if(i[12] ==3 and i[10] ==0.0):             # n=3 & phi0 = 0
            i.append("multi")
            i.append(i[11]/4.184)
            i.append(-3*i[11]/4.184)
            i.append(0.0)
            i.append(4*i[11]/4.184)
            i.append(0.0)
            i.append("multi")

        if(i[12] ==3 and i[10] ==180.0):             # n=3 & phi0 = 180.0
            i.append("multi")
            i.append(i[11]/4.184)
            i.append(3*i[11]/4.184)
            i.append(0.0)
            i.append(-4*i[11]/4.184)
            i.append(0.0)
            i.append("multi")



        if(i[12] ==4 and i[10] ==0.0):             # n=4 & phi0 = 0
            i.append("multi")
            i.append(2*i[11]/4.184)
            i.append(0.0)
            i.append(-8*i[11]/4.184)
            i.append(0.0)
            i.append(8*i[11]/4.184)
            i.append("multi")

        if(i[12] ==4 and i[10] ==180.0):             # n=4 & phi0 = 180.0
            i.append("multi")
            i.append(0.0)
            i.append(0.0)
            i.append(8*i[11]/4.184)
            i.append(0.0)
            i.append(8*i[11]/4.184)
            i.append("multi")



        if(i[12] ==5 and i[11] ==0.0):             # n =5 & k_phi = 0.0 --> always multi/harmonic 
            i.append("multi")
            i.append(0.0)
            i.append(0.0)
            i.append(0.0)
            i.append(0.0)
            i.append(0.0)
            i.append("multi")

        if(i[12] ==5 and i[11]!=0.0):               # se n=5 & k_phi !=0 --> charmm
            i.append("charmm")
            i.append(i[11]/4.184)                         # add k_phi
            i.append(i[12])                         # add n 
            i.append(i[10])                         # add phi0
            i.append(0.0)                           # add weight factor always equals 0
            i.append("charmm")



        if(i[12] ==6 and i[11] ==0.0):             # n =6 & k_phi = 0.0 --> always multi/harmonic 
            i.append("multi")
            i.append(0.0)
            i.append(0.0)
            i.append(0.0)
            i.append(0.0)
            i.append(0.0)
            i.append("multi")

        if(i[12] ==6 and i[11]!=0.0):               # se n=6 & k_phi !=0 --> charmm
            i.append("charmm")
            i.append(i[11]/4.184)                         # add k_phi
            i.append(i[12])                         # add n 
            i.append(i[10])                         # add phi0
            i.append(0.0)                           # add weight factor always equals 0
            i.append("charmm")


    if(not i[5].startswith("torsion")):            
        if(i[4] != 2):  # We are dealing with dih4 and dih9 => (phi0 --> i[9], kphi --> i[10], n --> i[11]) 

            if(i[11] ==0 and i[9] ==0.0):               # n=0 & phi0=0
                i.append("multi")
                i.append(2*i[10]/4.184)                       #i[10] is the energy column
                i.append(0.0)
                i.append(0.0)
                i.append(0.0)                           # We add the string MULTI before and after the coefficients for later easier use.  
                i.append(0.0)
                i.append("multi")

            if(i[11] ==0 and i[9] ==180.0):             # n=0 & phi0=180
                i.append("multi")
                i.append(0.0)
                i.append(0.0)
                i.append(0.0)
                i.append(0.0)
                i.append(0.0)
                i.append("multi")


 
            if(i[11] ==1 and i[9] ==0.0):              # n=1 & phi0 = 0
                i.append("multi")
                i.append(i[10]/4.184)
                i.append(i[10]/4.184)
                i.append(0.0)
                i.append(0.0)
                i.append(0.0)
                i.append("multi")

            if(i[11] ==1 and i[9] ==180.0):            # n=1 & phi0 = 180.0
                i.append("multi")
                i.append(i[10]/4.184)
                i.append(-1*i[10]/4.184)
                i.append(0.0)
                i.append(0.0)
                i.append(0.0)
                i.append("multi")



            if(i[11] ==2 and i[9] ==0.0):             # n=2 & phi0 = 0
                i.append("multi")
                i.append(0.0)
                i.append(0.0)
                i.append(2*i[10]/4.184)
                i.append(0.0)
                i.append(0.0)
                i.append("multi")

            if(i[11] ==2 and i[9] ==180.0):             # n=2 & phi0 = 180.0
                i.append("multi")
                i.append(2*i[10]/4.184)
                i.append(0.0)
                i.append(-2*i[10]/4.184)
                i.append(0.0)
                i.append(0.0)
                i.append("multi")



            if(i[11] ==3 and i[9] ==0.0):             # n=3 & phi0 = 0
                i.append("multi")
                i.append(i[10]/4.184)
                i.append(-3*i[10]/4.184)
                i.append(0.0)
                i.append(4*i[10]/4.184)
                i.append(0.0)
                i.append("multi")

            if(i[11] ==3 and i[9] ==180.0):             # n=3 & phi0 = 180.0
                i.append("multi")
                i.append(i[10]/4.184)
                i.append(3*i[10]/4.184)
                i.append(0.0)
                i.append(-4*i[10]/4.184)
                i.append(0.0)
                i.append("multi")



            if(i[11] ==4 and i[9] ==0.0):             # n=4 & phi0 = 0
                i.append("multi")
                i.append(2*i[10]/4.184)
                i.append(0.0)
                i.append(-8*i[10]/4.184)
                i.append(0.0)
                i.append(8*i[10]/4.184)
                i.append("multi")

            if(i[11] ==4 and i[9] ==180.0):             # n=4 & phi0 = 180.0
                i.append("multi")
                i.append(0.0)
                i.append(0.0)
                i.append(8*i[10]/4.184)
                i.append(0.0)
                i.append(8*i[10]/4.184)
                i.append("multi")



            if(i[11] ==5 and i[10] ==0.0):             # n =5 & k_phi = 0.0 --> always multi/harmonic 
                i.append("multi")
                i.append(0.0)
                i.append(0.0)
                i.append(0.0)
                i.append(0.0)
                i.append(0.0)
                i.append("multi")

            if(i[11] ==5 and i[10]!=0.0):               # se n=5 & k_phi !=0 --> charmm
                i.append("charmm")
                i.append(i[10]/4.184)                         # add k_phi
                i.append(i[11])                         # add n 
                i.append(i[9])                          # add phi0
                i.append(0.0)                           # add weight factor always equals 0
                i.append("charmm")



            if(i[11] ==6 and i[10] ==0.0):             # n =6 & k_phi = 0.0 --> always multi/harmonic 
                i.append("multi")
                i.append(0.0)
                i.append(0.0)
                i.append(0.0)
                i.append(0.0)
                i.append(0.0)
                i.append("multi")

            if(i[11] ==6 and i[10]!=0.0):               # se n=6 & k_phi !=0 --> charmm
                i.append("charmm")
                i.append(i[10]/4.184)                         # add k_phi
                i.append(i[11])                         # add n 
                i.append(i[9])                          # add phi0
                i.append(0.0)                           # add weight factor always equals 0
                i.append("charmm")


        else:           # We are dealing with dih2 => (phi0 -->i[9], kphi --> i[10])  [no multiplicity "n" present in dih2]
            continue    # Nothing to do, as dih2 will be defined as harmonic, and there is no series expansion. 

print("\nlist of dihedrals updated... 75% completed\n")


         ####################################################################
############ (E) CMAP (Torsional Correction Map) ONLY FOR CHARMM.ff        ##########################################
          ###################################################################

'''
   Only in Charmm.ff, the all-atom topology file presents also the [ cmap ] section. 
   There are 5 atom names that define the two torsional angles. Atoms 1-4 define "phi", while atoms 2-5 define "psi"
   The corresponding atom types are then matched to the correct CMAP type in the cmap.itp file that contains the correction maps. 

   In this part we need to know the number of cmaps and which cmpas are between in CANVAS model 
   We have to take in account that the cmap among five survived atoms could be only atomistic.

   Basically, only one scenario is possible: 
   
   1. Five survived atoms, in the fully atomistic representation have an Torsional Correction Map, 
      therefore the same atomistic angle must be reproduced in the Mul-Res representation. 

   The strategy applied is the same of BOND part (A) and ANGLE part (B). 
'''
if(Flag_cmaps == True):   # i.e. the "original_cmaps_protein_fully_at" list is not-empty, i.e. the charmm.ff is employed  

    #E.1- After reading the topology file, and which quitets of atoms have  in the fullyAT representation
    #     (original_cmaps_protein_fully_at = toop[7]), we check if five survived atoms keep this property. 

    cmaps_prot_fully_at = []

    for i in original_cmaps_protein_fully_at:
        result=all(x in list_survived_only_num for x in [i[0],i[1],i[2],i[3],i[4]])
        if(result==True):
            cmaps_prot_fully_at.append(i)         # list survived atoms that have one or more cmaps (original index of fullyAT representation) 


    #E.2- Appending, for each quintet of survived atoms the "at_type" (i[8] of list_survived_atoms)   
    #           [a1, a2, a3, a4, a5, at_type-1, at_type-2, at_type-3, at_type4, at_type-5]

    cmaps_prot_fully_at = [i + [k[8]] for i in cmaps_prot_fully_at for k in list_survived_atoms if(k[3] == i[0])]
    cmaps_prot_fully_at = [i + [k[8]] for i in cmaps_prot_fully_at for k in list_survived_atoms if(k[3] == i[1])]
    cmaps_prot_fully_at = [i + [k[8]] for i in cmaps_prot_fully_at for k in list_survived_atoms if(k[3] == i[2])]
    cmaps_prot_fully_at = [i + [k[8]] for i in cmaps_prot_fully_at for k in list_survived_atoms if(k[3] == i[3])]
    cmaps_prot_fully_at = [i + [k[8]] for i in cmaps_prot_fully_at for k in list_survived_atoms if(k[3] == i[4])]

    #E.3- After reading the cmaptypes in "cmap.itp" file (cmaptypes = readcmap(f_cmap)), we append the func, grid1, grid2, and 
    #     the next 58 internal lists consisting of 576 values for the Torsional Correction Map for each quintet of survived atoms. 
    #
    #     Basically cmaptypes and cmaps_prot_fully_at look like as follow: 
    #
    #               cmaps_prot_fully_at : [[1,2,3,4,5, func, X,Y,Z,K,H],[2,3,5,8,9, func, A,B,C,D,E],[...]]
    #               cmaptypes           : [[X,Y,Z,K,H, func, grid1, grid2, string1, string2, .... string58],[...]] 
    #									
    #					  where func=1, grid1=24, grid2=24\, string1=1st line of float-values
    #

    for i in cmaps_prot_fully_at:
        for j in cmaptypes:
            if((j[0]==i[6] and j[1]==i[7] and j[2]==i[8] and j[3]==i[9] and j[4]==i[10]) or \
               (j[0]==i[10] and j[1]==i[9] and j[2]==i[8] and j[3] ==i[7] and j[4]==i[6])):

                i.append(int(j[6]))       # Append grid1, i.e. "24"
                i.append(str(j[7]))       # Append grid2, i.e. "24\"
                for count in range(8,66):     
                    i.append(str(j[count]))  # Append the 58 string-list  


    #E.4- Reconverting the original index (of fullyAT representation) in the new index after decimated the most of atoms in CANVAS model
    #     and substitute the original value (1st, 2nd, 3rd, 4th, and 5th columns) 

    for i in list_survived_atoms:
        for j in cmaps_prot_fully_at:
            if(i[3] == j[0]):
                j[0] = i[9]

    for i in list_survived_atoms:
        for j in cmaps_prot_fully_at:
            if(i[3] == j[1]):
                j[1] = i[9]

    for i in list_survived_atoms:
        for j in cmaps_prot_fully_at:
            if(i[3] == j[2]):
                j[2] = i[9]

    for i in list_survived_atoms:
        for j in cmaps_prot_fully_at:
            if(i[3] == j[3]):
                j[3] = i[9]

    for i in list_survived_atoms:
        for j in cmaps_prot_fully_at:
            if(i[3] == j[4]):
                j[4] = i[9]


    #E.5- Adding the CG/AT-type for each bead and atoms. 

    cmaps_prot_fully_at = [j + [i[7]] for j in cmaps_prot_fully_at for i in list_survived_atoms if(i[9] == j[0])]
    cmaps_prot_fully_at = [j + [i[7]] for j in cmaps_prot_fully_at for i in list_survived_atoms if(i[9] == j[1])]
    cmaps_prot_fully_at = [j + [i[7]] for j in cmaps_prot_fully_at for i in list_survived_atoms if(i[9] == j[2])]
    cmaps_prot_fully_at = [j + [i[7]] for j in cmaps_prot_fully_at for i in list_survived_atoms if(i[9] == j[3])]
    cmaps_prot_fully_at = [j + [i[7]] for j in cmaps_prot_fully_at for i in list_survived_atoms if(i[9] == j[4])]


    #E.6- Adding the new at_type (A{$i}) in case of CG beads. The at_type does not change for AT atoms. 

    for i in list_survived_atoms:
        for j in cmaps_prot_fully_at:
            if(i[9] == j[0]):
                j[6] = i[16]

    for i in list_survived_atoms:
        for j in cmaps_prot_fully_at:
            if(i[9] == j[1]):
                j[7] = i[16]

    for i in list_survived_atoms:
        for j in cmaps_prot_fully_at:
            if(i[9] == j[2]):
                j[8] = i[16]

    for i in list_survived_atoms:
        for j in cmaps_prot_fully_at:
            if(i[9] == j[3]):
                j[9] = i[16]

    for i in list_survived_atoms:
        for j in cmaps_prot_fully_at:
            if(i[9] == j[4]):
                j[10] = i[16]

    #E.7- Creating a new list called "cmaps_mult_res_prot" (in analogy with bonds, pairs, angles, and dihedrals) which is simply cmaps_prot_fully_at, 
    #     and sorting the list. 

    cmaps_mult_res_prot = cmaps_prot_fully_at
    cmaps_mult_res_prot = sorted(cmaps_mult_res_prot, key=itemgetter(0,1,2,3,4))  # sorting 1st, then 2nd, 3rd, 4th, and finally 5th column.

    print("\nlist of Correction Torsional Map (cmap) for charmm forcefield updated... 85% completed\n")

else: # i.e. the "original_cmaps_protein_fully_at" list is empty, i.e. the charmm.ff is not employed
    cmaps_mult_res_prot = []

  #########################################
####  I. CREATING posre_mul_res.itp FILE ##############
  ######################################### 

""" 
    Creating 'posre_mul_res.itp' file, whose purpose is to apply a position restraining force on the heavy atoms of the protein 
    (anything that is not a hydrogen). Movement is permitted, but only after overcoming a substantial energy penalty.
    The utility of position restraints is that they allow us to equilibrate our solvent around our protein, 
    without the added variable of structural changes in the protein.
"""

f_posre = open("posre_mul_res.itp", "w")
write_posre_new_file(list_survived_atoms)
f_posre.close()


  #########################################
####  II. CREATING ffbonded_new.itp FILE ##############
  #########################################

"""
   Creating the file of bonded parameters (bond- angle- dihedraltypes) in which is involved 
   at least one CG bead. 
"""

fb_new = open("ffbonded_new.itp", "w")
write_ffbonded_new_file(fb_new, bonds_mult_res_prot, angles_mult_res_prot, dihedrals_mult_res_prot)
fb_new.close()

  ##############################################
####  III. CREATING ffnonbonded_new.itp FILE ##############
  ##############################################

"""
   Creating the file of nobonded parameters (epsilon and sigma) for each CG bead.
"""

fnb_new = open("ffnonbonded_new.itp", "w")
write_ffnonbondednew_file(list_survived_atoms)
fnb_new.close()

  ####################################
####  IV. CREATING mult.itp FILE    ##############
  ####################################

"""
   Creating 'mult.itp' file. that is the topology file for our system that includes bonds, angles, and dihedrals, and cmaps (only for Charmm.ff)  
"""

f_mult_itp = open("mult.itp", "w")
write_mult_itp_file(f_mult_itp, list_survived_atoms, bonds_mult_res_prot, pairs_mult_res_prot, angles_mult_res_prot, dihedrals_mult_res_prot,\
                    k_strong, cmaps_mult_res_prot)
f_mult_itp.close() 

  #######################################
####  V. CREATING cmap_new.itp FILE    ##############
  #######################################

if(Flag_cmaps == True):   # i.e. the "original_cmaps_protein_fully_at" list is not-empty, i.e. the charmm.ff is employed  
    f_cmap_new = open("cmap_new.itp", "w") 
    write_new_cmap_file(f_cmap_new, cmaps_mult_res_prot)
    f_cmap_new.close()

  ########################################
####  V. CREATING topol_new.top FILE    ##############
  ########################################

"""
   Creating the file 'topol_new.top' that includes all the files written in the previous parts 
   i.e. mult.itp, ffnonbonded_new.itp, and ffbonded_new.itp. Moreover, we add the atomistic forcefield, the 
   topology for ions and tip3p water, and finally a summary of used molecules (usually, Proteins, SOL and ions)
   Since so far system is only solvated, we include only the number of protein and SOL molecules. The number of 
   ions is includes automatically by the system after neutralizing the system with ions (later) 
"""

ft_new = open("topol_new.top", "w")
write_new_topology_file(ft_new, path_ff_new, count_SOL, Flag_cmaps)
ft_new.close() 

print("\ntopology file with all inclusions created and updated... 89% completed\n")


  #################################################
####  VI. CREATING VdW radius and charge FILE    ##############
  #################################################
"""
   Creating two files:
       1- One file containing the VdW radius (sigma/2) and its charge for each atom present in the Mult-res configuration; 
 
       2- A second file containing the Tkl script for reading the file (1) and for representing each atom in VDW with radius 
          given by the value of sigma/2, and an atom coloration given by the Beta Coloring Method: positive charges are 
          coloured with shades of blue, neutral charges in white, and negative charges with shades of red. 
"""

# VI.a- 1st File 

f_rc = open("radius_charge-allatom.txt", "w")
write_radius_charges(f_rc, list_survived_atoms) 
f_rc.close() 


# VI.b- 2nd File 

f_tcl_rc = open("radius_charge-allatom.tcl", "w") 
write_tcl_radius_charges(f_tcl_rc) 
f_tcl_rc.close()

  ##########################################################
####  VII. CREATING CA-list FILE & rmsf tkl script FILE   ##############
  ##########################################################
"""
   Creating two files: 
    1- One file containing the list of CA atomistic number for CANVAS model, useful per RMSD and RMSF of the initial Mul-Res configuration. 

    2- A second file containing the Tkl script for performing the RSMF of only the CA (atomistic, medium-grained and coarse-grained) 
       of Mul-Res trajectory. 
"""

# VII.a- 1st File 

f_CA = open("CA_list.txt", "w")
write_CA(f_CA, CA_list)
f_CA.close() 


# VII.b- 2nd File 

f_tcl_rmsf = open("rmsf.tcl", "w") 
write_rmsf(f_tcl_rmsf, CA_list) 
f_tcl_rmsf.close() 

  ###############################################
####  VI. ADDING IONS TO SOLVATED .gro FILE  ##############
  ###############################################

"""
   Creating, 'ions.mdp' file that needs to add ions to the system in two steps in Gromacs. 
"""

if(solvate_Flag == True):
    # VI.a- Creation of ions.mdp file 

    f_ions = open("ions.mdp", "w+")
    write_ions_mdp(f_ions) 
    f_ions.close()


    # VI.b- Adding ions to the system. Topology file is updated automatically. 

    os.system("{} grompp -f ions.mdp -c newsolvated.gro -p topol_new.top -o ions.tpr > /dev/null 2>&1".format(gmx_command))

    if os.path.isfile("ions.tpr"):
        os.system("echo SOL | {} genion -s ions.tpr -o solvated_ions.gro -p topol_new.top -pname NA -nname CL -neutral > /dev/null 2>&1".format(gmx_command))

        print("\nSolvated file with ions completed... 90% completed\n")

    if not os.path.isfile("ions.tpr"):
        print("WARNING: Something went wrong. It is not possible put ions for neutralizing the system.")
        print("         Maybe the system is too big, or the number of atomtypes is very high, or it is required more RAM.")
        print("         Please, try the launch the following command on terminal for understanding wht the problem occours:\n")
        print("         {} grompp -f ions.mdp -c newsolvated.gro -p topol_new.top -o ions.tpr \n".format(gmx_command))
        print("         As the file with ions has not been created, the solvated file (with no ions) will be employed for writing the lammps file...\n")

        print("\nSolvated file with no ions completed... 90% completed\n")


    print("\nSolvated file with ions completed... 90% completed\n")


  #####################################################################
####  VII. CREATING minimization, equilibration, and run .mdp FILES  ##############
  #####################################################################

"""
   Creating four files whose purpose is to minimize the system (em.mdp), then equilibrate it (nvt.mdp and npt.mdp) 
   and finally for doing a long run (run.mdp). In all cases, the value of coulomb and VdW cutoffs are identical and equals to 2.5 * max_sigma.  
"""

if(FlagLammps == False): 
    # VII.a- Creating em.mdp file  

    f_em = open("em.mdp","w")
    write_em_mdp(f_em, max_sigma) 
    f_em.close()
    
    
    # VII.b- Creating nvt.mdp file 
    
    f_nvt = open("nvt.mdp", "w")
    write_nvt_mdp(f_nvt, max_sigma)
    f_nvt.close()
    
    
    # VII.c- Creating npt.mdp file
    
    f_npt = open("npt.mdp", "w")
    write_npt_mdp(f_npt, max_sigma)
    f_npt.close()
     
    
    # VII.d- Creating run.mdp file
    
    f_run = open("run.mdp", "w")
    write_run_mdp(f_run, max_sigma)
    f_run.close()

print("\nmdp files created... 98% completed\n")

  ############################################################
####  VIII. CREATING running directory and copy files in it ##############
  ############################################################

"""
     o Copying the file residuetypes.dat from gromacs path, and updating it. 
     o Creating running directory 
     o Copying all files required to simulate the system in the running directory 
"""

# VIII.a- Copy residuetypes.dat from Gromacs path in our simulating folder. Moreover, we add in the last row, 
#         the new residue name for all CG beads (i.e. **MUL**) with the appropriate specification **Protein**.

gromacs_folder = ''.join(path_ff.split(forcefield))  
residuetypes   = gromacs_folder + "residuetypes.dat"

os.system("cp {} .".format(residuetypes))
os.system("echo 'MUL "" "" "" "" Protein' >> residuetypes.dat")  # only way to add also 4 extra spaces between MUL and Protein. 


# VIII.b- Creating running directory (and remove it if already present) 

if(os.path.isdir("run_simulation") == True):
    os.system("rm -rf run_simulation")

os.system("mkdir run_simulation")


# VIII.c- Copy all files required to simulate the system in the running directory 
if(FlagLammps == False): 
    os.system("mv em.mdp nvt.mdp npt.mdp run.mdp run_simulation")
    os.system("mv residuetypes.dat ffnonbonded_new.itp ffbonded_new.itp run_simulation")
    os.system("mv mult.itp posre_mul_res.itp topol_new.top run_simulation")
    os.system("mv *gro other-files/") 
 
    if(Flag_cmaps == True):   				 # i.e. the "original_cmaps_protein_fully_at" list is not-empty, i.e. the charmm.ff is employed 
        os.system("mv cmap_new.itp run_simulation")      # therefore "cmap_new.itp" exists: 

    
os.system("rm \#*")
    
os.system("mkdir other-files")
os.system("mv *txt *mdp *tpr other-files/")
    
os.system("mkdir analysis")
os.system("mv other-files/radius_charge-allatom.txt analysis/")  # Moving this file because this file has to be used in the analysis section. 
os.system("mv other-files/CA_list.txt analysis/")                # Moving this file because this file has to be used in the analysis section.
os.system("mv *tcl analysis/")


#### CREATING FILES FOR LAMMPS ######

if(FlagLammps == True):
    print("\nCreating input files for simulating system in LAMMPS... 99% completed\n")  
    # L1- Finding correct path of water and ions topologies inside the main forcefield folder.   

    tip3pfile = check_tip3p(path_ff_new)
    ionsfile  = check_ions(path_ff_new) 

    
    # L2- Opening and Reading topology for ions, namely the "ions.itp" file

    f_ions_itp = open(ionsfile, "r")
    
    top_ions               = readtop_ions(f_ions_itp)
    atoms_topol_ions       = top_ions

    # L3- Opening and Reading topology for water, that is "tip3p.itp" file 

    f_water_itp = open(tip3pfile, "r")
    
    toop_wat = readtop_water(f_water_itp)
    
    atoms_topol_water       = toop_wat[0]
    original_bonds_water    = toop_wat[1]
    original_angles_water   = toop_wat[2]

    # L4- Reading solvated file if present (incl. ions for neutralizing charge)
    if(solvate_Flag == True):  

        if os.path.isfile("solvated_ions.gro"):
            f_water_ions = open("solvated_ions.gro", "r")
        else:
            f_water_ions = open("newsolvated.gro", "r")
        
        gro_ions = readgro_all(f_water_ions)
        
        atoms_ions       = gro_ions[0]    #atoms = [res_number, res_name, at_name, at_number, pos_x, pos_y, pos_z]
        n_atoms_ions     = gro_ions[1] 

    # L5- Updating "list_survived_atoms" with all the properties for water and ions, i.e.: 
    #     [res_number, res_name, at_name, at_number, pos_x, pos_y, pos_z, cg/at, at_type, new_at_num, mass, \
    #     charge, sigma, epsilon, one/four/atom, atomic_number, new_at_new, new_res_num ]   

    if(solvate_Flag == True):            
        ## L5.a- changing X,Y,Z coordinates for the protein using the new solvated with ions coordinate file 
        for i in atoms_ions:
            for j in list_survived_atoms:
                if(j[9] == i[3] and i[1]!="SOL" and i[1]!="NA" and i[1]!="CL"):    # if new_attype == attype of solvated file-
                    j[4] = i[4]                   
                    j[5] = i[5]
                    j[6] = i[6]

    if(solvate_Flag == False): 
        for i in atoms_box:
            for j in list_survived_atoms:
                if(j[9] == i[3]):
                    j[4] = i[4]
                    j[5] = i[5]
                    j[6] = i[6] 


        
    if(solvate_Flag == True):      
        ## L5.b- Appending [res_number, res_name, at_name, at_number, pos_x, pos_y, pos_z, "at"] for ions and water 
        for i in atoms_ions:
            if(i[1]=="SOL" or i[1]=="NA" or i[1]=="CL"):
                list_survived_atoms.append([i[0], i[1], i[2], i[3], i[4], i[5], i[6], "at"])
        
        
        ## L5.c- Appending "at_type" (i[8]) for water and ions  
        for j in atoms_topol_water:
            for i in list_survived_atoms:
                if(j[3]==i[1]):                   			# If "SOL" is found in both lists 
                    if(j[4] ==i[2]):
                        i.append(j[1])           		        # j[1] is the column relative at the at_type of each atom 
        
        for j in atoms_topol_ions:
            for i in list_survived_atoms:
                if(j[3] == i[1]):
                    i.append(j[1])
        
        
        ## L5.d- Appending new_at_num (i[9]) for water and ions (it is the same of the at_number found in solvated_ions.gro file)       
        count = 1
        for i in list_survived_atoms:
            if(len(i) == 9):
                i.append(count)
            count = count + 1
        
        
        ## L5.e- Appending mass of water and ions (i[10]) (dict_mass = ffnonbon[2]) 
        for k,v in dict_mass.items():
            for i in list_survived_atoms:
                if(len(i) == 10 and k==i[8]):
                    i.append(v)
        
        
        ## L5.f- Appending charge for water and Ions (i[11]) 
        for j in atoms_topol_water:
            for i in list_survived_atoms:
                if(len(i) == 11 and j[4]==i[2]):       
                    i.append(j[6])
        
        for j in atoms_topol_ions:
            for i in list_survived_atoms:
                if(len(i) == 11 and j[4]==i[2]):       
                    i.append(j[6])
        
        
        ## L5.g- Appending sigma for Water and Ions (i[12]) 
        for k,v in dict_sigma.items():
            for i in list_survived_atoms:
                if(len(i) == 12 and k == i[8]):
                    i.append(v)                
        
        
        ## L5.h- Appending epsilon for Water and Ions (i[13]) 
        for k,v in dict_epsilon.items():
            for i in list_survived_atoms:
                if(len(i) == 13 and k==i[8]):                                    # If the at_type is the same in dict_epsilon and list_other_atoms
                    i.append(v)                                                  # append the original value of epsilon of the specific at_type   
        
                                                                      
        ## L5.i- Appending "atom" for Water and Ions (i[14])
        for i in list_survived_atoms:
            if(len(i) == 14):
                i.append("atom")
        
        
        ## L5.j- Appending the atomic number (H = 1, C =6) for Water and Ions (i[15])
        for k,v in dict_atnum.items():
            for i in list_survived_atoms:
                if(len(i) == 15 and i[8]==k):   			         # If the at_type (i[8]) is the same 
                    i.append(v)                     				 # Append "v" that is the atomic number 
        
        
        ## L5.k- Appending the new at_number for Water and Ions (i[16]) (it is the same of at_type i[8] beacuse they are atomistic) 
        for i in list_survived_atoms:
           if(len(i) == 16):
               i.append(i[8])
        
        
        ## L5.l- Appending the "new res_number" for Water and Ions (that is the same of old res_number i[0] by construction) 
        for i in list_survived_atoms:
            if(len(i) == 17):
                i.append(i[0])


    # L6- Calculating the "Number of Atomtypes" namely only the total number of atomtypes of survived atoms 
    #     plus the atomtypes of the Water (HW and OW usally) and Ions (Na, Cl, etc.). A list of atomtypes is then created. 
      
    ## L6.a- Creating list of atomtypes
    list_at_types = []
      
    for i in list_survived_atoms:
        list_at_types.append(i[16])
    
    list_at_types = list(dict.fromkeys(list_at_types))   # remove duplicates without changing the order of the elements found before. 
         						 # SET operations, on the other hand, is not ordered and the elemenets change the order 
    
    ## L6.b- Creating a dictionary such that each atomtype has a number that goes between 1 and N.  
    dict_at_types = {j+1 : list_at_types[j] for j in range(0,len(list_at_types))}
    
    
    ## L6.c- Appending, the number of atomtype (NUM_OF_ATTYPE) for each suvived atom 
    for i in list_survived_atoms:
        for k,v in dict_at_types.items():
            if(v == i[16]):
                i.append(k)
    
    ## L6.d- Calculating the total number of atomtypes, needed the be declared in the lammps .lmp file 
    number_atom_types = len(dict_at_types)


        
    # L7- Number of Bondtypes. Updating the part relative to the bonds for water (OW-HW1 and OW-HW2), 
    #     appending the result in "bonds_mult_res_prot"

    ## L7.a- Usually the at_type of the water molecules are HW and OW, but in charmm36m, the latter are 
    ##       denoted as HT and OT. In order to guarantee no mistakes, in this step type_OW and type_OW strings are computed. 
    
    if(solvate_Flag == True):     
        for i in list_survived_atoms:
            if(i[1] == "SOL"): 
                if(i[2].startswith("O")): 
                    type_OW = i[8]
        
                if(i[2].startswith("H")):
                    type_HW = i[8]     
                    break 

        ## L7.b- Updating "bonds_mult_res_prot" with indexes that involve OW-HW1 and OW-HW2. 
        for i in list_survived_atoms:
            if(i[8] == type_OW):
        
                bonds_SOL1 = [i[9], i[9] + 1]     			# bond between OW e HW1 (indexes j e j+1) 
                bonds_SOL2 = [i[9], i[9] + 2]    	 		# bond between OW e HW2 (indexes j e j+2) 
        
                bonds_mult_res_prot.append(bonds_SOL1)
                bonds_mult_res_prot.append(bonds_SOL2)
        
        ## L7.c- Appending the at_type of each couple of Water atoms (i[2] of list_survived_atoms):
        ##       By construction, the 1st element is always OW, the second one is HW
        for k in bonds_mult_res_prot:
            if(len(k) == 2):
                k.append(type_OW)      			        	# at_type in bonds_mult_res_prot
                k.append(type_HW)      				        # at_type in bonds_mult_res_prot
        
        
        ## L7.d- Appending "at" and "at" because, we know that the water atoms are always atomistic 
        for i in bonds_mult_res_prot:
            if(len(i) == 4):
                i.append("at")
                i.append("at")
        
        
        ## L7.e- Appending kb and b0 for the OW-HW bond (it can be found in "bondtypes")  
        for i in bonds_mult_res_prot:
            for j in bondtypes:
                if(len(i) == 6):
                    if((j[0]==i[2] and j[1]==i[3]) or (j[0]==i[3] and j[1]==i[2])):
        
                        i.append(float(j[3]*836.8))      		# kb value (in Gromacs units)
                        i.append(float(j[2]/10))         		# b0 value (in Gromacs units)
        
        
    ## L7.f- Appending the corresponding value of bond_type for each bond, and counting them.
    ##       For simplicity, each bond that involves the protein will have a different bondtype_number 
    ##       On the other hand, the OW-HW bondtype_number for Water, must be considered only one. 

    count = 1
    for i in bonds_mult_res_prot:
        if(solvate_Flag == True): 
            if(i[2]==type_OW or i[2]==type_HW or i[3]==type_OW or i[3]==type_HW):   # For OW-HW bond, the bondtype_number is not increased 
                i.append(count)    
            else:
                i.append(count)
                count = count +1

        if(solvate_Flag == False): 
            i.append(count) 
            count = count + 1 

    ## L7.g- Calculating the total number of bondtypes, needed the be declared in the lammps .lmp file
    if(solvate_Flag ==True): 
        number_bond_types = count
   
    else: 
        number_bond_types = count - 1

 
    # L8- Number of Angletypes: Updating the part relative to the angles for water (HW1-OW-HW1), 
    #     appending the result in "angles_mult_res_prot"

    if(solvate_Flag == True): 
        ## L8.a- Updating "angles_mult_res_prot" with indexes that involve HW1-OW-HW2  
        for i in list_survived_atoms:
            if(i[8] == type_OW):
                angles_SOL = [i[9] + 1, i[9], i[9] + 2]     			# bond between OW e HW1 (indici j e j+1) 
                angles_mult_res_prot.append(angles_SOL)
        
        
        ## L8.b- Appending the at_type of each triplet of Water atoms (i[2] of list_survived_atoms):
        ##       By construction, the 1st element is always HW, the second one is OW, the third one is HW
        for k in angles_mult_res_prot:
            if(len(k) == 3):
                k.append(type_HW)      						# at_type in bonds_mult_res_prot
                k.append(type_OW)      						# at_type in bonds_mult_res_prot
                k.append(type_HW)
        
        ## L8.c- Appending cth (Angle Energy) and th0 (Theta Angle) for the HW-OW-HW (water) angle (it can be found in "angletypes")
        for i in angles_mult_res_prot:
            for j in angletypes:
                if(len(i) == 6):
                    if((j[0]==i[3] and j[1]==i[4] and j[2]==i[5]) or (j[0]==i[5] and j[1]==i[4] and j[2]==i[3])):
                        i.append(float(j[5]*8.368))      			# cth value (in Gromacs units)
                        i.append(float(j[4]))            			# th0 value (in Gromacs units)
        
        
        ## L8.d- Appending "at","at","at" because, we know that the water atoms are always atomistic 
        for i in angles_mult_res_prot:
            if(len(i) == 8):
                i.append("at")
                i.append("at")
                i.append("at")
        
    ## L8.e- Appending the corresponding value of angle_type for each angle, and counting them.
    ##       For simplicity, each angle that involves the protein will have a different angletype_number 
    ##       On the other hand, the HW-OW-HW angletype_number for Water, must be considered only one. 
    count = 1
    for i in angles_mult_res_prot:
        if(solvate_Flag == True): 
            if(i[2]==type_OW or i[2]==type_HW or i[3]==type_OW or i[3]==type_HW or i[4]==type_OW or i[4]==type_HW):   # In case of HW-OW-HW angle, "count" 
                i.append(count)							                                      # is not increased 		 
            else:
                i.append(count)
                count = count +1

        if(solvate_Flag == False): 
            i.append(count) 
            count = count + 1 

    
    ## L8.f- Calculating the total number of angletypes, needed the be declared in the lammps .lmp file
    if(solvate_Flag == True): 
        number_angle_types = count
    else:
        number_angle_types = count - 1     

   
 
    # L9- Number of Dihedraltypes & Impropers: No Dihedrals are present in Water, therefore only protein can have Dihedrals.
    #     No updates needed. However, dihedrals with functional form equals to 9 and 4 will be treated as DIHEDRALS in Lammps, 
    #     while dihedrals with functional form equals to 2 will be treated as IMPROPERS. 
   
    ## L9.a- Splitting count for dih2 and dih4/9, and increase count for dih2 and dih4/9 
    count_dih2 = 1
    count_dih4_dih9 = 1 
    
    for i in dihedrals_mult_res_prot:
        if(i[4] == 2): 
            i.append(count_dih2)
            count_dih2 = count_dih2 + 1
        if(i[4] == 4 or i[4]==9):
            i.append(count_dih4_dih9)
            count_dih4_dih9 = count_dih4_dih9 + 1 
    
    ## L9.b- Calculating the total number of dihedraltypes and impropertypes, needed the be declared in the lammps .lmp file
    number_dihedral_types = count_dih4_dih9 -1         # "-1" because after reading the last dih, the final count is increased by 1 with respect the value  
    number_improper_types = count_dih2 -1 



    # L10- Number of Atoms: It is given by the sum of survived atoms, water atoms and ions. It is the lenght of "list_survived_atoms"
    number_of_atoms   = len(list_survived_atoms)
    
    # L11- Number of Bonds: It is given by the sum of all bonds (for protein and water). It is the lenght of "bonds_mult_res_prot"
    number_of_bonds   = len(bonds_mult_res_prot)
    
    # L12- Number of Angles: It is given by the sum of all angles (protein and water). It is the lenght of "angles_mult_res_prot"
    number_of_angles  = len(angles_mult_res_prot)
    
    # L13- Number of Dihedrals: It is given by the sum of dihedrals with funct=4 and funct=9
    #      Number of Impropers: It is given by the sum of dihedrals with funct=2 
    number_of_dihedrals = 0
    number_of_impropers = 0
    
    for i in dihedrals_mult_res_prot: 
        if(i[4] == 9 or i[4] == 4):
            number_of_dihedrals = number_of_dihedrals + 1
        else:
            number_of_impropers = number_of_impropers + 1

    
    # L14- Box Sizes: Mimimum and Maximum corrdinates are read from "solvated_ions.gro". In general, 
    #                 min_x, min_y, min_z are, indeed, the minimum value of x,y,z respectively.
    # 
    #                 On the other hand: max_x = min_x + Lx     where Lx,Ly,Lz is the size lenght.  
    #                                  : max_y = min_y + Ly 
    #                                  : max_z = min_z + Lz

    if(solvate_Flag == True):  
        min_x = gro_ions[6]
        min_y = gro_ions[7]
        min_z = gro_ions[8]
        
        max_x = gro_ions[9]
        max_y = gro_ions[10]
        max_z = gro_ions[11]

    if(solvate_Flag == False):  
        min_x = gro_box[6]
        min_y = gro_box[7]
        min_z = gro_box[8]
   
        max_x = gro_box[9]
        max_y = gro_box[10]
        max_z = gro_box[11]

    
    # L15- Opening and Writing .lmp file in Lammps 

    f_lmp = open("Multiple-Res-Protein-in-Water.lmp", "w")

    if(solvate_Flag == False):
        type_OW = "OW"
        type_HW = "HW"   # We need then when writing lmp_file and input_file, but they will be not used since there are no water molecules    

        write_lmp_file(f_lmp, number_of_atoms, number_of_bonds, number_of_angles, number_of_dihedrals, number_of_impropers, number_atom_types, \
                       number_bond_types, number_angle_types, number_dihedral_types, number_improper_types, min_x, max_x, min_y, max_y, \
                       min_z, max_z, list_survived_atoms, bonds_mult_res_prot, angles_mult_res_prot, dihedrals_mult_res_prot, atoms_box, \
                       type_OW, type_HW)

    if(solvate_Flag == True):
        write_lmp_file(f_lmp, number_of_atoms,number_of_bonds, number_of_angles, number_of_dihedrals, number_of_impropers, number_atom_types, \
                       number_bond_types, number_angle_types, number_dihedral_types, number_improper_types, min_x, max_x, min_y, max_y, \
                       min_z, max_z, list_survived_atoms, bonds_mult_res_prot, angles_mult_res_prot, dihedrals_mult_res_prot, atoms_ions, \
                       type_OW, type_HW)
    
    f_lmp.close() 


    # L16- Opening and Writing .input file for Lammps
 
    f_input = open("Multiple-Res-Protein-in-Water.input", "w")
    
    write_input_lammps_file(f_input, dihedrals_mult_res_prot, Flag_cmaps, number_of_bonds, number_of_angles, number_of_dihedrals, \
                            number_of_impropers, number_atom_types, list_survived_atoms, pairs_mult_res_prot, dict_at_types, \
                            bonds_mult_res_prot, angles_mult_res_prot, type_OW, type_HW, solvate_Flag)
    
    f_input.close() 


if(FlagLammps == False): 
    if(solvate_Flag == True): 
        os.system("mv solvated_ions.gro run_simulation/")
    else: 
        os.system("mv newbox.gro No_Solvated_ions.gro")
        os.system("mv No_Solvated_ions.gro run_simulation/")

    os.system("mv *gro other-files/") 

if(FlagLammps == True): 
    os.system("mv *lmp *input run_simulation") 
    os.system("mv *itp *gro *dat *top other-files/") 


 
print("\nNo errors...100% completed\n")

REAL_end = datetime.now()

print("The time for executing the code is: ", (REAL_end - REAL_start).total_seconds())

f_time = open("time.txt", "w") 
f_time.write("time = {} sec".format((REAL_end - REAL_start).total_seconds()))

os.system("mv time.txt other-files/")
