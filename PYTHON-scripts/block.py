import sys 
import re 
import argparse
import itertools 
import os 
import itertools   		      # It needs for removing duplicates in list of lists.

from operator import itemgetter
from math import sqrt 
from datetime import datetime

# --------------------------------------------

PYTHONPATH = os.path.abspath(os.getcwd())

REAL_start = datetime.now()

spl_word = "canvas" # We find CANVAS, cut the entire path until CANVAS and add /lib in order to find our libraries. 

python_modules_path = PYTHONPATH.split(spl_word)[0] + spl_word + "/lib"

sys.path.append(python_modules_path)


from read_gromacs import * 
from inp_out import * 

# 2. Input Arguments -------------------------------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, add_help=False) 

group_in=parser.add_argument_group("Required Arguments") 

group_in.add_argument('option', help = argparse.SUPPRESS)					                          # Mandatory
group_in.add_argument('-g', '--gro', dest='grofile', action='store', metavar = 'FILE', help = argparse.SUPPRESS)          # Mandatory
group_in.add_argument('-l', '--list', dest='listfile', metavar = 'FILE',help = argparse.SUPPRESS)                         # Mandatory 

group_in.add_argument('-h', '--help', action='help', help = argparse.SUPPRESS)                                            # Optional
group_in.add_argument('-D', '--diameter', dest='DiameterValue', metavar = 'VALUE', type=float, help = argparse.SUPPRESS)  # Optional 

# ------------------------------------------------------------------------------------------------------------------------------------------

print("\n####################################################\n")
print("'{}' running...".format(os.path.basename(sys.argv[0]))) 
print("---------------------------------\n")

# 3. Printing help message if the script does not have any arguments  
if len(sys.argv)==1:
    print_usage_main_block()
    sys.exit(1)


# 4. Printing help message if the script does not present valid arguments 
found = False
tasks = ["choice1", "choice2", "choice3"]

for tk in tasks: 
    if(sys.argv[1] == tk):
            found = True;

if(found == False):        
    
    if(sys.argv[1].strip() == "--usage") or (sys.argv[1].strip() == "-u"):
        print_usage_main_block()
        quit() 

    if(sys.argv[1].strip() == "--help" or sys.argv[1].strip() == "-h"):
        print_help_main_block()
        quit() 

    print("Error. {} not in the list of accepted tasks\n".format(sys.argv[1]))
    print_usage_main_block()
    quit() 

check_argv_errors_block()


# 5. Parsing arguments, printing error and help message if the script presents more than one task, and checking if mandatory files are present  
try:
    args  = parser.parse_args()
except SystemExit:
    print("\nOnly one task is allowed. Check that each flag (e.g. -g) is followed by its specific argument")
    print("Look below for more information.\n")
    print_help_block() 
    sys.exit()  

option   = args.option           # Mandatory 
grofile  = args.grofile          # Mandatory
listfile = args.listfile         # Mandatory

diameter = args.DiameterValue    # Optional 

mandatory_files_present_block(grofile, listfile) 

## ----------------------------------------------------------------------------------------------------------------------------------



# 6. Choice1 selected... 
if(option == "choice1"):

    # 6.1 - Checking if grofile is present  
    if not os.path.isfile(grofile):
        print("Error while opening the file. '{}' does not exist\n".format(grofile))
        quit()
    
    gro_f = open(grofile, "r") 
    check_empty_file(grofile) 
   
    gro      = readgro_all(gro_f)

    atoms    = gro[0] 
    n_atoms  = gro[1]

    print("\no '{}' correctly read... 10% completed\n".format(os.path.basename(grofile))) 

    # 6.2 - Opening and reading the file showing the list of residue(s) (not atoms), that correspond at the center of atomistic region,
    #       from which the program will trace the atomistic sphere. 
    if not os.path.isfile(listfile):
       print("Error while opening the file. '{}' does not exist\n".format(listfile))
       quit()

    f = open(listfile, "r")

    list_AT_res      = []
    list_radii       = []
   
    count_line = 1  
    for line in f:
        check_empty_rows_block(line, listfile)
        line = line.split()

        if(len(line) != 2):       #N_lines = 2 
            print("*choice1* selected: the input file '{}' cannot be read correctly.".format(listfile))
            print("\nError. Each line must contain 2 columns: the {}th row has {} columns".format(count_line, len(line)))
            print("\nBe sure that '{}' is organized in two different columns\n".format(listfile))
            print(" • 1st column       INT                List of residues numbers from which the atomistic sphere will be traced.")
            print(" • 2nd column       INT/FLOAT          Atomistic radius of each residue in nm.") 
            print("\nLook below for further information\n")
            print_help_block() 
            quit() 

        else:
            if not line[0].isdigit():  #1st column = INT
                print("*choice1* selected: the input file '{}' cannot be read correctly.".format(listfile))
                print("\nError. The 1st column must be an integer number. '{}' at {}th row is not". format(line[0], count_line))
                print("\nBe sure that '{}' is organized in two different columns\n".format(listfile))
                print(" • 1st column       INT                List of residues numbers from which the atomistic sphere will be traced.")
                print(" • 2nd column       INT/FLOAT          Atomistic radius of each residue in nm.") 
                print("\nLook below for further information\n")
                print_help_block()
                quit()

            if not isIntFloat(line[1]): #2nd column = INT/FLOAT
                print("*choice1* selected: the input file '{}' cannot be read correctly.".format(listfile))
                print("\nError. The 2nd column must be a Float or Integer number. '{}' at {}th row is not". format(line[1], count_line)) 
                print("\nBe sure that '{}' is organized in two different columns\n".format(listfile))
                print(" • 1st column       INT                List of residues numbers from which the atomistic sphere will be traced.")
                print(" • 2nd column       INT/FLOAT          Atomistic radius of each residue in nm.") 
                print("\nLook below for further information\n")
                print_help_block() 
                quit()  

        res    = int(line[0])
        radius = float(line[1]) 
  
        list_AT_res.append(res)
        list_radii.append(radius) 

        count_line = count_line + 1 

    print("\no '{}' correctly read... 20% completed\n".format(os.path.basename(listfile)))

    # 6.3 - Creation of a list corresponding at the C-alpha atoms, kept from the list of the atomistic residues (list_AT_residues)  
    list_AT_atoms = [j for i in list_AT_res for j in atoms if (j[0] == i)]
 
    list_other_atoms = [x for x in atoms if not x in list_AT_atoms]


    # 6.4 - Computing the distance between the atomistic atoms (list_AT_atoms) and the other atoms (list_other_atoms): 
    #
    #       If distance < atomistic_radius             => list_fullyAT_residues = list of residues that require an atomistic description 
    #
    #       If  at_radius < distance < bb Diameter (D) => list_bb_residues =  list of residues that require an hybrid description 
    #
    #       If distance > bb Diameter (D)              => list_CG_residue = list of residues that require a CG description    
    #

    if(diameter is not None):  
        D = diameter 
        if(D < 0): 
            print("A negative value of diameter do not make any sense. Change it.\n\n")
            print_help_block()
            quit() 

    else:   
        D = 1.0   # 1 nm: default value of the diameter of the transition area modelled with only the backbone atoms.  

    print("\no Diameter D of hybrid region: {} nm... 30% completed\n".format(D))

    D = float(D)

    
    print("\no Computing distances... 35% completed\n") 

    list_AT_residues = []
    list_bb_residues = []
    list_CG_residues = []

    count = 0 

    for i in list_AT_atoms:
        for j in list_other_atoms:

            x   = j[4]-i[4]
            y   = j[5]-i[5]
            z   = j[6]-i[6]

            dist = sqrt(x**2 + y**2 + z**2)

            if(dist < float(list_radii[count])):      # at the beginning: count=0  and  we take list_radii[0], and so on until the endfile
                list_AT_residues.append(j[0])         # append res number 

            if(dist >= float(list_radii[count]) and dist < float(list_radii[count]) + D):   # 1.0 nm is the default radius of backbone part.
                list_bb_residues.append(j[0])

            if(dist >= float(list_radii[count]) + D):   
                list_CG_residues.append(j[0]) 

        if(i[2] == "O" or i[2] == "OC2"): 			      # It means that the next residue has been reached (OC2 in case of terminal res) 
            count = count + 1			                      # then, count can be increased by one 



    list_AT_residues = list_AT_residues + list_AT_res  # Adding in list_AT_residue the atomistic residues writtein in file.   


    # 6.5 - Removing repetitions because some atoms have the same residue number and sorting out the three list in ascending order.

    list_AT_residues = list(set(list_AT_residues))
    list_bb_residues = list(set(list_bb_residues))
    list_CG_residues = list(set(list_CG_residues))

    list_AT_residues.sort() 
    list_bb_residues.sort() 
    list_CG_residues.sort() 

    # 6.6 - Summarizing, there are the following list of residues:
    #
    #        list_AT_residues            => It contains the list of residues that has to be modelled atomistically,
    #                                       therefore all the atoms of a specific amino acid will be retained (they survive) 
    #
    #        list_bb_residues            => It contains the list of residues that has to be modelled as CG but retaining 
    #                                       only CA, C, O and N of a specific amino acid. 
    #
    #        list_CG_residues            => It contains the list of residue that has to be modelled as CG, but at difference 
    #                                       with the previous case, we retain only the CA of a specific amino acid.  
    #
    #        Assuming that two or more resolutions (AT, bb or CG) overlap each other, the higher resolution is kept: 
    #
    #        AT > BB > CG 
    #
    #        In order to do it easily, we transform the previous three list as set: in this way, than by applying just the difference 
    #        between sets, it will be possible to know only the residues that will be present ONLY in each set. The operation mentioned 
    #        between sets is, indeed, the difference. In other words: 
    #
    #        set_CG_residues = set_CG_residues - set_bb_residues - set_AT_residues   => this is the set of residues that belongs to the 
    # 											CG only, and not at AT and backbone 
    #
    #        set_bb_residues = set_bb_residues - set_AT_residues                     => this is the set of residues that belongs to the 
    # 										        backbone only, and not at AT. 
    #
    
    set_AT_residues  = set(list_AT_residues)
    set_bb_residues  = set(list_bb_residues)
    set_CG_residues  = set(list_CG_residues)

    set_CG_residues  = set_CG_residues - set_bb_residues - set_AT_residues
    set_bb_residues  = set_bb_residues - set_AT_residues

    list_bb_residues = sorted(list(set_bb_residues))
    list_CG_residues = sorted(list(set_CG_residues)) 

    print("\no List of atomistic, hybrid, and coarse-grained residues computed... 70% completed\n") 


    # 6.7 -  Now, we have the list of all residues and, for each of them, we know wich degree of resolution we have: at, cg (CA) or cg (backbone).
    #        Therefore, starting from each list, we trasform them respectively in list_fullyAT_atoms, list_bb_atoms and list_CG_atoms. 
    #
    #        o  The the first list we take all the atoms with the same res_number 
    # 
    #        o  For the second list, on the other hand, we need to take only the atoms with "at_name" equals to CA, O, N and C (or OC1 in case of 
    #           terminal residues where O is not present, but OC1 and OC2)  
    #         
    #        o  For the third list, we need to retain only the CA, namely when the at_name is equals to CA. 

    list_AT_atoms = [j[3] for i in list_AT_residues for j in atoms if(j[0] == i)] 
    list_bb_atoms = [j[3] for i in list_bb_residues for j in atoms if((j[0] == i) and (j[2]=="CA" or j[2]=="O" or j[2]=="N" or j[2]=="C" or j[2]=="OC1"))]
    list_CG_atoms = [j[3] for i in list_CG_residues for j in atoms if((j[0] == i) and (j[2] == "CA"))]

    # 6.8  - Let us write the file "list-atoms-at-cg.txt" containing the list of survived atoms. 
    #        In particular the atomistic atoms requires the caption "at"; on the other hand, backbone and cg atoms requires the caption "cg" 

    print("\no Writing file containing the list of survived atoms... 90% completed\n") 

    f_list = open("list-atoms-opt-1.txt", "w+")

    for i in list_AT_atoms:
        f_list.write("{:6d}   {:2s}\n".format(i, "at"))      # j[3] is the column relative of the atomistic number.
                           	                             # "at" is the corresponding caption  
    for i in list_bb_atoms: 
        f_list.write("{:6d}   {:2s}\n".format(i, "cg"))      # j[3] is the column relative of the atomistic number.
                                                             # "cg" is the corresponding caption
    for i in list_CG_atoms:
        f_list.write("{:6d}   {:2s}\n".format(i, "cg"))      # j[3] is the column relative of the atomistic number.
                                                             # "cg" is the corresponding caption
    f_list.close()

    os.system("sort -n -o temp list-atoms-opt-1.txt")        # "sort -n" is better than "sort -t" that it does not work properly 
    os.system("mv temp list-atoms-opt-1.txt")

    print("\no '{}' containing the list of survived atoms has been written... 95% completed.\n".format(f_list.name))

    print("\no No Errors... 100% completed\n")

## --------------------------------------------------------------------------------------------------------------------------------------------------


#7. Choice2 selected...  
elif(sys.argv[1].strip() == "choice2"):

    # 7.1 - Checking if grofile is present 
    if not os.path.isfile(grofile):
        print("Error while opening the file. '{}' does not exist\n".format(grofile))
        quit()

    gro_f    = open(grofile, "r")
    check_empty_file(grofile)

    gro      = readgro_all(gro_f)

    atoms    = gro[0]
    n_atoms  = gro[1]

    print("\no '{}' correctly read... 10% completed.\n".format(os.path.basename(grofile)))

    # 7.2 - Opening and reading the file containing the list of residue(s) (not atoms), that will be described atomistically 

    if not os.path.isfile(listfile):
       print("Error while opening the file. '{}' does not exist\n".format(listfile))
       quit()

    f = open(listfile, "r")

    list_AT_residues = []

    count_line = 1 
    for line in f:
        check_empty_rows_block(line, listfile)
        line = line.split() 
 
        if(len(line) != 1): 
            print("\n*Choice2* selected: the input file '{}' cannot be read".format(listfile))
            print("\nError. Each line must contain only 1 column: the {}th row has {} columns".format(count_line, len(line)))
            print("\nBe sure that '{}' is organized in only one column containing the list of residue numbers modelled atomistically\n".format(listfile))
            print("\nLook below for further information\n")
            print_help()
            quit()
      
        else: 
            if not line[0].isdigit():
                print("*choice2* selected: the input file '{}' cannot be read correctly.".format(listfile))
                print("\nError. The unique column must contain an integer number. '{}' at {}th row is not". format(line[0], count_line))
                print("\nBe sure that '{}' is organized in only one column containing the list of residue numbers modelled atomistically\n".format(listfile))
                print("\nLook below for further information\n")
                print_help()
                quit()

        
        res = int(line[0]) 
        list_AT_residues.append(res)

        count_line = count_line + 1
    
    print("\no '{}' correctly read... 20% completed.\n".format(os.path.basename(listfile)))
      
    # 7.3 - Creation of a list corresponding at all the atoms described atomistically, taken from the list of the atomistic residues (list_AT_residues)  
    #       and the list of the other atoms (list_other_atoms)

    list_AT_atoms = [j for i in list_AT_residues for j in atoms if(j[0] == i)]

    list_other_atoms = [x for x in atoms if not x in list_AT_atoms]

    # 7.4 - Computing the distance between the  atoms described atomistically (list_AT_atoms) and the other atoms (list_other_atoms): 
    #
    #       If distance  < bb Diameter (D)    => list_bb_residues =  list of residues that require an hybrid description 
    #
    #       If distance  >= bb Diameter (D)    => list_CG_residue = list of residues that require a CG description    
    #

    if(diameter is not None):  
        D = diameter 
        if(D < 0): 
            print("A negative value of diameter do not make any sense. Change it.\n\n")
            print_help_block()
            quit() 

    else:   
        D = 1.0   # 1 nm: default value of the diameter of the transition area modelled with only the backbone atoms.  

    print("\no Diameter D of hybrid region: {} nm... 30% completed\n".format(D))

    D = float(D)

    print("\no Computing distances... 35% completed\n")

    list_bb_residues = []
    list_CG_residues = []

    for i in list_other_atoms: 
        for j in list_AT_atoms: 

            x   = j[4]-i[4]
            y   = j[5]-i[5]
            z   = j[6]-i[6]

            dist = sqrt(x**2 + y**2 + z**2)

            if(dist < D):
                list_bb_residues.append(i[0]) 

            if(dist >=D):
                list_CG_residues.append(i[0]) 

    # 7.5 - Removing repetitions because some atoms have the same residue number and sorting out the three list in ascending order.

    list_AT_residues = list(set(list_AT_residues))
    list_bb_residues = list(set(list_bb_residues))
    list_CG_residues = list(set(list_CG_residues))

    list_AT_residues.sort()
    list_bb_residues.sort()
    list_CG_residues.sort()


    # 7.6  - Summarizing, the following list of residues have been constructed: 
    #
    #        list_AT_residues            => It contains the list of residues that has to be modelled atomistically,
    #                                       therefore all the atoms of a specific amino acid will be retained (they survive) 
    #
    #        list_bb_residues            => It contains the list of residues that has to be modelled as CG but retaining 
    #                                       only CA, C, O and N of a specific amino acid. 
    #
    #        list_CG_residues            => It contains the list of residue that has to be modelled as CG, but at difference 
    #                                       with the previous case, we retain only the CA of a specific amino acid.  
    #
    #        Assuming that two or more resolutions (AT, bb or CG) overlap each other, the higher resolution is kept: 
    #
    #        AT > BB > CG 
    #
    #        In order to do it easily, we transform the previous three list as set: in this way, than by applying just the difference 
    #        between sets, it will be possible to know only the residues that will be present ONLY in each set. The operation mentioned 
    #        between sets is, indeed, the difference. In other words: 
    #
    #        set_CG_residues = set_CG_residues - set_bb_residues - set_AT_residues   => this is the set of residues that belongs to the 
    #                                                                                   CG only, and not at AT and backbone 
    #
    #        set_bb_residues = set_bb_residues - set_AT_residues                     => this is the set of residues that belongs to the 
    #                                                                                   backbone only, and not at AT. 
    #

    set_AT_residues  = set(list_AT_residues)
    set_bb_residues  = set(list_bb_residues)
    set_CG_residues  = set(list_CG_residues)

    set_CG_residues  = set_CG_residues - set_bb_residues - set_AT_residues
    set_bb_residues  = set_bb_residues - set_AT_residues

    list_bb_residues = sorted(list(set_bb_residues))
    list_CG_residues = sorted(list(set_CG_residues))

    print("\no List of atomistic, hybrid, and coarse-grained residues computed... 70% completed\n")

    # 7.7 -  Trasforming each list of residues in list of atoms (list_fullyAT_atoms, list_bb_atoms and list_CG_atoms). 
    #
    #        o  list_AT_atoms =>  all the atoms with the same res_number 
    #
    #        o  list_bb_atoms =>  backbone atoms, i.e. with "at_name" equals to CA, O, N and C (or OC1 in case of 
    #                             terminal residues where O is not present, but OC1 and OC2)
    #
    #        o  list_CG_atoms =>  only C-alpha atoms, i.e. with "at_name" equals to CA. 

    list_AT_atoms = [j[3] for i in list_AT_residues for j in atoms if(j[0] == i)]
    list_bb_atoms = [j[3] for i in list_bb_residues for j in atoms if((j[0]==i) and (j[2]=="CA" or j[2]=="O" or j[2]=="N" or j[2]=="C" or j[2]=="OC1"))]
    list_CG_atoms = [j[3] for i in list_CG_residues for j in atoms if((j[0] == i) and (j[2] == "CA"))]


    # 7.8  - Writing file "list-atoms-at-cg.txt" containing the list of survived atoms: 
    #        the atomistic atoms requires the caption "at"; on the other hand, backbone and cg atoms requires the caption "cg". 

    print("\no Writing file containing the list of survived atoms... 90% completed\n")

    f_list = open("list-atoms-opt-2.txt", "w+")

    for i in list_AT_atoms:
        f_list.write("{:6d}   {:2s}\n".format(i, "at"))      # j[3] is the column relative of the atomistic number.
                                                             # "at" is the corresponding caption  
    for i in list_bb_atoms:
        f_list.write("{:6d}   {:2s}\n".format(i, "cg"))      # j[3] is the column relative of the atomistic number.
                                                             # "cg" is the corresponding caption
    for i in list_CG_atoms:
        f_list.write("{:6d}   {:2s}\n".format(i, "cg"))      # j[3] is the column relative of the atomistic number.
                                                             # "cg" is the corresponding caption
    f_list.close()

    os.system("sort -n -o temp list-atoms-opt-2.txt")        # "sort -n" is better than "sort -t" that it does not work properly 
    os.system("mv temp list-atoms-opt-2.txt")

    print("\no '{}' containing the list of survived atoms has been written... 95% completed.\n".format(f_list.name))

    print("\no No Errors... 100% completed\n")


## -----------------------------------------------------------------------------------------------------------------------------------


# 8. Choice3 selected...
elif(sys.argv[1].strip() == "choice3"):


    #8.1 - Opening grofile and reading it.
    if not os.path.isfile(grofile):
        print("Error while opening the file. '{}' does not exist\n".format(grofile))
        quit()

    gro_f    = open(grofile, "r")
    check_empty_file(grofile)

    gro      = readgro_all(gro_f)

    atoms    = gro[0]
    n_atoms  = gro[1]

    print("\no '{}' correctly read... 10% completed\n".format(os.path.basename(grofile)))

    # 8.2 - Opening and reading the file containing the list of residue(s) (not atoms) in the fullyAT (1st column) part 
    #       and the hybrid one (2nd column) 

    if not os.path.isfile(listfile):
       print("Error while opening the file. '{}' does not exist\n".format(listfile))
       quit()

    f = open(listfile, "r")

    list_AT_residues = []
    list_bb_residues = [] 

    count_line = 1 
    for line in f:
        check_empty_rows_block(line, listfile)
        line = line.split()     # number of columns must be == 2 
				# if Ncols == 2, both columns must contain integer numbers (isdigit check for int numbers, not float) 
				# If Ncols == 1, the column must contain integer numbers 
 
        if(len(line) < 1 or len(line) > 2): 
            print("\n*Choice3* selected: the input file '{}' cannot be read".format(listfile))
            print("\nError. Each line must contain 1 or 2 columns: the {}th row has {} columns:".format(count_line, len(line)))
            print("\nBe sure that each row of '{}' is organized in only one or two columns\n") 
            print("\nLook below for further information\n")
            print_help_block()
            quit()

        else: 
            if(len(line) == 1):
                if not line[0].isdigit():
                    print("*choice3* selected: the input file '{}' cannot be read correctly.".format(listfile))
                    print("\nError. The 1st column must be an integer number. '{}' at {}th row is not". format(line[0], count_line))
                    print("\nBe sure that each row of '{}' is organized in one or two columns:\n".format(listfile))
                    print(" • 1st column       INT                List of residue numbers modelled atomistically")
                    print(" • 2nd column       INT                List of residue numbers that require an hybrid treatment") 
                    print("\nLook below for further information\n")
                    print_help_block()
                    quit()
                else: 
                    AT_res  = int(line[0])
                    list_AT_residues.append(AT_res)
                    count_line = count_line + 1 

            if(len(line) == 2):
                if not line[1].isdigit():
                    print("*choice3* selected: the input file '{}' cannot be read correctly.".format(listfile))
                    print("\nError. The 2nd column must be an integer number. '{}' at {}th row is not". format(line[1], count_line))
                    print("\nBe sure that each row of '{}' is organized in one or two columns:\n".format(listfile))
                    print(" • 1st column       INT                List of residue numbers modelled atomistically")
                    print(" • 2nd column       INT                List of residue numbers that require an hybrid treatment")
                    print("\nLook below for further information\n")
                    print_help_block()
                    quit()
                else: 
                    AT_res  = int(line[0])
                    bb_res  = int(line[1])
                    list_AT_residues.append(AT_res)
                    list_bb_residues.append(bb_res)
                    count_line = count_line + 1 


    print("\no '{}' correctly read... 30% completed\n".format(os.path.basename(listfile)))    
    
    # 8.3 - Creation of a list corresponding at the residues modelled as CG beads (list_CG_residues)  
    #       as we know both the list of AT and bb residues without repetitions, and we remove repetitions (using SET) and sort it. 
 
    list_CG_residues = [x[0] for x in atoms if (not x[0] in list_AT_residues and not x[0] in list_bb_residues)]

    set_CG_residues  = set(list_CG_residues)

    list_CG_residues = sorted(list(set_CG_residues)) 


    print("\no List of atomistic, hybrid, and coarse-grained residues computed... 70% completed\n")

    # 8.4 -  Trasforming each list of residues in list of atoms (list_fullyAT_atoms, list_bb_atoms and list_CG_atoms). 
    #
    #        o  list_AT_atoms =>  all the atoms with the same res_number 
    #
    #        o  list_bb_atoms =>  backbone atoms, i.e. with "at_name" equals to CA, O, N and C (or OC1 in case of 
    #                             terminal residues where O is not present, but OC1 and OC2) 
    #
    #        o  list_CG_atoms =>  only C-alpha atoms, i.e. with "at_name" equals to CA. 

    list_AT_atoms = [j[3] for i in list_AT_residues for j in atoms if(j[0] == i)]
    list_bb_atoms = [j[3] for i in list_bb_residues for j in atoms if((j[0]==i) and (j[2]=="CA" or j[2]=="O" or j[2]=="N" or j[2]=="C" or j[2]=="OC1"))]
    list_CG_atoms = [j[3] for i in list_CG_residues for j in atoms if((j[0] == i) and (j[2] == "CA"))]   


    # 8.5 -  Let us write the file "list-atoms-at-cg.txt" containing the list of survived atoms. 
    #        In particular the atomistic atoms requires the caption "at"; on the other hand, backbone and cg atoms requires the caption "cg" 

    print("\no Writing file containing the list of survived atoms... 90% completed\n")

    f_list = open("list-atoms-opt-3.txt", "w+")

    for i in list_AT_atoms:
        f_list.write("{:6d}   {:2s}\n".format(i, "at"))      # j[3] is the column relative of the atomistic number.
                                                             # "at" is the corresponding caption  
    for i in list_bb_atoms:
        f_list.write("{:6d}   {:2s}\n".format(i, "cg"))      # j[3] is the column relative of the atomistic number.
                                                             # "cg" is the corresponding caption
    for i in list_CG_atoms:
        f_list.write("{:6d}   {:2s}\n".format(i, "cg"))      # j[3] is the column relative of the atomistic number.
                                                             # "cg" is the corresponding caption
    f_list.close()

    os.system("sort -n -o temp list-atoms-opt-3.txt")        # "sort -n" is better than "sort -t" that it does not work properly 
    os.system("mv temp list-atoms-opt-3.txt")

    print("\no '{}' containing the list of survived atoms has been written... 95% completed.\n".format(f_list.name))

    print("\no No Errors... 100% completed\n")   

REAL_end = datetime.now()
print("The time for executing the code(BLOCK) is: ", (REAL_end - REAL_start).total_seconds()) 
