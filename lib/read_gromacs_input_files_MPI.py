import os 
from inp_out import *
from read_gromacs import * 

# FUNCTION: Reading Coordinate input file 
def read_coordinate_grofile(grofile):

    ### Checking if grofile is present  
    if not os.path.isfile(grofile):
        print("Error while opening the file. '{}' does not exist\n".format(grofile))
        quit()

    ### Opening and reading grofile 
    f = open(grofile, "r")
    check_empty_file(grofile)

    groo = readgro_all(f)

    print("\nOriginal coordinate file correctly read... 2% completed\n")
    
    return groo


# FUNCTION: Reading List of Survived Atoms
def read_list_survived_atoms(listatomsfile):

    ### Checking if listatomsfile is present 
    if not os.path.isfile(listatomsfile):
        print("Error while opening the file. '{}' does not exist\n".format(listatomsfile))
        quit()
        
    ### Opening and reading the file, and printing errors 
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
        if(k<0 or k>3000000):#n_atoms):
            print("Error. The list of survived atoms '{}' must have index between 1 and {}".format(listatomsfile, 300000))#n_atoms))
            print("\nLook below for further information\n")
            print_help_main_CANVAS()
            quit()
    
    f_list_survived.seek(0)                        # read again file if everything went fine.
    n_lines_f = sum(1 for line in f_list_survived) # compute n lines of listatomsfile
    
    if(n_lines_f != len(dict_survived)):
        print("Error. The list of survived atoms '{}' contains, at least, two identical indexes at 1st column".format(listatomsfile))#Avoid repeated elements
        print("\nLook below for further information\n")
        print_help_main_CANVAS()  
        quit()
        
    ### Number of survived atoms  
    n_survived  = len(dict_survived)     
        
    print("\nlist of survived atoms correctly read... 4% completed\n")

    return dict_survived, n_survived
    #l1.put((dict_survived, n_survived))



# FUNCTION: Reading Topology File 
def read_TopologyFile(topfile):  # (topofile, t1)  

    ### Checking if topfile is present  
    if not os.path.isfile(topfile):
        print("Error while opening the file. '{}' does not exist\n".format(topfile))
        quit()
    
    ### Opening and reading topology file
    ft = open(topfile, "r")
    check_empty_file(topfile)
    
    toop = readtop(ft)
    
    original_cmaps_protein_fully_at      = toop[7]  # In case of amber-forcefield, this list is empty
        
    Flag_cmaps = True
    if not original_cmaps_protein_fully_at: # If CMAP list is empty:
       Flag_cmaps = False
    
    print("\nFully-atomistic topology file topol.top correctly read... 5% completed\n")

    return toop, Flag_cmaps
