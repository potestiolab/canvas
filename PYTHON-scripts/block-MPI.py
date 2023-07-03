"""
The scope of this program is to create a file containing the list of survived atoms each having having the label "at" or "cg". The former is used for those atoms that will be treated atomistically, whereas "cg" label is used for those atoms that will be treat as bead (Medium-Grained or Coarse-Grained). 
"""

                   #################################################
#####################  1. Importing libraries, parsing arguments  ####################################
                   #################################################


# 1.1 Importing main libraries
import sys 
import re 
import argparse
import itertools 
import os 
import itertools   		      # It needs for removing duplicates in list of lists.

from operator import itemgetter
from math import sqrt 

from datetime import datetime

from functools import partial
from multiprocessing import Process, Pool, cpu_count, Queue


# --------------------------------------------

REAL_start = datetime.now() 


# 1.2 Finding the path of the  main folder (usually "canvas") after searching for 'PYTHON-scripts' folder. 
#     Then, add /lib in order to find our libraries. 

desired_folder_name = "PYTHON-scripts"
current_directory = os.getcwd()
desired_path = None

while True:
    if desired_folder_name in os.listdir(current_directory):
        desired_path = current_directory
        break
    elif current_directory == os.path.dirname(current_directory):
        print("ERROR. 'PYTHON-script' folder has not been found. Please, check it out...\n")
        quit()
    else:
        current_directory = os.path.dirname(current_directory)

python_modules_path = desired_path + "/lib"
sys.path.append(python_modules_path)


# 1.3 Importing user-libraries 
from read_gromacs import * 
from inp_out import * 


# 1.4. Input Arguments -------------------------------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, add_help=False) 

group_in=parser.add_argument_group("Required Arguments") 

group_in.add_argument('option', help = argparse.SUPPRESS)					                          # Mandatory
group_in.add_argument('-g', '--gro', dest='grofile', action='store', metavar = 'FILE', help = argparse.SUPPRESS)          # Mandatory
group_in.add_argument('-l', '--list', dest='listfile', metavar = 'FILE',help = argparse.SUPPRESS)                         # Mandatory 

group_in.add_argument('-h', '--help', action='help', help = argparse.SUPPRESS)                                            # Optional
group_in.add_argument('-D', '--diameter', dest='DiameterValue', metavar = 'VALUE', help = argparse.SUPPRESS)              # Optional 
group_in.add_argument('-n', '--ncpu', dest='NumberCpu', metavar='INT', type=int, help = argparse.SUPPRESS)                # Optional 

# ------------------------------------------------------------------------------------------------------------------------------------------



                   ##################################################
#####################  2. Writing Functions employed in this code  ####################################
                   ##################################################


# 2.1 - Knowing the list of residues having atomistic description, the function returns the list of all the atom_numbers that belong to those residues. 
def from_ATresnum_to_ATatomnum(list_AT_residues, j):  
    for i in list_AT_residues:
        if(j[0]==i):
            return j[3]

# 2.2 - Knowing the list of residues having medium-grained description, the function returns the list of the atom_numbers that correspond at 
#       CA, O, N, and C (or OC1 in case of terminal residues where O is not present), that are the backbone atoms for each residue. 
def from_MGresnum_to_MGatomnum(list_bb_residues, j):  
    for i in list_bb_residues:
        if((j[0]==i) and (j[2]=="CA" or j[2]=="O" or j[2]=="N" or j[2]=="C" or j[2]=="OC1")):
            return j[3]


# 2.3 - Knowing the list of residues having coarse-grained description, the function returns the list of the atoms_numbers that correspond at 
#       CA, i.e. only the Carbon-alpha atom for each residue. 
def from_CGresnum_to_CGatomnum(list_CG_residues, j):    #update_6p7_cg     #update_7p7_CG     #update_8p4_CG
    for i in list_CG_residues:
        if((j[0]==i) and (j[2] == "CA")):
            return j[3]


# 2.4 - Checking if a residue will described atomistically, medium-grained or coarse-grained. It is done by knowing the diameter (distance) D 
#       with respect all the atoms that are described atomistically. Computing the distance between the atoms described atomistically (list_AT_atoms) 
#       and the other atoms (list_other_atoms) it turns out that:
#
#        If distance  < bb Diameter (D)     => list_bb_residues =  list of residues that require a MG description 
#        If distance  >= bb Diameter (D)    => list_CG_residue  = list of residues that require a CG description
#
#      This function is employed only when 'choice2' task is employed.  
list_bb_residues = []
list_CG_residues = []

def checking_MG_CG_residues_choice2(list_other_atoms, D, j):                   #update_7p4
    for i in list_other_atoms:

        x   = j[4]-i[4]
        y   = j[5]-i[5]
        z   = j[6]-i[6]

        dist = sqrt(x**2 + y**2 + z**2)

        if(dist < D):
            list_bb_residues.append(i[0])

        if(dist >=D):
            list_CG_residues.append(i[0])

    return list_bb_residues, list_CG_residues



# 2.5 - Computing the distance between the atomistic atoms (list_AT_atoms) and the other atoms (list_other_atoms) 
#       and creating a list of atomistic, medium-grained, and coarse-grained residues. 
#   
#         If distance < atomistic_radius             => list_fullyAT_residues = list of residues that require an atomistic description         
#         If  at_radius < distance < bb Diameter (D) => list_bb_residues      = list of residues that require a MG description 
#         If distance > bb Diameter (D)              => list_CG_residue       = list of residues that require a CG description    

def checking_AT_MG_CG_residues_choice1(list_AT_atoms, list_other_atoms, D, list_AT_res):
    list_AT_residues = []
    list_bb_residues = []
    list_CG_residues = []

    count = 0


    ### Avoid creating MPI function because the presence of "count" creates issues 
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

        if(i[2] == "O" or i[2] == "OC2"):                             # It means that the next residue has been reached (OC2 in case of terminal res) 
            count = count + 1                                         # then, count can be increased by one 


    list_AT_residues = list_AT_residues + list_AT_res  # Adding in list_AT_residue the atomistic residues writtein in file.   

    return list_AT_residues, list_bb_residues, list_CG_residues 


                   ##################################################
#####################  3. Checking Errors in this code  ####################################
                   ##################################################


if __name__ == '__main__':

    # 3.1 Printing on terminal that this code is running 
    print("\n####################################################\n")
    print("'{}' running...".format(os.path.basename(sys.argv[0]))) 
    print("---------------------------------\n")
    
    # 3.2 Printing help message if the script does not have any arguments  
    if len(sys.argv)==1:
        print_usage_main_block()
        quit()
    

    # 3.3 Printing help message if the script does not present valid arguments: in particular, the first argument must be 'choice1', 'choice2' or 'choice3' 
    checking_accepted_tasks_block()
    check_argv_errors_block()


    # 3.4 Printing error and help message if code presents not allowed arguments
    checking_valid_arguments_block(parser)

   
    # 3.5 Parsing arguments
    args     = parser.parse_args()

    option   = args.option           # Mandatory 
    grofile  = args.grofile          # Mandatory
    listfile = args.listfile         # Mandatory    
    diameter = args.DiameterValue    # Optional 
    ncpu     = args.NumberCpu        # Optional 

    
    # 3.6 Checking if mandatory files are present 
    mandatory_files_present_block(grofile, listfile) 


    # 3.7 Printing, on terminal which option has been set among 'choice1', 'choice2' or 'choice3'.
    print("\no '{}' option set. 10% completed.\n".format(sys.argv[1]))

   
    # 3.8 Checking if grofile (coordinate file) is actually found and that it is not empty 
    checking_file_found_block(grofile)
    check_empty_file(grofile)

    
    # 3.9 Reading grofile
    gro_f = open(grofile, "r")

    gro      = readgro_all(gro_f)
    atoms    = gro[0]
    n_atoms  = gro[1]

    print("\no '-g/--gro {}' set. Coordinate file correctly read... 15% completed.\n".format(os.path.basename(grofile)))  # Print ONLY Filename


    # 3.10 Checking if listfile  is actually found and that it is not empty (for choice1 or choice2). 
    #      Indeed, when using choice3, listfile could also be empty: in such case the biomolecules is modelled all CG.
    checking_file_found_block(listfile)

    if(sys.argv[1].strip() == "choice1" or sys.argv[1].strip() == "choice2"):
        check_empty_file(listfile)


    # 3.11 Reading listfile and checking that the format is correct according with the choice done (choice1, choice2, or choice3):
    #
    #      If 'choice1' --> listfile is organized in two columns: 
    #                        o 1st: list of atomistic residue(s) (INT) that correspond at the center of atomistic region from which an atomistic sphere is traced.    #                        o 2nd: radius (INT/FLOAT) of each atomistic sphere traced.  
    #                   
    #      If 'choice2' --> listfile is organized is one column corresponding at the list of all atomistic residues (INT) that will be modelled atomistically. 
    #
    #      If 'choice3' --> listfile is organized in two columns:
    #                        o 1st: list of all atomistic residues (INT) that will be modelled atomistically. 
    #                        o 2nd: list of all medium-grained residues (INT) that will be modelled as medium-grained. 
    

    if(sys.argv[1].strip() == "choice1"):
        list_choice1 = checking_errors_listfile_choice1_mand_arg(listfile)
        list_AT_res  = list_choice1[0]
        list_radii   = list_choice1[1]        

    if(sys.argv[1].strip() == "choice2"):
        list_AT_residues = checking_errors_listfile_choice2_mand_arg(listfile)

    if(sys.argv[1].strip() == "choice3"):
        list_choice3     = checking_errors_listfile_choice3_mand_arg(listfile)
        list_AT_residues = list_choice3[0]
        list_bb_residues = list_choice3[1]
    
    print("\no '-l/--list {}' set. listfile for {} task correctly read... 25% completed\n".format(os.path.basename(listfile), sys.argv[1]))

    # 3.12 Checking if the optional argument 'diameter' is set. An error occurs if a string or negative (incl. zero) is inserted. 
    #      Only an integer or float number is accepted. Moreover, this argument is accepted only if 'choice1' or 'choice2' option is set.  
    D = checking_errors_diameter_opt_arg(diameter)  


    # 3.13 Checking if the optional argument 'ncpu' is set. The program returns an error if the user askes for a number of cores 
    #      higher than the maximum allowed. If 'ncpu' is not specified, the maximum number of cores is employed. 
    #      This argument can be set whatever the choice (choice1, choice2, or choice3).  
    ncpu_employed = checking_errors_ncpu_opt_arg(ncpu)  



                  ############################
#####################  4. CHOICE1 selected  ####################################
                   ###########################
   
    
    
    # 4. Choice1 selected... 
    if(option == "choice1"):
    

        # 4.1 - Creation of a list corresponding at the C-alpha atoms, kept from the list of the atomistic residues (list_AT_residues)        
        list_AT_atoms    = [j for i in list_AT_res for j in atoms if (j[0] == i)]
        list_other_atoms = [x for x in atoms if not x in list_AT_atoms]
    

        # 4.2 - Creating list of atomistic, MG, and CG residues according with a distance criterium ('checking_AT_MG_CG_residues_choice1' function)          
        res_list         = checking_AT_MG_CG_residues_choice1(list_AT_atoms, list_other_atoms, D, list_AT_res)
        list_AT_residues = res_list[0]
        list_bb_residues = res_list[1]
        list_CG_residues = res_list[2]    

    
        # 4.3 - Removing repetitions because some atoms have the same residue number and sorting out the three list in ascending order.
        list_AT_residues = list(set(list_AT_residues))
        list_bb_residues = list(set(list_bb_residues))
        list_CG_residues = list(set(list_CG_residues))
    
        list_AT_residues.sort() 
        list_bb_residues.sort() 
        list_CG_residues.sort() 
    
        # 4.4 - Summarizing, there are the following list of residues:
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
    
        print("\no List of atomistic, medium-grained, and coarse-grained residues computed... 70% completed\n") 
    
        # 4.5 -  The list of all residues and relative resolution is known (list_AT_residues, list_bb_residues, list_CG_residues). 
        #        Thus, the next step is to create the lists of the the atoms that will be retained for each residue according with its resolution. 
        #
        #        o list_AT_residues <--> list_AT_atoms: all the atoms will be reatined 
        #        o list_bb_residues <--> list_bb_atoms: all the backbone atoms will be retained. 
        #        o list_CG_residues <--> list_CG_atoms: only the CA atom will be retained.  


        if(ncpu is not None): 
            pool       = Pool(ncpu)
        else:
            pool = Pool() 
        

        temp       = partial(from_ATresnum_to_ATatomnum, list_AT_residues)  
        upd_6p7_AT = pool.map(temp, iterable=atoms) 
        
        temp       = partial(from_MGresnum_to_MGatomnum, list_bb_residues)  
        upd_6p7_BB = pool.map(temp, iterable=atoms)
        
        temp       = partial(from_CGresnum_to_CGatomnum, list_CG_residues) 
        upd_6p7_CG = pool.map(temp, iterable=atoms)
    
        pool.close()
        pool.join() 
    
        list_AT_atoms = [i for i in upd_6p7_AT if i is not None]
        list_bb_atoms = [i for i in upd_6p7_BB if i is not None]
        list_CG_atoms = [i for i in upd_6p7_CG if i is not None]
    
    
        # 4.6 - Let us write the file "list-atoms-at-cg.txt" containing the list of survived atoms. 
        #       In particular the atomistic atoms requires the caption "at"; on the other hand, backbone and cg atoms requires the caption "cg" 
    
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
    

                  ############################
#####################  5. CHOICE2 selected  ####################################
                   ###########################

    
    elif(sys.argv[1].strip() == "choice2"):
    
    
        # 5.1 - Creation of a list corresponding at the C-alpha atoms, kept from the list of the atomistic residues (list_AT_residues)  
        list_AT_atoms    = [j for i in list_AT_residues for j in atoms if(j[0] == i)]
        list_other_atoms = [x for x in atoms if not x in list_AT_atoms]
    
   
        # 5.2 - Creating list of atomistic, MG, and CG residues according with a distance criterium ('checking_AT_MG_CG_residues_choice2' function) (MPI)  
        if(ncpu is not None):        
            pool = Pool(ncpu)
        else:
            pool = Pool()
    
        temp    = partial(checking_MG_CG_residues_choice2, list_other_atoms, D) 
        upd_7p4 = pool.map(temp, iterable=list_AT_atoms)
    
        pool.close()
        pool.join()
    
        list_bb_residues = [i[0] for i in upd_7p4]
        list_CG_residues = [i[1] for i in upd_7p4] 
    
        list_bb_residues = list(itertools.chain(*list_bb_residues))  #merge list of lists in a unique list
        list_CG_residues = list(itertools.chain(*list_CG_residues))
    
    

        # 5.3 - Removing repetitions because some atoms have the same residue number and sorting out the three list in ascending order.
    
        list_AT_residues = list(set(list_AT_residues))
        list_bb_residues = list(set(list_bb_residues))
        list_CG_residues = list(set(list_CG_residues))
    
        list_AT_residues.sort()
        list_bb_residues.sort()
        list_CG_residues.sort()
    
    
        # 5.4  - Summarizing, the following list of residues have been constructed: 
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
    
        print("\no List of atomistic, medium-grained, and coarse-grained residues computed... 70% completed\n")


        # 5.5 -  The list of all residues and relative resolution is known (list_AT_residues, list_bb_residues, list_CG_residues). 
        #        Thus, the next step is to create the lists of the the atoms that will be retained for each residue according with its resolution. 
        #
        #        o list_AT_residues <--> list_AT_atoms: all the atoms will be reatined 
        #        o list_bb_residues <--> list_bb_atoms: all the backbone atoms will be retained. 
        #        o list_CG_residues <--> list_CG_atoms: only the CA atom will be retained.  

    
        if(ncpu is not None):        
            pool = Pool(ncpu)
        else:
            pool = Pool()


        temp       = partial(from_ATresnum_to_ATatomnum, list_AT_residues)  #OLD --> #temp = partial(update_7p7_AT, list_AT_residues)
        upd_7p7_AT = pool.map(temp, iterable=atoms)

        temp       = partial(from_MGresnum_to_MGatomnum, list_bb_residues)  #OLD -->  #temp = partial(update_7p7_BB, list_bb_residues)
        upd_7p7_BB = pool.map(temp, iterable=atoms)

        temp       = partial(from_CGresnum_to_CGatomnum, list_CG_residues)  #OLD -->  #temp = partial(update_7p7_CG, list_CG_residues)
        upd_7p7_CG = pool.map(temp, iterable=atoms) 

    
        pool.close()
        pool.join()
    
        list_AT_atoms = [i for i in upd_7p7_AT if i is not None]
        list_bb_atoms = [i for i in upd_7p7_BB if i is not None]
        list_CG_atoms = [i for i in upd_7p7_CG if i is not None]
    
    
        # 5.6  - Writing file "list-atoms-at-cg.txt" containing the list of survived atoms: 
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
    
    
                  ############################
#####################  6. CHOICE3 selected  ####################################
                   ###########################
    
    
    elif(sys.argv[1].strip() == "choice3"):
    
    
        
        # 6.1 - Creation of a list corresponding at the residues modelled as CG beads (list_CG_residues)  
        #       as we know both the list of AT and bb residues without repetitions, and we remove repetitions (using SET) and sort it. 
        
        list_CG_residues = [x[0] for x in atoms if (not x[0] in list_AT_residues and not x[0] in list_bb_residues)]  
        set_CG_residues  = set(list_CG_residues)
        list_CG_residues = sorted(list(set_CG_residues)) 
    
    
        print("\no List of atomistic, medium-grained, and coarse-grained residues computed... 70% completed\n")


        # 6.2 -  The list of all residues and relative resolution is known (list_AT_residues, list_bb_residues, list_CG_residues). 
        #        Thus, the next step is to create the lists of the the atoms that will be retained for each residue according with its resolution. 
        #
        #        o list_AT_residues <--> list_AT_atoms: all the atoms will be reatined 
        #        o list_bb_residues <--> list_bb_atoms: all the backbone atoms will be retained. 
        #        o list_CG_residues <--> list_CG_atoms: only the CA atom will be retained.  


        if(ncpu is not None):        
            pool       = Pool(ncpu)
        else:
            pool = Pool()
    
        temp       = partial(from_ATresnum_to_ATatomnum, list_AT_residues)     #OLD --> #temp = partial(update_8p4_AT, list_AT_residues)
        upd_8p4_AT = pool.map(temp, iterable=atoms)
    
        temp       = partial(from_MGresnum_to_MGatomnum, list_bb_residues)     #OLD --> #temp = partial(update_8p4_BB, list_bb_residues)
        upd_8p4_BB = pool.map(temp, iterable=atoms)
    
        temp       = partial(from_CGresnum_to_CGatomnum, list_CG_residues)     #OLD --> #temp = partial(update_8p4_CG, list_CG_residues)
        upd_8p4_CG = pool.map(temp, iterable=atoms)
    
        pool.close()
        pool.join()
    
        list_AT_atoms = [i for i in upd_8p4_AT if i is not None]
        list_bb_atoms = [i for i in upd_8p4_BB if i is not None]
        list_CG_atoms = [i for i in upd_8p4_CG if i is not None]
    
    
        # 6.3 -  Let us write the file "list-atoms-at-cg.txt" containing the list of survived atoms. 
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
    
