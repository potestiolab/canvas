import os 
import sys

# Function: it writes a file containing a list in sequence of the CA atom number (it can be useful when using VMD) 
def write_CA(f_CA, CA_list): 
    for i in CA_list:
        f_CA.write("{} ".format(i))   # print CA numbers on the same line 
     

# Function: it checks if X is Integer or Float (Note that isdigit returns True only if X is Integer) 
def isIntFloat(x):
    try:
        float(x)
        return True
    except ValueError:
        return False

# Function: it retuns the number of columns 
def columns(line):
    line = line.split()
    lenght = len(line)

    return lenght


# Function: check if mandatory files are present 
def mandatory_files_present_CANVAS(grofile, listatomsfile, topfile):
    if (grofile is None):
        print("\nError. The Coordinate file is missing. Look below for further help.\n")
        print_help_main_CANVAS()
        quit()

    if (listatomsfile is None):
        print("\nError. The file containing the list of survived atoms is missing. Look below for further help.\n")
        print_help_main_CANVAS()
        quit()

    if (topfile is None): 
        print("\nError. The topology file is missing. Look below for further help. \n") 
        print_help_main_CANVAS()
        quit() 

# Function: check if mandatory files are present 
def mandatory_files_present_block(grofile, listfile): 
    if (grofile is None):
        print("\nError. The Coordinate file missing. Look below for further help.\n")
        print_help_block()
        quit() 

    if (listfile is None):
        print("\nError. The List AT-bb-CG is missing. Look below for further help.\n")
        print_help_block()
        quit() 



# Function: it checks for not accepted flags for CANVAS.py program
def check_argv_errors_CANVAS():
    for i in range(1, len(sys.argv)):
        if(i%2 != 0): # even arguments 

            if(sys.argv[i][0] == '-' and len(sys.argv[i]) == 1):
                print("Error. '-' is not accepted as flag. Use, for example, '-g' instead of '-'. Look below for further help.\n");
                print_help_main_CANVAS()
                quit()

            if(sys.argv[i][0] == '-' and sys.argv[i][1] != '-' and  len(sys.argv[i]) > 2):
                print("Error. Each flag must contain '-' plus ONLY ONE letter. Example: -g. Look below for further help.\n");
                print_help_main_CANVAS()
                quit()

            if(sys.argv[i][0] == '-' and sys.argv[i][1] == '-' and  len(sys.argv[i]) == 2):
                print_help_main_CANVAS()
                quit()


# Function: it checks for not accepted flags for block.py program
def check_argv_errors_block():
    for i in range(1, len(sys.argv)):
        if(i%2 == 0):

            if(sys.argv[i][0] == '-' and len(sys.argv[i]) == 1):
                print("Error. '-' is not accepted as flag. Use, for example, '-g' instead of '-'. Look below for further help.\n");
                print_help_block()
                quit()

            if(sys.argv[i][0] == '-' and sys.argv[i][1] != '-' and  len(sys.argv[i]) > 2):
                print("Error. Each flag must contain '-' plus ONLY ONE letter. Example: -g. Look below for further help.\n");
                print_help_block()
                quit()

            if(sys.argv[i][0] == '-' and sys.argv[i][1] == '-' and  len(sys.argv[i]) == 2):
                print_help_block()
                quit()


# Function: it checks if the file is empty
def check_empty_file(file):
    if(os.path.getsize(file) == 0):
        print("Error. Your file is empty. Please, fill out {} with significant data or use another file".format(file))
        quit()


# Function: it checks if empty rows are present in a file in CANVAS.py program 
def check_empty_rows_CANVAS(line, file):
    if len(line.strip()) == 0:
        print("In '{}' there is, at least, an empty row. Check it out...\n".format(file))
        print("\nLook below for further information\n")
        print_help_main_CANVAS()
        quit()

# Function: it checks if empty rows are present in a file in block.py program 
def check_empty_rows_block(line, file):     
    if len(line.strip()) == 0:
        print("In '{}' there is, at least, an empty row. Check it out...\n".format(file))
        print("\nLook below for further information\n")
        print_help_block()
        quit() 



# Function: it prints the main usage of CANVAS.py program 
def print_usage_main_CANVAS():
    print("Usage: python3 {} -g <Coordinate FILE> -l <List survived atoms FILE> -t <Topology FILE> [-d <Domain division FILE>] [-r Y] [-n <nCPU>] [-c lammps] [-s n]".format(sys.argv[0]))
    print("   or: python3 {} --gro <Coordinate FILE> --list <List survived atoms FILE> --top <Topology FILE> [--dom <Domain division FILE>] [--resc14 Y] [--ncpu <nCPU>] [--code lammps] [--solvate n]".format(sys.argv[0]))

    print("\nTry python3 {} -h or {} --help for more information.\n".format(sys.argv[0], sys.argv[0]))


# Function: it prints the main usage of block.py program 
def print_usage_main_block(): 
    tasks = ["choice1", "choice2", "choice3"]
   
    for tk in tasks: 
        print("Usage: python3 {} {} [OPTIONS]".format(sys.argv[0], tk))

    print("\nTry python3 {} -h or {} --help for more information.\n".format(sys.argv[0], sys.argv[0]))
 

# Function: it prints the main help of CANVAS.py program 
def print_help_main_CANVAS():  
    print("Usage: python3 {} -g <Coordinate FILE> -l <List survived atoms FILE> -t <Topology FILE> [-d <Domain division FILE>] [-r Y] [-n <nCPU>] [-c lammps] [-s n]".format(sys.argv[0]))
    print("   or: python3 {} --gro <Coordinate FILE> --list <List survived atoms FILE> --top <Topology FILE> [--dom <Domain division FILE>] [--resc14 Y] [--ncpu <nCPU>] [--code lammps] [--solvate n]".format(sys.argv[0]))

    print("\n-----------------------------------------------------------------------------------------------------")

    print("*{}* requires the following inputs:\n".format(sys.argv[0]));
    print("   Coordinate FILE                MANDATORY          File of atom Coordinates in .gro format\n")
    print("   List survived atoms FILE       MANDATORY          File containing the list of survived atoms organized in two columns:")
    print("                                                        -------------------------------------------- ")
    print("                                                        | at_num 1 (int) | caption 1 ('at' or 'cg') |")
    print("                                                        | at_num 2       | caption 2                |")
    print("                                                        | at_num 3       | caption 3                |")
    print("                                                        | at_num 4       | caption 4                |")
    print("                                                        | ....           | ....                     |")
    print("                                                        | at_num N       | caption N                |")
    print("                                                        ---------------------------------------------  \n")
    print("   Topology FILE                  MANDATORY          File containing the topology of fullyAT representation\n") 
    print("   Domain Division FILE           OPTIONAL           File with the list of atoms divided in domains")
    print("                                                     The number of columns is equal to the number of domains")
    print("                                                     Each row contains the at_numbers of the atoms that belong")
    print("                                                     to the same block separated by spaces:") 
    print("                                                        ----------------------------------------------------")
    print("                                                        | AT_Num-1   AT_Num-2   AT_Num-3   .....  AT_Num-N  |")  
    print("                                                        | AT_Num-10  AT_Num-11  AT_Num-12  .....  AT_Num-M  |")
    print("                                                        | .....      .....      .....      .....  .....     |")
    print("                                                        ----------------------------------------------------- \n")
    print("  rescaled14                      OPTIONAL           String containing the word 'Y' or 'y'. Other strings")
    print("                                                     are not allowed. If '-r/--resc14 Y' is set, than all pairs")
    print("                                                     where ONLY one CG bead is kept if in the corresponding")
    print("                                                     all-atom representation that pair is present.")
    print("                                                     Keep attention, as it might create artifacts in MD simulation\n")
    print("  codestring                      OPTIONAL           String containing 'lammps' word. Other strings are not")
    print("                                                     allowed. If '-c/--code lammps' is set, then this program")
    print("                                                     produces the input files needed for simulating the CANVAS")
    print("                                                     model in lammps. If '-c/--code gromacs' is set, or in case")
    print("                                                     the user ignore this flag, this program returns the input files")
    print("                                                     for simulating the CANVAS model in gromacs\n")
    print("  solvatestring                   OPTIONAL           String containing 'n' or 'y' letters. Other string are not")
    print("                                                     allowed. If '-s/--solvate n' is set, then this code")
    print("                                                     do not solvate the system. Use it in case of implicit solvent simulations")
    print("                                                     If '-s/--solvate y' is set, or in case the user ignore this flag,")
    print("                                                     this program solvates and neutralizes (with NA and CL ions) the system\n")     
    print("  NumberCpu                       OPTIONAL           Integer number corresponding at the number of cpus that you would like to employ")
    print("                                                     for creating the CANVAS model of a biomolecule.") 
    print("                                                     If '-n/--ncpu <nCPU>' is set, the code will be parallelized by employing nCPU cores")  
    print("                                                     If '-n/--npu <nCPU>' is not set, the code will be parallelized by employing")
    print("                                                     the maximum number of allowed cores in your laptop/cluster\n") 

    print("-----------------------------------------------------------------------------------------------------");
    print("Hereafter the list of flags:\n");

    print("   -g  --gro                      FILE               Coordinate FILE (gro format)")
    print("   -l  --list                     FILE               List Survived Atoms FILE (txt or dat format)")
    print("   -t  --top                      FILE               Topology FILE") 
    print("  [-d] [--dom]                    FILE               Domain Division FILE (txt or dat format")   
    print("  [-r] [--resc14]                 STR                If string is 'Y', rescaled VdW and Coulomb is applied for CG beads") 
    print("  [-c] [--code]                   STR                If string is 'lammps' the program returns the input files for LAMMPS simulation")  
    print("  [-s] [--solvate]                STR                If string is 'n' the system will sot solvated.") 
    print("  [-n] [--ncpu]                   INT                Integer number corresponding at the number of cores employed for parallelizing the code")  

    print("  [-h] [--help]                                      Give this help list\n")

    print("Report bugs to <raffaele.fiorentini@unitn.it>\n")


# Function: it prints the main help of block.py program
def print_help_main_block(): 
    tasks = ["choice1", "choice2", "choice3"]

    for tk in tasks:
        print("Usage: python3 {} {} [OPTIONS]".format(sys.argv[0], tk))

    print("\n-----------------------------------------------------------------------------------------------------")
    print("Please, choose one of the following tasks:\n")
    print("   *choice1*                       You know one or more central atomistic residues from which the system")
    print("                                   will be divided in three regions. In this specific case an atomistic ")
    print("                                   sphere with radius R (around the central residue) will be traced.")
    print("                                   Then, around the atomistic sphere will be traced a transition area ")
    print("                                   with diameter D in which only the backbone atoms are retained.")
    print("                                   Finally, the remainder will be modelled more coarse-grained,")
    print("                                   where only the C-alpha atoms are kept.\n")
    print("   *choice2*                       You know which which residues will be modelled atomistically.")
    print("                                   Around them will be traced an hybrid region with diameter D ")
    print("                                   where only the backbone atoms are retained. ")
    print("                                   Note that in *choice1* the atomistic region is not defined ")
    print("                                   (you know only one or more central residues) and it will be traced with radius R \n")
    print("   *choice3*                       You know in advance which residues will be modelled atomistically, ")
    print("                                   and which residues will be modelled retaining only the backbone atoms. ")
    print("                                   In this case you do not need to specify any diameter D ")
    print("                                   since the knowledge of AT-bb-CG division is already completed. \n")
    print("----------------------------------------------------------------------------------------------------")
    print("Hereafter the list of OPTIONS:\n");

    print("  -g   --gro         FILE          Coordinate FILE (gro format)")
    print("  -l   --list        FILE          List AT-bb-CG FILE (txt or dat format)")
    print(" [-d] [--diameter]   INT/FLOAT     Value of diameter of hybrid region")
    print(" [-h] [--help]                     Give this help list\n")

    print("Try: python3 {} <TASK> for more information about the mandatory options of a specific task\n".format(sys.argv[0]))

    print("Report bugs to <raffaele.fiorentini@unitn.it>\n")



# Function: it prints the help of each task chosen in "block.py"  
def print_help_block(): 
    if(sys.argv[1]=="choice1" or sys.argv[1]=="choice2"):
        print("usage: python3 {} {} -g <Coordinate FILE> -l <List AT-bb-CG FILE> [-D <diameter hybrid region>] [-n <nCPU>]".format(sys.argv[0], sys.argv[1]))
        print("   or: python3 {} {} --gro <Coordinate FILE> --list <List AT-bb-CG FILE> [-diameter <diameter hybrid region>] [--ncpu <nCPU>]"\
               .format(sys.argv[0],sys.argv[1])) 

    else: 
        print("usage: python3 {} {} -g <Coordinate FILE> -l <List AT-bb-CG FILE> [-n <nCPU>]".format(sys.argv[0], sys.argv[1])) 
        print("   or: python3 {} {} --gro <Coordinate FILE> --list <List AT-bb-CG FILE> [--npu <nCPU>]".format(sys.argv[0], sys.argv[1]))

    print("\n-----------------------------------------------------------------------------------------------------")
    print("The *{}* task requires the following inputs:\n".format(sys.argv[1]));  
    print("   Coordinate FILE             MANDATORY          File of atom Coordinates in .gro format\n")

    if(sys.argv[1] == "choice1"): 
        print("   List At-bb-CG FILE          MANDATORY          File with a list of central atomistic(s) residue(s)")
        print("                                                  and corresponding atomistic(s) radius(ii) organized in two columns:")
        print("                                                  -------------------------------------------  ")
        print("                                                  | residue1 (int)  |   radius1 (int/float) |  ")
        print("                                                  | residue2        |   radius2             |  ")
        print("                                                  | ........        |   .......             |  ")
        print("                                                  | residueN        |   radiusN             |  ")
        print("                                                  -------------------------------------------  \n")

        print("  [Diameter hybrid region]     OPTIONAL           Value (in nm) of diameter of the hybrid region")
        print("                                                  Default value: 1.0 nm\n")
        print("  [NumberCpu]                  OPTIONAL           Integer number corresponding at the number of cpus that you would like to employ")
        print("                                                  for creating the CANVAS model of a biomolecule.")
        print("                                                  If '-n/--ncpu <nCPU>' is set, the code will be parallelized by employing nCPU cores")
        print("                                                  If '-n/--npu <nCPU>' is not set, the code will be parallelized by employing")
        print("                                                  the maximum number of allowed cores in your laptop/cluster\n")

            
    if(sys.argv[1] == "choice2"):
        print("   List At-bb-CG FILE          MANDATORY          File organized in only one column")
        print("                                                  that contain the list of atomistic residues:")
        print("                                                  ------------------  ")
        print("                                                  | residue1 (int) |  ")
        print("                                                  | residue2       |  ")
        print("                                                  | ........       |  ")
        print("                                                  | residueM       |  ")
        print("                                                  ------------------  \n")

        print("  [Diameter hybrid region]     OPTIONAL           Value (in nm) of diameter of the hybrid region")
        print("                                                  Default value: 1.0 nm\n")
        print("  [NumberCpu]                  OPTIONAL           Integer number corresponding at the number of cpus that you would like to employ")
        print("                                                  for creating the CANVAS model of a biomolecule.")
        print("                                                  If '-n/--ncpu <nCPU>' is set, the code will be parallelized by employing nCPU cores")
        print("                                                  If '-n/--npu <nCPU>' is not set, the code will be parallelized by employing")
        print("                                                  the maximum number of allowed cores in your laptop/cluster\n")


    if(sys.argv[1] == "choice3"): 
        print("   List At-bb-CG FILE          MANDATORY          File organized in two columns")
        print("                                                  that contains the list of atomistic residues")
        print("                                                  and the list of residues that require an hybrid resolution:")
        print("                                                  -------------------------------------  ")
        print("                                         	 | res1_AT (int)   |  res1_hy (int)  |  ")
        print("                                         	 | res2_AT         |  res2_hy        |  ")
        print("                                          	 | res3_AT         |  ........       |  ")
        print("                                         	 | res4_AT         |  resM_hy        |  ")
        print("                                          	 | ........        |                 |  ")
        print("                                          	 | resN_AT         |                 |  ")
        print("                                          	 -------------------------------------  \n")

        print("  [NumberCpu]                  OPTIONAL           Integer number corresponding at the number of cpus that you would like to employ")
        print("                                                  for creating the CANVAS model of a biomolecule.")
        print("                                                  If '-n/--ncpu <nCPU>' is set, the code will be parallelized by employing nCPU cores")
        print("                                                  If '-n/--npu <nCPU>' is not set, the code will be parallelized by employing")
        print("                                                  the maximum number of allowed cores in your laptop/cluster\n")


    print("-----------------------------------------------------------------------------------------------------");
    print("Hereafter the list of flags:\n");

    print("   -g  --gro                   FILE               Coordinate FILE (gro format)")
    print("   -l  --list                  FILE               List AT-bb-CG FILE (txt or dat format)")

    if(sys.argv[1]=="choice1" or sys.argv[1]=="choice2"):
       print("  [-D] [--diameter]            INT/FLOAT          Value of diameter of hybrid region")

    print("  [-n] [--ncpu]                INT                Integer number corresponding at the number of cores employed for parallelizing the code")
    print("  [-h] [--help]                                   Give this help list\n")

    print("Report bugs to <raffaele.fiorentini@unitn.it>\n")
