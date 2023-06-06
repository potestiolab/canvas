import os 
import sys

from multiprocessing import cpu_count

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



# Function: it checks if an argument of argparse (that will be treated always as string) is integer, float, or string. The function returns "x" 
def check_Int_Float_Str(x):
    if(isIntFloat(x)):
        try:
            x = int(x)
        except:
            x = float(x)
    return x



# Function: Checking if a mandatory file actually found or not. 
def checking_file_found_block(FileName):
    if not os.path.isfile(FileName):
        print("\n####################################################################################")
        print("ERROR. Error while opening the file. '{}' does not exist.\n".format(FileName))
        print("####################################################################################\n\n")
        print_usage_main_block()
        quit()


# Function: check if mandatory files are present 
def mandatory_files_present_CANVAS(grofile, listatomsfile, topfile):
    if (grofile is None):
        print("\n####################################################################################")
        print("ERROR. The Coordinate file is missing.")
        print("       Look below for further help.")
        print("####################################################################################\n\n")
        print_help_main_CANVAS()
        quit()

    if (listatomsfile is None):
        print("\n####################################################################################")
        print("ERROR. The file containing the list of survived atoms is missing.")
        print("       Look below for further help.")
        print("####################################################################################\n\n")
        print_help_main_CANVAS()
        quit()

    if (topfile is None): 
        print("\n####################################################################################")
        print("ERROR. The topology file is missing.")
        print("       Look below for further help.") 
        print("####################################################################################\n\n")
        print_help_main_CANVAS()
        quit() 


# Function: check if mandatory files are present in 'block.py' code 
def mandatory_files_present_block(grofile, listfile): 
    if (grofile is None):
        print("\n####################################################################################")
        print("ERROR. The Coordinate file is missing.")
        print("       Look below for further help.")
        print("####################################################################################\n\n")
        print_help_block()
        quit() 

    if (listfile is None):
        print("\n####################################################################################")
        print("ERROR. The List AT-bb-CG is missing.")
        print("       Look below for further help.")
        print("####################################################################################\n\n")
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

# Function: Printing help message if 'block.py' does not present valid arguments: 
#           in particular, the first argument must be 'choice1' or 'choice2', or 'choice3' according the user choice.
def checking_accepted_tasks_block():
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

        print("\n####################################################################################")
        print("ERROR. '{}' not in the list of accepted tasks!\n".format(sys.argv[1]))
        print("       Be sure that either 'bin' or 'density' task is the first argument. Other tasks are not allowed.")
        print("       Look below for more information.")
        print("####################################################################################\n\n")
        print_usage_main_block()
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



# Function: Parsing arguments, printing error and help message if 'Hs-Hk-plot.py' presents not allowed arguments, and checking if mandatory files are present 
def checking_valid_arguments_block(parser):
    try:        
        args  = parser.parse_args()
    except SystemExit:
        print("\n####################################################################################")
        print("ERROR. Arguments with no flag are not allowed. Check that each flag (e.g. -f) is followed by its specific argument")
        print("       Look below for more information.")
        print("####################################################################################\n\n")
        print_help_block()
        quit()


# Function: Checking if the optional argument "ncpu" is set. The program returns an error if the user askes for a number of cores 
#           higher than the maximum allowed. If 'ncpu' is not specified, the maximum number of cores is employed.  
def checking_errors_ncpu_opt_arg(ncpu):
    if(ncpu is not None):
        ncpu = check_Int_Float_Str(ncpu)  # it returns "ncpu" accordingly with its nature: integer, float, or string
        if(not isinstance(ncpu, int)): # if ncpu is not integer (i.e. float or string)
                        print("\n####################################################################################")
                        print("ERROR. '-n/--ncpu {}' set, but not recognized. Only an integer number higher than 0 is permitted.".format(ncpu))
                        print("       Please, insert an integer number higher than 0 (-n/--ncpu <INT>)")
                        print("       or ignore this flag, that means that the maximum number of cores are employed (-n/--ncpu {})".format(cpu_count()))
                        print("       Look below for further help.")
                        print("####################################################################################\n\n")
                        print_help_block()
                        quit()
        else:
            if(ncpu <= 0):
                print("\n####################################################################################")
                print("ERROR. '-n/--ncpu {}' set, but not recognized. Only an integer number higher than 0 is allowed.".format(ncpu))
                print("       A negative number or zero are meaningless.")
                print("       Please, insert an integer number (-n/--ncpu <INT>),")
                print("       or ignore this flag leaving that means that the maximum number of cores are employed (-n/--npu {})".format(cpu_count()))
                print("       Look below for further help.")
                print("####################################################################################\n\n")
                print_help_block()
                quit()

            if(ncpu > cpu_count()):
                print("\n####################################################################################")
                print("ERROR. '-n/--ncpu {} set, but not recognized. The max number of cores available in your cluster/laptop is {}.".format(ncpu, cpu_count()))
                print("       You are asking for {} cores. Please, take care of it and repeat!".format(ncpu))
                print("       Look below for further help.")
                print("####################################################################################\n\n")
                print_help_block()
                quit()

        ncpu_employed = ncpu
        print("\no '-n/--ncpu {}' set. The number of cores employed is {}... 40% completed. \n".format(ncpu_employed, ncpu_employed))

    if(ncpu is None):
        ncpu_employed = cpu_count()
        print("\no '-n/--ncpu {}' set. The number of cores employed is {}... 40% completed. \n".format(ncpu_employed, ncpu_employed))
 
    return ncpu_employed



# Function: Checking if the optional argument "diameter" is set.  
#           As first, this argument can be set ONLY if the 'choice1' or 'choice2' option is set, otherwise an error occurs. 
#           If 'choice1' or 'choice2' option is set, moreover, an error occurs if a string or negative (incl. zero) is inserted. 
#           Only an integer or float number is accepted. 
#           Consider that, after defining arguments, the latter are treated always as STRING, thus we have to check if we are deal with INT, FLOAT or STRINGS.    

def checking_errors_diameter_opt_arg(diameter): 
    if(diameter is not None):   # If diameter exists... 
        if(sys.argv[1]=='choice1' or sys.argv[1]=='choice2'):      # If 'choice1' or 'choice2' option is set  
            diameter = check_Int_Float_Str(diameter)               # it returns "diameter" accordingly with its nature: integer, float, or string
 
            if(isinstance(diameter, str)):                         # if diameter is string 
                print("\n####################################################################################")
                print("ERROR. '-D/--diameter {}' set, but not recognized. Only an integer or float number higher than 0 is permitted.".format(diameter))
                print("       Please, insert an integer or float number higher than 0 (-D/--diameter <INT/FLOAT>)")
                print("       or ignore this flag leaving the default value of 1.0 (-D/--diameter 1.0)")
                print("       Look below for further help.")
                print("####################################################################################\n\n")
                print_help_block()
                quit()

            else:
                if(diameter <= 0):
                    print("\n####################################################################################")
                    print("ERROR. '-D/--diameter {}' set, but not recognized. Only an integer or float number higher than 0 is allowed.".format(diameter))
                    print("       A negative number or zero are meaningless.")
                    print("       Please, insert an integer or float number (-D/--diameter <INT/FLOAT>)")
                    print("       or ignore this flag leaving the default value of 1.0 (-D/--diameter 1.0).")
                    print("       Look below for further help.")
                    print("####################################################################################\n\n")
                    print_help_block()
                    quit()

            print("\no '-D/--diameter {}' set. The diameter 'D' of medium-grained region: {} nm... 30% completed. \n".format(diameter, diameter))


        if(sys.argv[1]=='choice3'):  # If 'choice3' option is set
            print("\n####################################################################################")
            print("ERROR. With 'choice3' option the optional argument 'diameter' (-D/--diameter <INT/FLOAT>) cannot be set.")
            print("       Please, do not use this flag if 'choice3' option is set.")
            print("       Look below for further help.")
            print("####################################################################################\n\n")
            print_help_block()
            quit()


    if(diameter is None):
        if(sys.argv[1]=='choice1' or sys.argv[1]=='choice2'):
            diameter = 1.0
            print("\no '-D/--diameter' not set. The diameter 'D' of medium-grained region is {} nm. Diameter = 1.0 nm (default value)... 30% completed. \n".format(diameter))

    return diameter



# Function: Checking if the mandatory argument "listfile" is correctly set for *choice1*.
#           This file has to be organized in two columns. The first one must be a vertical list of residue numbers (INT NUMBERS) (not atom_numbers) 
#           from which the atomistic sphere is traced. The radius of this sphere is defined in the second column: it can be 
#           integer or float, and it is the radius value in nm.  

def checking_errors_listfile_choice1_mand_arg(listfile):  
    f = open(listfile, "r")

    list_AT_res      = []
    list_radii       = []

    count_line = 1
    for line in f:
        check_empty_rows_block(line, listfile)
        line = line.split()

        if(len(line) != 2):       #N_lines = 2 
            print("\n####################################################################################")
            print("ERROR. *choice1* selected: the input file '{}' cannot be read correctly.".format(listfile))
            print("       Each line must contain 2 columns: the {}th row has {} columns".format(count_line, len(line)))
            print("       Be sure that '{}' is organized in two different columns:".format(listfile))
            print("       • 1st column       INT                List of residues numbers from which the atomistic sphere will be traced.")
            print("       • 2nd column       INT/FLOAT          Atomistic radius of each residue in nm.")
            print("       Look below for further information.")
            print("####################################################################################\n\n")
            print_help_block()
            quit()

        else:
            if not line[0].isdigit():  #1st column = INT
                print("\n####################################################################################")
                print("ERROR. *choice1* selected: the input file '{}' cannot be read correctly.".format(listfile))
                print("        The 1st column must be an integer number. '{}' at {}th row is not integer.". format(line[0], count_line))
                print("        Be sure that '{}' is organized in two different columns:".format(listfile))
                print("        • 1st column       INT                List of residues numbers from which the atomistic sphere will be traced.")
                print("        • 2nd column       INT/FLOAT          Atomistic radius of each residue in nm.")
                print("        Look below for further information.")
                print("####################################################################################\n\n")
                print_help_block()
                quit()

            if not isIntFloat(line[1]): #2nd column = INT/FLOAT
                print("\n####################################################################################")
                print("ERROR. *choice1* selected: the input file '{}' cannot be read correctly.".format(listfile))
                print("       The 2nd column must be a Float or Integer number. '{}' at {}th row is a string.". format(line[1], count_line))
                print("       Be sure that '{}' is organized in two different columns:".format(listfile))
                print("       • 1st column       INT                List of residues numbers from which the atomistic sphere will be traced.")
                print("       • 2nd column       INT/FLOAT          Atomistic radius of each residue in nm.")
                print("       Look below for further information.")
                print("####################################################################################\n\n")
                print_help_block()
                quit()

        res    = int(line[0])
        radius = float(line[1])

        list_AT_res.append(res)
        list_radii.append(radius)

        count_line = count_line + 1

    return list_AT_res, list_radii


# Function: Checking if the mandatory argument "listfile" is correctly set for *choice2*.
#           This file has to be organized in one column containing a vertical list of all residue numbers (INT NUMBERS) (not atom_numbers) 
#           modelled atomistically.  
def checking_errors_listfile_choice2_mand_arg(listfile):
    f = open(listfile, "r")

    list_AT_residues = []

    count_line = 1
    for line in f:
        check_empty_rows_block(line, listfile)
        line = line.split()

        if(len(line) != 1):
            print("\n####################################################################################")
            print("ERROR. *Choice2* selected: the input file '{}' cannot be read".format(listfile))
            print("       Each line must contain only 1 column: the {}th row has {} columns".format(count_line, len(line)))
            print("       Be sure that '{}' is organized in only one column containing the list of residue numbers modelled atomistically.".format(listfile))
            print("       Look below for further information.")
            print("####################################################################################\n\n")
            print_help_block()
            quit()

        else:
            if not line[0].isdigit():
                print("\n####################################################################################")
                print("ERROR. *choice2* selected: the input file '{}' cannot be read correctly.".format(listfile))
                print("       The unique column must contain an integer number. '{}' at {}th row is not integer.". format(line[0], count_line))
                print("       Be sure that '{}' is organized in only one column containing the list of residue numbers modelled atomistically.".format(listfile))
                print("       Look below for further information.")
                print("####################################################################################\n\n")
                print_help_block()
                quit()


        res = int(line[0])
        list_AT_residues.append(res)

        count_line = count_line + 1

    return list_AT_residues



# Function: Checking if the mandatory argument "listfile" is correctly set for *choice3*.
#           This file has to be organized in two columns: 
#           o 1st: It contains a vertical list of all residue numbers (INT NUMBERS) (not atom_numbers) that require a fully-atomistic description.
#           o 2nd: It contains a vertical list of all residue numbers (INT NUMBERS) (not atom_numbers) that require a medium-grained description.   
def checking_errors_listfile_choice3_mand_arg(listfile):

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
            print("\n####################################################################################")
            print("ERROR. *Choice3* selected: the input file '{}' cannot be read".format(listfile))
            print("       Each line must contain 1 or 2 columns: the {}th row has {} columns:".format(count_line, len(line)))
            print("       Be sure that each row of '{}' is organized in only one or two columns.")
            print("       Look below for further information.")
            print("####################################################################################\n\n")
            print_help_block()
            quit()

        else:
            if(len(line) == 1):
                if not line[0].isdigit():
                    print("\n####################################################################################")
                    print("ERROR. *choice3* selected: the input file '{}' cannot be read correctly.".format(listfile))
                    print("       The 1st column must be an integer number. '{}' at {}th row is not an integre number.". format(line[0], count_line))
                    print("       Be sure that each row of '{}' is organized in one or two columns:".format(listfile))
                    print("       • 1st column       INT                List of residue numbers modelled atomistically")
                    print("       • 2nd column       INT                List of residue numbers that require an hybrid treatment")
                    print("       Look below for further information.")
                    print("####################################################################################\n\n")
                    print_help_block()
                    quit()
                else:
                    AT_res  = int(line[0])
                    list_AT_residues.append(AT_res)
                    count_line = count_line + 1

            if(len(line) == 2):
                if not line[1].isdigit():
                    print("\n####################################################################################")
                    print("ERROR. *choice3* selected: the input file '{}' cannot be read correctly.".format(listfile))
                    print("       The 2nd column must be an integer number. '{}' at {}th row is not an integer number.". format(line[1], count_line))
                    print("       Be sure that each row of '{}' is organized in one or two columns:".format(listfile))
                    print("       • 1st column       INT                List of residue numbers modelled atomistically")
                    print("       • 2nd column       INT                List of residue numbers that require an hybrid treatment")
                    print("       Look below for further information.")
                    print("####################################################################################\n\n")
                    print_help_block()
                    quit()
                else:
                    AT_res  = int(line[0])
                    bb_res  = int(line[1])
                    list_AT_residues.append(AT_res)
                    list_bb_residues.append(bb_res)
                    count_line = count_line + 1

    return list_AT_residues, list_bb_residues 
