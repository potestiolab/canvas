import os 

# Function that checks for the path where GROMACS is installed
def check_for_GRXRC(path_ff, forcefield):
    prefix_folder = ''.join(path_ff.split(forcefield))                            # Only way to make a difference between two strings.

    prefix_folder = ''.join(prefix_folder.split("/share/gromacs/top"))

    path_GMXRC = prefix_folder + "bin/GMXRC"

    Flag = False

    while(Flag == False):

        if(os.path.isfile(path_GMXRC)==True):
            print("\nThe file GMXRC was correctly found in {}... 28% completed\n".format(path_GMXRC))
            Flag = True

        else:
            print("\nI did not find the path where GMXRC is present")
            print("Thus, write directly the path of the executable file 'gmx' or 'gmx_mpi' is present")
            print("A possible path could be /usr/bin/gmx. Check it out")
            path_GMXRC = input("Please, type now the entire path of 'gmx' or 'gmx_mpi' (included 'gmx' or 'gmx_mpi'): \n")
            if(os.path.isfile(path_GMXRC)==True):
                print("\nThe executable file 'gmx' or 'gmx_mpi' was correctly found in {}... 28% completed\n".format(path_GMXRC))
                Flag = True
   
 
    return path_GMXRC 



# Function that checks if the gmx command is in the same folder of GMXRC or elsewhere. 
def check_for_gmx_command(path_ff, forcefield):         
    prefix_folder = ''.join(path_ff.split(forcefield))
  
    prefix_folder = ''.join(prefix_folder.split("/share/gromacs/top"))

    if(os.path.isfile(prefix_folder + "bin/gmx")==True): 
        gmx_command = "gmx"
        print("\n'gmx' was correctly found... 29% completed\n")
    
    elif(os.path.isfile(prefix_folder + "bin/gmx_mpi")==True):
        gmx_command = "gmx_mpi" 
        print("\n'gmx_mpi' was correctly found... 29% completed\n")
    
    else: 
        print("\nI did not find 'gmx' command for launching commands in GROMACS")
        print("'gmx' or 'gmx_mpi' is usually present in the same directory of GMXRC:")
        print("PREFIX_FOLDER/bin/\n")
        print("where:")
        print("PREFIX_FOLDER is the path where gromacs is installed")
        print("Check also if you have renamed the 'gmx' command or if it is not in the same directory of GMXRC\n")
        print("Please check it out and then relunch this python script\n\n")
        quit()

    return gmx_command 


# Function that checks if the forcefield chosen is present in gromacs working directory and it returns its path. 
# If not, the user is asked to write the path where the forcefield is present. 
def check_forcefield_dir(forcefield, path_ff):
    print("\nThe forcefield used is '{}'".format(forcefield))

    if(os.path.isdir(path_ff) == True):
        print("The force field used has been found in {}".format(path_ff))
        path_ff_new = path_ff

    else:
        print("I did not find the path where '{}' is present".format(forcefield))
        print("The path has usually the following shape:\n")
        print("PREFIX_FOLDER/share/gromacs/top/{}\n".format(forcefield))
        print("where:")
        print("PREFIX_FOLDER is the path where gromacs is installed\n")
        print("Check if the forcefield directory is in your current directory\n") 
        path_ff_new = input("Please, insert now the entire path: \n" )

    if(os.path.isdir(path_ff_new) == False):
        print("The path chosen is uncorrect. Please, repeat again\n")
        quit()

    else:
        print("The path chosen is correct... 7% completed\n")

    if(path_ff_new[-1] != "/"):
        path_ff_new = path_ff_new + "/"

    return path_ff_new


# Function that checks if "ffnonbonded.itp" files is present in the original path_ff (path_ff) or in the new one (path_ff_new) and it returns its path. 
def check_ffnonbonded(path_ff_new): 
    ffnonbondedfile = path_ff_new + "ffnonbonded.itp"
    return ffnonbondedfile

# Function that checks if "ffbonded.itp" files is present in the original path_ff (path_ff) or in the new one (path_ff_new) and it returns its path.
def check_ffbonded(path_ff_new): 
    ffbondedfile    = path_ff_new + "ffbonded.itp"
    return ffbondedfile

# Function that checks if "cmap.itp" files is present in the original path_ff (path_ff) or in the new one (path_ff_new) and it returns its path.
def check_cmap(path_ff_new): 
    cmapfile         = path_ff_new + "cmap.itp"
    return cmapfile

# Function that checks if "tip3p.itp" file is present in the original path_ff (path_ff) or in the new one (path_ff_new) and it returns its path.
def check_tip3p(path_ff_new):
    tip3pfile       = path_ff_new + "tip3p.itp"
    return tip3pfile 

# Function that checks if "ions.itp" file is present in the original path_ff (path_ff) or in the new one (path_ff_new) and it returns its path.
def check_ions(path_ff_new):
    ionsfile       = path_ff_new + "ions.itp"
    return ionsfile
