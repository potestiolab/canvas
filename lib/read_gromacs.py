# Function that reads the coordinate gromacs file (.gro). It returns: 
# - column[0] 		   --> the list of lists.
# - column[1] 		   --> the total number of atoms
# - column[2], [3] and [4] --> box sizes i.e. Lx, Ly and Lz. 
def readgro_all(f):        
                           
    lista = [] 

    i = 1
    title   = f.readline()

    i = 2
    n_atoms = f.readline()
    n_atoms = int(n_atoms)

    min_x = 10000
    min_y = 10000
    min_z = 10000

    old_res_number = 0 
    max_res_number = 0 

    old_at_number = 0 
    max_at_number = 0

    for i in range(n_atoms):
        line = f.readline()

        split = line.split() 

        if(len(split)==6 or len(split)==5):    # 5 or 6 beacuse when the index is higher than 9999, two column are merged in .gro file 

            res_number = int(line[0:5])
            res_name   = str(line[5:10]).strip()     
            at_name    = str(line[10:15]).strip()
            at_number  = int(line[15:20])
            pos_x      = float(line[20:28])
            pos_y      = float(line[28:36])
            pos_z      = float(line[36:44])

            ## This part is required because sometimes the res number starts again from 0 (in case 9999 is reached) 
            ## or from 1 (in case a chain is repeated). In this way the NEW_RES_NUMBER continues and does not start again from 1 (or 0).  
            if(res_number < old_res_number): 
                if(res_number == 0):
                    max_res_number = max_res_number + old_res_number + 1
                else:
                    max_res_number = max_res_number + old_res_number

            old_res_number = res_number
            new_res_number = res_number + max_res_number


            ## In analogy with the res_number, also the at_number could from 0 or 1 in case 99'999 is reached or 
            ## a chain is repeated. In the way the NEW_AT_NUMBER continues and does not start again from 1 (or 0). 
            if(at_number < old_at_number):
                if(at_number == 0):
                    max_at_number = max_at_number + old_at_number + 1
                else:
                    max_at_number = max_at_number + old_at_number

            old_at_number = at_number
            new_at_number = at_number + max_at_number



            lis=[new_res_number, res_name, at_name, new_at_number, pos_x, pos_y, pos_z]

            lista.append(lis)

            if(pos_x < min_x):
                min_x = pos_x 

            if(pos_y < min_y):
                min_y = pos_y 

            if(pos_z < min_z):
                min_z = pos_z 


        if(len(split)==9 or len(split)==8): 

            res_number = int(line[0:5])
            res_name   = str(line[5:10]).strip()     
            at_name    = str(line[10:15]).strip()
            at_number  = int(line[15:20])
            pos_x      = float(line[20:28])
            pos_y      = float(line[28:36])
            pos_z      = float(line[36:44])
            vel_x      = float(line[44:52])
            vel_y      = float(line[52:60])
            vel_z      = float(line[60:68])

            vel_x      = vel_x / 100
            vel_y      = vel_y / 100
            vel_z      = vel_z / 100

            ## This part is required because sometimes the res number starts again from 0 (in case 9999 is reached) 
            ## or from 1 (in case a chain is repeated). In this way the NEW_RES_NUMBER continues and does not start again from 1 (or 0).  
            if(res_number < old_res_number):
                if(res_number == 0):
                    max_res_number = max_res_number + old_res_number + 1
                else:
                    max_res_number = max_res_number + old_res_number

            old_res_number = res_number
            new_res_number = res_number + max_res_number


            ## In analogy with the res_number, also the at_number could from 0 or 1 in case 99'999 is reached or 
            ## a chain is repeated. In the way the NEW_AT_NUMBER continues and does not start again from 1 (or 0). 
            if(at_number < old_at_number):
                if(at_number == 0):
                    max_at_number = max_at_number + old_at_number + 1
                else:
                    max_at_number = max_at_number + old_at_number
            
            old_at_number = at_number
            new_at_number = at_number + max_at_number


            lis=[new_res_number, res_name, at_name, new_at_number, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z]

            lista.append(lis)

            if(pos_x < min_x):
                min_x = pos_x 

            if(pos_y < min_y):
                min_y = pos_y 

            if(pos_z < min_z):
                min_z = pos_z 

    min_x   = min_x * 10
    min_y   = min_y * 10
    min_z   = min_z * 10

    i=1

    box=f.readline()
    box=box.split()
    
    Lx = float(box[0])
    Ly = float(box[1])
    Lz = float(box[2])

    Lx = Lx * 10
    Ly = Ly * 10
    Lz = Lz * 10

    max_x = min_x + Lx
    max_y = min_y + Ly
    max_z = min_z + Lz

    return lista, n_atoms, Lx, Ly, Lz, title, min_x, min_y, min_z, max_x, max_y, max_z  


# Function that reads the topology file 
def readtop(ft):

    for line in ft:

        if 'Data prefix' in line:
             
            lin = line.replace(" ", "")    # remove empty spaces
            lin = lin.replace("\t", "")   # remove tabs 
            lin = lin.strip()                 
            break
  
        else:
            lin = ""

    ft.seek(0)
    for line in ft:
    
        if '; Include forcefield parameters' in line:
            break

    for line in ft: 
        lin2 = line.replace(" ", "")
        lin2 = lin2.replace("\t", "")
        lin2 = lin2.strip() 
        
        break 

    forcefield = lin2[9:-15]

    path_ff = lin[12:]  + "/share/gromacs/top/" + lin2[9:-15] 



    old_res_number = 0
    max_res_number = 0

    for line in ft: 
        line = line.strip() 
        if(line == "[ atoms ]"): 
            break 


    # Reading ATOMS: The program has stopped to "[ atoms ]" line, therefore we continue from this point reading all atoms
    atoms_topol = []

    for line in ft: 
        lin = line.strip() 
        if(lin == "[ bonds ]"): 
            break 

        if(line !="\n" and not line.isspace() and line[0] != ";" and line[0] != "#"):  # If the row is not empty, do not have whitespace (not isspace()),
                                                                                       #  the 1st letter of line is not ";" and "#" then we read the line 
            splt = line.split() 

            at_number  = int(splt[0])
            at_type    = str(splt[1])    # atomistic type; very important 
            res_number = int(splt[2])
            res_name   = str(splt[3])
            at_name    = str(splt[4])
            cgnr       = int(splt[5])
            charge     = float(splt[6])
            mass       = float(splt[7])

            ## This part is required because sometimes the res number starts again from 0 (in case 9999 is reached) 
            ## or from 1 (in case a chain is repeated). In this way the NEW_RES_NUMBER continues and does not start again from 1 (or 0).  
            if(res_number < old_res_number):
                if(res_number == 0):
                    max_res_number = max_res_number + old_res_number + 1
                else:
                    max_res_number = max_res_number + old_res_number

            old_res_number = res_number
            new_res_number = res_number + max_res_number


            top=[at_number, at_type, new_res_number, res_name, at_name, cgnr, charge, mass]
            atoms_topol.append(top)


    # Reading BONDS: The program has stopped to "[ bonds ]" line, therefore we continue from this point reading all bonds
    original_bonds_protein_fully_at = []

    for line in ft: 
        lin = line.strip()
        if(lin == "[ pairs ]"):
            break


        if(line !="\n" and not line.isspace() and line[0] != ";" and line[0] != "#"):
 
            splt = line.split() 

            first_atom   = int(splt[0])
            second_atom  = int(splt[1])
            func         = int(splt[2])

            bon = [first_atom, second_atom]
            original_bonds_protein_fully_at.append(bon)


    # Reading PAIRS: The program has stopped to "[ pairs ]" line, therefore we continue from this point reading all pairs 
    original_pairs_protein_fully_at = []

    for line in ft:
        lin = line.strip()
        if(lin == "[ angles ]"):
            break

        if(line !="\n" and not line.isspace() and line[0] != ";" and line[0] != "#"):

            splt = line.split() 

            first_atom  = int(splt[0])
            second_atom = int(splt[1])
            func        = int(splt[2])

            par = [first_atom, second_atom]
            original_pairs_protein_fully_at.append(par)


    # Reading ANGLES: The program has stopped to "[ angles ]" line, therefore we continue from this point reading all angles 
    original_angles_protein_fully_at = []

    for line in ft: 
        lin = line.strip()
        if(lin == "[ dihedrals ]"):
            break

        if(line !="\n" and not line.isspace() and line[0] != ";" and line[0] != "#"):
        
            splt = line.split() 

            first_atom   = int(splt[0])
            second_atom  = int(splt[1])
            third_atom   = int(splt[2])
            func         = int(splt[3])

            a = [first_atom, second_atom, third_atom, func]
            original_angles_protein_fully_at.append(a)

    # Reading DIHEDRALS: The program has stopped to "[ dihedrals ]" line, therefore we continue from this point reading all dihedrals (2,4,9,tors)
    original_dihedrals_all = []
    Flag_cmap = False

    for line in ft:
        lin = line.strip()
        if(lin == "[ cmap ]"):    # In case of Charmm.ff also "[ cmap ]" section is present 
            Flag_cmap = True
            break

        if(lin[:18] =="; Include Position"):
            break

        if(line != '\n' and line[0] != ';' and line[0] !="#" and line[0]!= "[" and not line.isspace()):
   
                splt = line.split() 

                first_atom   = int(splt[0])
                second_atom  = int(splt[1])
                third_atom   = int(splt[2])
                fouth_atom   = int(splt[3])
                func         = int(splt[4])

                if(len(splt) == 5): # dih2, dih4 and dih9 (no tors) 
                    a = [first_atom, second_atom, third_atom, fouth_atom, func]
                    original_dihedrals_all.append(a)
 
                if(len(splt) > 5):  # tors
                   string       = str(splt[5])
                   a = [first_atom, second_atom, third_atom, fouth_atom, func, string]
                   original_dihedrals_all.append(a)


    # Reading TORSIONAL CORRECTION MAP ([ cmap ] section): This section is present only in Charmm.ff. 
    #         There are 5 atom names that define the two torsional angles. Atoms 1-4 define "phi", while atoms 2-5 define "psi" 
    #         The corresponding atom types are then matched to the correct CMAP type in the cmap.itp file that contains the correction maps. 
    original_cmap_protein_fully_at = [] 

    if Flag_cmap == True:  # It means the "cmap" section exists 
    

        for line in ft:
            lin = line.strip() 
            if(lin[:18] =="; Include Position"):
                break 

            if(line != '\n' and line[0] != ';' and line[0] !="#" and not line.isspace()):

                splt = line.split()
                
                first_atom   = int(splt[0])
                second_atom  = int(splt[1])
                third_atom   = int(splt[2])
                fouth_atom   = int(splt[3])
                fifth_atom   = int(splt[4]) 
                func         = int(splt[5])

                c = [first_atom, second_atom, third_atom, fouth_atom, fifth_atom, func]
                original_cmap_protein_fully_at.append(c) 

    return atoms_topol, original_bonds_protein_fully_at, original_angles_protein_fully_at, original_dihedrals_all, original_pairs_protein_fully_at,\
           path_ff, forcefield, original_cmap_protein_fully_at


# Function that reads the ffnonbonded.itp file 
def read_ffnonbonded(f_nb):

    for line in f_nb:

        line =line.strip()

        if(line =="[ atomtypes ]"):
            break

    dict_sigma   = {}
    dict_epsilon = {}
    dict_mass    = {}
    dict_atnum   = {}
  
    pairtypes    = []

    for line in f_nb:
        lin = line.strip() 
        if(lin == "[ pairtypes ]"): 
            break 

        if(line !="\n" and not line.isspace() and line[0] != ";" and line[0] != "#"):    # If the row is not empty, do not have whitespace (not isspace()),
											 #  the 1st letter of line is not ";" and "#" then we read the line 
            splt = line.split()

            atomtypes   = str(splt[0])
            atom_number = int(splt[1]) 
            mass_attype = float(splt[2])
            sigma       = float(splt[5])       
            epsilon     = splt[6]   # In charmm27.ff we could find 0.0; --> we take only 0.0 (without semicolon ";") 
            if ";" in epsilon:
                epsilon = float(epsilon.split(";")[0])
            else:
                epsilon = float(epsilon)

            dict_sig    = {atomtypes : sigma}
            dict_eps    = {atomtypes : epsilon} 
            dict_m      = {atomtypes : mass_attype}
            dict_a      = {atomtypes : atom_number}

            dict_sigma.update(dict_sig)       # In case of same elements, e.g. OT in case of HEAVY_H or not, the second element overwrites the first one 
            dict_epsilon.update(dict_eps)     # This is good as, we do not consider heavy hydrogen (always first element), and we take only the second	   
            dict_mass.update(dict_m)          # i.e. not heavy hydrogen (in case of OT, the mass is 16 and not 7.9 for HEAVY_H) if using charmm forcefield  
            dict_atnum.update(dict_a) 


    # The program has stopped to "[ pairtypes ]" line, therefore we continue from this point reading all pairtypes.
    for line in f_nb: 
        if(line[0]!=";" and line[0]!="#" and line!="\n" and not line.isspace()):         # If the row is not empty, do not have whitespace (not isspace()),
											 # the 1st letter of line is not ";" and "#",  then we read the line 
            splt = line.split() 

            first_atom   = str(splt[0])
            second_atom  = str(splt[1])
            func         = int(splt[2])
            sigma14      = float(splt[3])                    # sigma14 for this couple of atoms 
            epsilon14    = float(splt[4])                    # epsilon14 for this couple of atoms 

            p_types =  [first_atom, second_atom, sigma14, epsilon14]
            pairtypes.append(p_types)
           
    return dict_sigma, dict_epsilon, dict_mass, dict_atnum, pairtypes   



# Function that reads the "ffbonded.itp file" 
def read_ffbonded(f_b):

    conv_factor_bonds = 1/836.8

    for line in f_b:
        line = line.strip()

        if(line =="[ bondtypes ]"):
            break

    bondtypes = [] 

    Flag_constr = False 

    for line in f_b:
        lin = line.strip() 
        if(lin == "[ constrainttypes ]"):    # Avoid reading constrainttypes in present (only in Amber.ff and not in Charmm)
            Flag_constr = True 
            break 
 
        if(lin =="[ angletypes ]"):
            break 

        if(line != '\n' and line[0] != ';' and line[0] !="#" and not line.isspace()): 
            splt = line.split() 

            first_atom   = str(splt[0])
            second_atom  = str(splt[1])
            func         = int(splt[2])       
            bond_len     = float(splt[3])                   # It is b0
            bond_const   = float(splt[4])                   # It is kb 

            bond_len     = bond_len * 10                    # Conversion nm to Angstrom
            bond_const   = bond_const * conv_factor_bonds   # Conversion kJ mol-1 nm-2 to kcal mol-1 A-2 (times 1/2 factor for lammps)

            b_types = [first_atom, second_atom, bond_len, bond_const]

            bondtypes.append(b_types)

    # If constrainttypes has been found ( Flag_constr == True), then the program stopped to "[ constrainttypes ]" line.
    # If constrainttypes has not been found (Flag_constr == False), then the program stopped to "[ angletypes ]" line. 

    if(Flag_constr == True):
        for line in f_b: 
            line = line.strip() 
            if(line == "[ angletypes ]"):
                break  

    # The program has stopped to "[ angletypes ]" line, therefore we continue from this point reading all angletypes.
    
    conv_fact_angle = 1/8.368
    angletypes = []


    for line in f_b:
        lin = line.strip() 

        if(lin == "[ dihedraltypes ]"):
            break 

        if(line != '\n' and line[0] != ';' and line[0] !="#" and not line.isspace()):
            splt = line.split() 

            first_atom   = str(splt[0])
            second_atom  = str(splt[1])
            third_atom   = str(splt[2])
            func         = int(splt[3])      
            theta        = float(splt[4])                # It is th0. This value is the same both for Gromacs and Lammps (deg) 
            theta_ener   = float(splt[5])                # It is cth

            theta_ener   = theta_ener * conv_fact_angle  # Conversion kJ mol-1 to kcal mol-1 (times 1/2 factor for Lammps)
             
            if(func == 5):            # There are two more parameters if func = 5. 
                ub0 = float(splt[6]) * 10                                   # Conversion nm to Angstrom
                kub = float(splt[7]) * conv_factor_bonds                    # Conversion kJ mol-1 nm-2 to kcal mol-1 A-2 (times 1/2 factor for lammps)

                a_types =  [first_atom, second_atom, third_atom, func, theta, theta_ener, ub0, kub] 

            elif(func == 1 or func == 2):  
                a_types = [first_atom, second_atom, third_atom, func, theta, theta_ener]

            angletypes.append(a_types)

    # The program has stopped to "[ dihedraltypes ]" line, therefore we continue from this point reading all dihedraltypes.
    # At most, it is possible to have dih4, dih9, torsion (part of dih9 only for amber99sb-ildn) and dih2. 
    # Usually:
    # 
    # o  *dih4* and *dih9* are found in every version of amber. 
    # o  *dih4*, *dih9* and *tors* are found ONLY in Amber99sb-ildn 
    # o  *dih9* and *dih2* are found in every version of Charmm forcefield.   
    conv_fact_dihedral  = 1/4.184
    conv_fact_dih2      = 1/8.368  # because in dih2 a factor 1/2 must be included. 
    dihedraltypes = []

    for line in f_b:
        if(line[0]!=";" and line!='\n' and not line.isspace() and line[0]!="[" and line[0]!= "#"):    # not line.isspace() ensures that 
												      # also whitespace are exluded
     
            splt = line.split()               # In case of dih2, multiplicity is not present. 

            first_atom    = str(splt[0])
            second_atom   = str(splt[1])
            third_atom    = str(splt[2])
            fourth_atom   = str(splt[3])
            func          = int(splt[4])
            phi_0         = float(splt[5])    # in ffbonded.itp is called phase   
            k_phi         = float(splt[6])    # in ffbonded.itp is called kd.  
            
            if(func == 2):     # i.e. function is 2 => dih2
                k_phi   = k_phi * conv_fact_dih2
                d_types = [first_atom, second_atom, third_atom, fourth_atom, func, phi_0, k_phi]
                dihedraltypes.append(d_types)

            elif(func == 4 or func == 9): 
                multiplicity  = int(splt[7])          # in ffbonded.itp is called pn: it is the multiplicity i.e. "n" as declared before. 
                k_phi   = k_phi * conv_fact_dihedral  # i.e. divided by 4.184: from kJ mol-1 to kcal mol-1 
                                                      # (we do not need the factore 1/2 this time) (for Lammps)
                d_types = [first_atom, second_atom, third_atom, fourth_atom, func, phi_0, k_phi, multiplicity]
                dihedraltypes.append(d_types)

        if(line[:7] =="#define"):  #torsion part: each interesting row starts with "#define"  --> Present only in Amber99sb-ildn
            splt = line.split()

            string_def    = str(splt[0])     # It stores the word "#define" 
            string_tors   = str(splt[1])
            phi_0         = float(splt[2])
            k_phi         = float(splt[3])
            multiplicity  = int(splt[4])

            k_phi = k_phi * conv_fact_dihedral
            d_types = [string_tors, phi_0, k_phi, multiplicity]
            dihedraltypes.append(d_types)
            
    #4th part# -> Let us split the dihedraltypes in dih2, dih4, dih9 and tors: 

    #             if column[4] of dihedraltypes is 4                => dih4 
    #             if column[4] of dihedraltypes is 9                => dih9 
    #             if column[4] of dihedraltypes is 2                => dih2 
    #             if column[0] of dihedraltypes start with TORSION  => tors

    dihedraltypes_dih4 = []
    dihedraltypes_dih9 = []
    dihedraltypes_dih2 = []
    dihedraltypes_tors = []


    for i in dihedraltypes:
      
        if(len(i) == 7): 
            if(i[4] == 2): 
                dihedraltypes_dih2.append(i) # In case of dih2, each internal list has 7 element: the 5th column i.e. i[4] == funct == 2. 
       
        if(len(i) == 8):       # In case of dih4 and dih9, each internal list has 8 elements: [at1, at2, at3, at4, funct, phi_0, k_phi, n]
            
            if(i[4] == 4):
                dihedraltypes_dih4.append(i)

            if(i[4] == 9):
                dihedraltypes_dih9.append(i)

        if(i[0][:7] == "torsion"):         # tors has internal list lenght equals to 4. 
            dihedraltypes_tors.append(i)


    return bondtypes, angletypes, dihedraltypes_dih4, dihedraltypes_dih9, dihedraltypes_tors, dihedraltypes_dih2


# Function that the cmap.itp file (only present in charmm.ff) 
def readcmap(f_cmap): 

    for line in f_cmap: 
        line = line.strip() 
        if(line == "[ cmaptypes ]"):
            break 

    cmaptypes = []   

    count = 0 

    for line in f_cmap: 
        lin = line.strip()
        if(lin[-6:-1] == "24 24"): 
            splt = line.split() 

            first_atom    = str(splt[0])
            second_atom   = str(splt[1])
            third_atom    = str(splt[2])
            fourth_atom   = str(splt[3])
            fifth_atom    = str(splt[4])  
            func          = int(splt[5]) 
            grid1         = int(splt[6])     # int, because the second-last term is 24 
            grid2         = str(splt[7])     # string, because the last term is 24\

            cm_types = [first_atom, second_atom, third_atom, fourth_atom, fifth_atom, func, grid1, grid2]
            cmaptypes.append(cm_types)

            for i in range(58): 
                string = f_cmap.readline()  
                cmaptypes[count].append(string) # the next 58 lines (57 lines with 10 elements and the last line with 6 elements) 
						# belong to cmaptypep just found (A,B,C,D,E,1,24,24\); thus we append the 58 string-line. 

            count = count + 1 

    return cmaptypes  


# Function that reads the topology for tip3p water 
def readtop_water(ft_w):

    for line in ft_w:

        line =line.strip()

        if(line =="[ atoms ]"):
            break

    atoms_topol_water = []

    for line in ft_w:
        if(line =="\n"):
                                        # If the row is empty, then we exit the loop beacuse we have read all the atoms in topology file...
            break

        if(line[0] !=";" and line[0] !="#"):             # If the 1st letter of line is ";" or "#" then skip the line itself, otherwise we split it. 
            splt = line.split()

            id_water   = int(splt[0])
            at_type    = str(splt[1])    # atomistic type; very important 
            res_number = int(splt[2])
            res_name   = str(splt[3])
            at_name    = str(splt[4])
            cgnr       = int(splt[5])
            charge     = float(splt[6])
            if(len(splt)>7): 
                mass       = float(splt[7])
            else:
                if(at_name == "OW"):
                    mass = float(16.00000) 
                if(at_name == "HW1" or at_name == "HW2"):
                    mass = float(1.00800)


            top=[id_water, at_type, res_number, res_name, at_name, cgnr, charge, mass]

            atoms_topol_water.append(top)


    #### Bonds 

    conv_factor_bonds = 1/836.8

    for line in ft_w:

        line = line.strip()

        if(line =="[ bonds ]"):
            break 

    original_bonds_water = []

    for line in ft_w:

        if not line.strip():             # if line is an empty row, even if it has whitespaces
            break 

        if(line[0]!=";"):       
            splt = line.split()
            first_atom   = int(splt[0])
            second_atom  = int(splt[1])
            func         = int(splt[2])     
            b0_water     = float(splt[3]) * 10 		        # Conversion nm to Angstrom
            kb_water     = float(splt[4]) * conv_factor_bonds   # Conversion kJ mol-1 nm-2 to kcal mol-1 A-2 (times 1/2 factor for lammps) 

    
            bon = [first_atom, second_atom, b0_water, kb_water]
            original_bonds_water.append(bon)


    #### Angles 

    conv_fact_angles = 1/8.368

    for line in ft_w:
        line = line.strip()

        if(line =="[ angles ]"):
            break

    original_angles_water = []

    for line in ft_w:

        if not line.strip():
            break
            
        if line[0]=="[":
            break

        if(line[0]!=";" and line[0]!="#"):
            splt = line.split()

            first_atom       = int(splt[0])
            second_atom      = int(splt[1])
            third_atom       = int(splt[2])
            func             = int(splt[3])
            theta_water      = float(splt[4])                        # theta value has the same units for gromacs and lammps      
            theta_ener_water = float(splt[5]) * conv_fact_angles     # Conversion kJ mol-1 nm-2 to kcal mol-1 A-2 (times 1/2 factor for lammps) 

            a = [first_atom, second_atom, third_atom, theta_water, theta_ener_water]
            original_angles_water.append(a)

    return atoms_topol_water, original_bonds_water, original_angles_water




# Function that reads the topology for IONS 
def readtop_ions(ft_ions):
    
    count = 0 
    for line in ft_ions:
        line = line.strip() 

        if(line == "[ atoms ]"): 
            count = count + 1           # In qeusto modo ci troviamo quanti volte "ATOMS" viene ripetuto, e sappamo quindi quanti ioni ci sono. 

    ft_ions.seek(0)

    atoms_topol_ions = []
    
    for i in range(count):

        for line in ft_ions:

            line =line.strip()

            if(line =="[ atoms ]"):
                break


        for line in ft_ions:
            if not line.strip():
                                        # If the row is empty, then we exit the loop beacuse we have read all the atoms in topology file...
                break

            if(line[0] !=";"):               # If the 1st letter of line is ";" then skip the line itself, otherwise we split it. 
                splt = line.split()

                id_ion   = int(splt[0])
                at_type    = str(splt[1])    # atomistic type; very important 
                res_number = int(splt[2])
                res_name   = str(splt[3])
                at_name    = str(splt[4])
                cgnr       = int(splt[5])
                charge     = float(splt[6])

                if(res_name == "SOD"): 
                    res_name = "NA" 
                elif(res_name == "CLA"):
                    res_name = "CL"

                if(at_name == "SOD"):
                    at_name = "NA" 
                elif(at_name == "CLA"):
                    at_name = "CL"

                top=[id_ion, at_type, res_number, res_name, at_name, cgnr, charge]
                atoms_topol_ions.append(top)
    
    return atoms_topol_ions 
