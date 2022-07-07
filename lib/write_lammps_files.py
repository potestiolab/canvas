import itertools
from operator import itemgetter

# Function that writes the .lmp file of Lammps 
def write_lmp_file(f_lmp, number_of_atoms, number_of_bonds, number_of_angles, number_of_dihedrals, number_of_impropers, number_atom_types, \
                   number_bond_types, number_angle_types, number_dihedral_types, number_improper_types, min_x, max_x, min_y, max_y, \
                   min_z, max_z, list_survived_atoms, bonds_mult_res_prot, angles_mult_res_prot, dihedrals_mult_res_prot, atoms_ions, \
                   type_OW, type_HW): 

    ## A- Writing title, number of: atoms, bonds, angles, dihedrals, and impropers 
    f_lmp.write("File ilmp Variable resolution Protein in fully atomistic water\n\n")
    
    f_lmp.write("{} atoms\n".format(number_of_atoms))
    f_lmp.write("{} bonds\n".format(number_of_bonds))
    f_lmp.write("{} angles\n".format(number_of_angles))
    f_lmp.write("{} dihedrals\n".format(number_of_dihedrals))
    f_lmp.write("{} impropers\n\n".format(number_of_impropers))


    ## B- Writing the number of atom-, bond-, angle-, and dihedral-types 
    f_lmp.write("{} atom types\n".format(number_atom_types))
    f_lmp.write("{} bond types\n".format(number_bond_types))
    f_lmp.write("{} angle types\n".format(number_angle_types))
    f_lmp.write("{} dihedral types \n".format(number_dihedral_types))
    f_lmp.write("{} improper types \n\n".format(number_improper_types))


    ## C- Writing Min and Max coordinates of box
    f_lmp.write("{:f}{:12.7f} xlo xhi\n".format(min_x, max_x))
    f_lmp.write("{:f}{:12.7f} ylo yhi\n".format(min_y, max_y))
    f_lmp.write("{:f}{:12.7f} zlo zhi\n\n".format(min_z, max_z))
    
    
    ## D- Writing Masses of each atomtype
    f_lmp.write("Masses\n\n")
    
    for count in range(1, number_atom_types + 1):
        for i in list_survived_atoms:
            if(count == i[-1]):
                f_lmp.write("{:<5d}{:>3.7f}  #{}\n".format(i[-1], i[10], i[16]))
                break
    
    
    ## E- Writing Bond Coefficients for each atomtype present. Recall that the water (HW-OW) bond has a unique bondtype 
    ##    This section is divided in four columns: 
    ##    i[0] --> $count where $count goes from 1 to N, and N is the maximum value of bondtype 
    ##    i[1] --> bond style, namely "harmonic"
    ##    i[2] --> Kb (energy) (Lammps Units)  
    ##    i[3] --> b0 (eql.distance) (Lammps Units, i.e. Angstrom) 
    
    f_lmp.write("\nBond Coeffs\n\n")
    
    for count in range(1, number_bond_types + 1):
        for i in bonds_mult_res_prot:
            if(count == i[-1]):                                                                       # namely the specific number of bondtype
                f_lmp.write("{} harmonic {:17.8e} {:17.8e}\n".format(i[-1], i[6]/836.8, i[7]*10))     # i[6]=Ek, i[7]=b0 both in Lammps Units. 
                break
    
    
    ## F- Writing Angle Coefficients for each atomtype: Recall that the water (HW-OW-HW) angle has a unique angletype 
    ##    If the functional form of angles is 1, this section is divided in four coulumns: 
    ##
    ##    i[0] --> $count where $count goes from 1 to N, and N is the maximum value of angletype 
    ##    i[1] --> angle style, namely "harmonic"
    ##    i[2] --> cth (in Lammps units). (angle energy)
    ##    i[3] --> th0 (eql. angle, in deg.) 
    ##
    ##    If the fucntional form of angles is 5, this section is divided in six columns: 
    ##
    ##     i[0] --> $count where $count goes from 1 to N, and N is the maximum value of angletype
    ##     i[1] --> angle style, namely "charmm"
    ##     i[2] --> k_phi (in Lammps units). (angle energy)
    ##     i[3] --> phi0 (eql. angle, in deg.) 
    ##     i[4] --> kub (in Lammps Units) (Urey Bradley component of angle energy divided by distance squared) [energy/distance^2]
    ##     i[5] --> ub (in Lammps Units) (Urey Bradley component of distance between atom 1 and atom 3) 
    
    f_lmp.write("\nAngle Coeffs\n\n")
    
    for count in range(1, number_angle_types + 1):
        for i in angles_mult_res_prot:    
            if(i[3] == 1):               							               # If funct == 1 
                if(count == i[-1]):     							               # namely the specific number of angletype 
                    f_lmp.write("{} harmonic {:17.8e} {:17.8e}\n".format(i[-1], i[7]/8.368, i[8]))  	       # Lammps Units 
                    break
    
            if(i[3] == type_OW or i[3] == type_HW):   							       # if we are dealing with water atoms: 
                if(count == i[-1]): 
                    f_lmp.write("{} harmonic {:17.8e} {:17.8e}\n".format(i[-1], i[6]/8.368, i[7]))
                    break
    
            if(i[3] == 5):  										       # if funct == 5
                if(count == i[-1]):     								       # namely the specific number of angletype
                    f_lmp.write("{} charmm {:17.8e} {:17.8e} {:17.8e} {:17.8e}\n".format(i[-1], i[7]/8.368, i[8], i[10]/836.8, i[9]*10))  # Lammps Units 
                    break
    
    
    ## G- Writing Dihedral Coefficients. According the multiplicity and the torsional angle value, the style 
    ##    could be multi/harmonic or charmm (only for dih4, dih9, and tors): 
    ##
    ##    $i  multi/harmonic  $1st_coeff  $2nd_coeff   $3rd_coeff   $4th coeff   $5th coeff
    ##
    ##    $i  charmm          $k_phi      $n           $phi0        0.0  
    
    if(number_dihedral_types != 0):
        f_lmp.write("\nDihedral Coeffs\n\n")
    
    for count in range(1, number_dihedral_types + 1):
        for i in dihedrals_mult_res_prot:
    
            if(count == i[-1]):      					                              # namaly the specific number of dihedraltype 
                if(i[-2] == "multi"):
                    f_lmp.write("{} multi/harmonic {:17.8e} {:17.8e} {:17.8e} {:17.8e} {:17.8e}\n".format(i[-1], i[-7], i[-6], i[-5], i[-4], i[-3]))
    
                if(i[-2] == "charmm"):
                    f_lmp.write("{} charmm {:17.8e} {:5d} {:4d} {:5d}\n".format(i[-1],i[-6],i[-5],int(i[-4]),int(i[-3])))    # i[7] has to be integer.  
    
    ## H- Writing the Improper Coefficients. If charmm.ff is employed, dihedrals with functional form equals to 2 is also used 
    ##    that Lammps traduces in impropers.  
    ##
    ##    E = k(X-X0)^2    X = chi greek symbol
    ##
    ##    This section is organized in four columns: 
    ##
    ##    i[0] --> $count where $count goes from 1 to N, and N is the maximum value of impropertype 
    ##    i[1] --> improper style, namely "harmonic"
    ##    i[2] --> Kb (energy prefactor) (Lammps Units)  
    ##    i[3] --> Chi0 (X0) (eql.value of improper angle) (in deg) 
    
    if(number_improper_types != 0):
        f_lmp.write("\nImproper Coeffs\n\n")
    
    for count in range(1, number_improper_types + 1):
        for i in dihedrals_mult_res_prot:
            if(count == i[-1]):      								       # namely the specific value of impropertype: 
                if(i[4] == 2):       								       # i[4]==2 means funct==2
                    f_lmp.write("{} harmonic {:17.8e} {:17.8e}\n".format(i[-1], i[10]/8.368, i[9]*10))     # i[10] = K, i[9] = X0 both in Lammps Units. 
                    break
    
    ## I- ATOMS: Writing everything about the atoms. This section is organized in 7 columns: 
    ##
    ##           i[0] --> $i, where i = [1, N_atoms]
    ##            i[1] --> molecule-tag, i.e. the new-res-name (i[17] of list_survived_atoms)
    ##            i[2] --> atomtype number (i[-1] of list_survived_atoms)
    ##            i[3] --> Charge Q of each survived atom (i[11] of list_survived_atoms)
    ##            i[4] --> NEW x-position (in Angstrom) (i[4] of list_survived_atoms) 
    ##	          i[5] --> NEW y-position (in Angstrom) (i[5] of list_survived_atoms) 
    ## 	          i[6] --> NEW z-position (in Angstrom) (i[6] of list_survived_atoms) 
    
    f_lmp.write("\nAtoms\n\n")
    
    count = 1
    for i in list_survived_atoms:
        f_lmp.write("{:6d}{:7d}{:7d} {:13.8f}  {:15.8f} {:14.8f} {:14.8f}\n".format(count, i[17], i[-1], i[11], i[4]*10, i[5]*10, i[6]*10))
        count = count + 1
    
    
    ## J- VELOCITIES: Writing the initial velocity for all atoms. In general, since there is no equilibration, they should be setted at 0.000.
    ##                In general, in case the solvated file has velocities, take in account that conversion is done:
    ##
    ##                [nm/psec] --> [Angstrom/fsec] 
    ##          
    ##                Therefore the factor is 100 beacuse "1nm = 10 Angstrom" and "1ps = 1000 fsec"
     
    f_lmp.write("\nVelocities\n\n")
    
    n_columns = len(atoms_ions[1])              # Take the first element e compute the lenght. If lenght = 7 => NO VELOCITIES 
                                                # Se lenght = 10, allora YES VELOCITIES 
    count = 1
    
    for i in atoms_ions:
        if(n_columns !=7):
            f_lmp.write("{:6d} {:14.8f} {:14.8f} {:14.8f}\n".format(count, i[7], i[8], i[9]))
            count = count + 1
    
        else:
            vx = 0.000
            vy = 0.000
            vz = 0.000
    
            f_lmp.write("{:6d} {:14.8f} {:14.8f} {:14.8f}\n".format(count, vx, vy, vz))
            count = count + 1
    
    ## K- BONDS: Writing the bonds (in terms of indexes) like to [ bonds ] section of Gromacs. 
    ##           This section is organized in 4 columns:
    ##
    ##           i[0] --> $i where i goes between 1 and the total number of bonds.
    ##           i[1] --> bondtype of reference (j[8] of bonds_mult_res_prot) 
    ##	         i[2] --> bond index of atom 1 (j[0] of bonds_mult_res_prot) 
    ##	         i[3] --> bond index of atom 2 (j[1] of bonds_mult_res_prot) 
    
    f_lmp.write("\nBonds\n\n")
    
    #count = 1
    #for i in bonds_mult_res_prot:
    #    if(count <= 99999):
    #        f_lmp.write("{:6d} {:7d}   {} {}\n".format(count, i[8], i[0], i[1]))
    #        count = count + 1
    # 
    #    else:
    #        f_lmp.write("{:6d} {:7d}   {} {}\n".format(count, i[8], i[0]+count, i[1]+count))   # This condition because otherwise count begins again from 0
    #        count = count + 1                                                                  

    count = 1 
    for i in bonds_mult_res_prot:
        f_lmp.write("{:6d} {:7d}   {} {}\n".format(count, i[8], i[0], i[1]))
        count = count + 1    
    

    ## L- ANGLES: Writing the Angles (in terms of indexes) like to [ angles ] section of Gromacs. 
    ##            This section is organized in 5 columns: 
    ##  
    ##            i[0] --> $i where i goes between 1 and the total number of angles. 
    ##            i[1] --> angletype of reference (j[11] or j[-1] of bonds_mult_res_prot) 
    ##            i[2] --> angle index of atom 1 (j[0] of angles_mult_res_prot) 
    ##            i[3] --> angle index of atom 2 (j[1] of angles_mult_res_prot)
    ##            i[4] --> angle index of atom 3 (j[2] of angles_mult_res_prot) 
             
    f_lmp.write("\nAngles\n\n")
    
    #count = 1
    #for i in angles_mult_res_prot:
    #    if(count <= 99999):
    #        f_lmp.write("{:6d} {:7d}   {} {} {}\n".format(count, i[-1], i[0], i[1], i[2]))
    #        count = count +1
    # 
    #    else:
    #        f_lmp.write("{:6d} {:7d}   {} {} {}\n".format(count, i[-1], i[0]+count, i[1]+count, i[2]+count))
    #        count = count +1

    count = 1
    for i in angles_mult_res_prot:
        f_lmp.write("{:6d} {:7d}   {} {} {}\n".format(count, i[-1], i[0], i[1], i[2]))
        count = count +1
    
    
    ## M- DIHEDRALS: Writing the Dihedrals (in terms of indexes) like to [ dihedrals ] section of Gromacs.
    ##               This section is organized in 6 columns: 
    ## 
    ##               i[0] --> $i where i goes between 1 and the total number of dihedrals. 
    ##               i[1] --> dihedraltype of reference (j[-1] of dihedrals_mult_res_prot with funct=4 or 9) 
    ##               i[2] --> dihedral index of atom 1 (j[0] of dihedrals_mult_res_prot with funct=4 or 9) 
    ##               i[3] --> dihedral index of atom 2 (j[1] of dihedrals_mult_res_prot with funct=4 or 9)
    ##               i[4] --> dihedral index of atom 3 (j[1] of dihedrals_mult_res_prot with funct=4 or 9) 
    ##               i[5] --> dihedral index of atom 4 (j[2] of dihedrals_mult_res_prot with funct=4 or 9) 
    
    if(number_dihedral_types != 0):
        f_lmp.write("\nDihedrals\n\n")
    
        #count = 1
        # 
        #for i in dihedrals_mult_res_prot:
        #    if(i[4]==9 or i[4]==4):
        #        if(count <= 99999):
        #            f_lmp.write("{:6d} {:7d}   {} {} {} {}\n".format(count, i[-1], i[0], i[1], i[2], i[3]))
        #            count = count + 1
        #        else:
        #            f_lmp.write("{:6d} {:7d}   {} {} {} {}\n".format(count, i[-1], i[0]+count, i[1]+count, i[2]+count, i[3]+count))
        #            count = count + 1

        count = 1 
        for i in dihedrals_mult_res_prot:
            if(i[4]==9 or i[4]==4):
                f_lmp.write("{:6d} {:7d}   {} {} {} {}\n".format(count, i[-1], i[0], i[1], i[2], i[3]))
                count = count + 1

    
    ## N- IMPROPERS: Writing the Impropers (in terms of indexes) like to [ dihedrals ] section of Gromacs. 
    ##               Recall that impropers in Lammps correspond to dihedrals with funct=2 in Gromacs. 
    ##               This section is organized in 6 columns: 
    ##
    ##               i[0] --> $i where i goes between 1 and the total number of impropers. 
    ##               i[1] --> dihedraltype of reference (j[-1] of dihedrals_mult_res_prot with funct=2) 
    ##               i[2] --> dihedral index of atom 1 (j[0] of dihedrals_mult_res_prot with funct=2) 
    ##               i[3] --> dihedral index of atom 2 (j[1] of dihedrals_mult_res_prot with funct=2)
    ##               i[4] --> dihedral index of atom 3 (j[1] of dihedrals_mult_res_prot with funct=2) 
    ##               i[5] --> dihedral index of atom 4 (j[2] of dihedrals_mult_res_prot with funct=2) 
    
    if(number_improper_types != 0):
        f_lmp.write("\nImpropers\n\n")
    
        #count = 1    
        #for i in dihedrals_mult_res_prot:
        #    if(i[4]==2):
        #        if(count <= 99999):
        #            f_lmp.write("{:6d} {:7d}   {} {} {} {}\n".format(count, i[-1], i[0], i[1], i[2], i[3]))
        #            count = count + 1
        #        else:
        #            f_lmp.write("{:6d} {:7d}   {} {} {} {}\n".format(count, i[-1], i[0]+count, i[1]+count, i[2]+count, i[3]+count))
        #            count = count + 1


        for i in dihedrals_mult_res_prot:
            if(i[4]==2):
                f_lmp.write("{:6d} {:7d}   {} {} {} {}\n".format(count, i[-1], i[0], i[1], i[2], i[3]))
                count = count + 1



def write_input_lammps_file(f_input, dihedrals_mult_res_prot, Flag_cmaps, number_of_bonds, number_of_angles, number_of_dihedrals, \
                            number_of_impropers, number_atom_types, list_survived_atoms, pairs_mult_res_prot, dict_at_types, \
                            bonds_mult_res_prot, angles_mult_res_prot, type_OW, type_HW, solvate_Flag):

    ## a- Checking if charmm forcefield is employed
    Flag_charmm = False
    
    for i in dihedrals_mult_res_prot:
        if(i[4] == 2):
            Flag_charmm = True
            break
    
    if(Flag_charmm == False):     # If dih2 are not found, let us give a look to cmap 
        if(Flag_cmaps == True):   # cmap are found, i.e. the charmm.ff is employed
            Flag_charmm = True
    
    
    ## b- Writing parameters...
    
    
    f_input.write("units real\n")          # Real Units  
    f_input.write("atom_style full\n\n")   # No HadResS
    f_input.write("dimension 3\n")         # 3 dimensions
    f_input.write("boundary p p p\n\n")    # Boundary conditions
    
    if(number_of_bonds!=0):
        f_input.write("bond_style       hybrid harmonic\n")                        # harmonic style for bonds

    if(number_of_angles!=0):
        if(Flag_charmm == False):    # in case of amber.ff 
            f_input.write("angle_style      hybrid harmonic\n")                    #   If amber.ff is employed only harmonic style for angles is possible 
        if(Flag_charmm == True and solvate_Flag == True):                          #   If charmm.ff is employed and the system is solvated the angle    
            f_input.write("angle_style      hybrid harmonic charmm\n")             # of water is harmonic, the remainder is charmm 
        if(Flag_charmm == True and solvate_Flag == False):                         #   If charmm.ff is employed and the system is NOT solvated only 
            f_input.write("angle_style      hybrid charmm\n")                      # charmm style is possible (as the water harmonic style is not present)

    if(number_of_dihedrals!=0):
        f_input.write("dihedral_style   hybrid charmm multi/harmonic\n")           # multi/harmonic (5 coeffs) and charmm (more than 5 coeffs) style.  
    
    if(number_of_impropers!=0):
        f_input.write("improper_style   hybrid harmonic\n")                        # harmonic style for impropers (dihedrals with funct = 2) 
    
    if(Flag_charmm == True):
        f_input.write("special_bonds    lj 0.0 0.0 1.0 coul 0.0 0.0 1.0\n\n")      # For charmm.ff no general scaling factor is present as for amber.ff 
    									           # Each pairtype will be read in [pairtypes] sections. 
    else: 
        f_input.write("special_bonds    lj 0.0 0.0 0.5 coul 0.0 0.0 0.8333\n\n")   # Special bonds for LJ and Coulomb 
                                                                                   # First neighbours  : No Coul & LJ  
                                                                                   # Second neighbours : No Coul & LJ 
                                                                                   # Third neighbours  : Rescaled Coul and LJ (0.8333 and 0.5) 
                                                                                   #                     [0.5 means factor of rescaling equals to 1/2]
    
    
    f_input.write("read_data        Multiple-Res-Protein-in-Water.lmp\n")          # Reading .lmp file  
    f_input.write("#read_restart    restart.Mult-Res-T300.50000\n\n")              # It is required in case of restart of simulation. 
    f_input.write("######################################################################\n\n")
    
    
    ## c- Writing variables... 
    
    f_input.write("variable         root index Mult-Res-T300\n")
    f_input.write("variable         Nrun equal 500000   #500 ps\n")      # 500'000 fsec i.e. 500 ps, i.e. 0.5 ns
    f_input.write("variable         Nf equal 100\n")
    f_input.write("variable         Ne equal 100\n")
    f_input.write("variable         Nr equal ${Nf}/${Ne}\n")
    f_input.write("variable         Ndump equal 10000\n")
    f_input.write("variable         Nrestart equal 10000\n")
    f_input.write("variable         Text equal 300.0\n")                 # External Temperature
    f_input.write("variable         Pext equal 10.0\n\n")                # External Pressure 
    f_input.write("######################################################################\n\n")
    f_input.write("pair_style       lj/cut/coul/dsf 0.2 10.0 12.0\n")    #LJ = cut,  #coul = with DSF 
    f_input.write("pair_modify      tail yes\n\n")
    
    
    
    ## d- Writing pair coefficients for Protein (Epsilon and Sigma). We know the epsilon and sigma values between the element I and itself.
    ##    We apply the mixing rules to know epsilon and sigma between the element I and the element J. 
    ##
    ##    This section is organized in the following way:          
    ##
    ##    pair_coeff $i $i $epsilon $sigma  --- where  $i goes from 1 to the NUMBER OF ATOMTYPE (not atoms)      
    ##
    ##    eps = i[13], sigma = i[12] of "list_survived_atoms". 
    
    for count in range(1, number_atom_types + 1):
        for i in list_survived_atoms:
            if(count == i[-1]):
                f_input.write("pair_coeff {} {}   {:9.7f}   {:9.7f}\n".format(i[-1], i[-1], i[13]/4.184, i[12]*10))
                break
    
    ## e - Writing Pair Coefficients for Protein (Eps14 and Sigma14) between the element I and the element J. 
    ##     In particular, in case of charmm.ff, there is a section dedicated to pairtypes in "ffbonded.itp" file.
    ##     Therefore, the mixing rules for I and J are applied for only the elements not listed in pairtypes. 
    ##     In case of Amber.ff, pairtype section is not present, thus the mixing rules must be applied after writing 
    ##     the previus section (L16.d)   
    ##
    ##     This section is organized as follows:  
    ##
    ##     pair_coeff $i $j $epsilon14 $sigma14 -- whose values depend on the corresponding atomtypes of $I and $J listed in [ pairtypes ] 
    ##
    ##     The value of Sigma14 and Eps14 is not always listed for the atomtypes of $I and $J. In such case, mixing rules are applied 
    ##     Knowing Epsilon and Sigma for both elements $I and $J.
    
    for i in pairs_mult_res_prot:
        for j in list_survived_atoms:
            if(j[7] == "at" or j[7] =="AT" or j[7] =="At"):   # Writing only atomTypes and n_atomType for AT atoms as CG beads 
                                                              # does not have sigma14 and eps14 for charmm.ff and will be deleted later. 

                if(j[8] == i[2]):             # If the atomType is the same, we append the number of atomtype in pairs_mult_res_prot of the 1st element
                    i.append(j[-1])
                    break  
    
    for i in pairs_mult_res_prot:
        for j in list_survived_atoms:
            if(j[7] == "at" or j[7] =="AT" or j[7] =="At"):   # Writing only atomTypes and n_atomType for AT atoms as CG beads 
                                                              # does not have sigma14 and eps14 for charmm.ff and will be deleted later. 

                if(j[8] == i[3]):             # If the atomType is the same, we append the number of atomtype in pairs_mult_res_prot of the 2nd element 
                    i.append(j[-1])
                    break 

    
    pairs_ij = [x for x in pairs_mult_res_prot if len(x)>8]                                   # Remove pairs where no sigma14 and Eps14 are not present. 
    pairs_ij = [[x[1],x[0], x[3],x[2], x[5],x[4], x[7],x[6], x[9],x[8]] if x[8]>x[9] else x for x in pairs_ij]  # If n_AT_Type1 > n_AT_Type2: exchange values
    
    pairs_ij = sorted(pairs_ij, key=itemgetter(8,9))                                                         # sorting 9th-10th columns, i.e. the n_AT_Type
    pairs_ij =[list(my_iterator)[0] for g, my_iterator in itertools.groupby(pairs_ij, lambda x: [x[-2],x[-1]])]
    
    
    for i in pairs_ij: 
       f_input.write("pair_coeff {} {}   {:9.7f}   {:9.7f}\n".format(i[8], i[9], i[7]/4.184, i[6]*10))						
    
     
    ## f- Writing Mixing Rules
    f_input.write("\npair_modify      mix arithmetic\n\n")   #Using the mixing rules 
    
    
    ## g- In case a restart is required: TODO before restarting the simulation...
    
    f_input.write("# In case you use a restart, you need to copy all bonds and angles coefficient from .lmp file\n")
    f_input.write("# Please recall that before each row you need to write the string bond_coeff, angle_coeff\n")
    f_input.write("# and dihedral_coeff according with the case we are considering: Example.\n")
    f_input.write("# bond_coeff 1 harmonic  5.97514340e+01    3.95488306e+00    and so on\n\n")

    ## h- Writing other parameters: Creating Group Water (water) if H20 molecules are present. 
    if(solvate_Flag == True):  
        for k,v in dict_at_types.items():
            if(v == type_OW):
                n_attype_OW = k
    
            if(v == type_HW):
                n_attype_HW = k

    
    ## i- Writing other parameters: Creating fullyAT Protein group (AA) 
    list_attypes_AA = []       
    
    for i in list_survived_atoms:
        if(i[7] == "at" and  i[1] != "SOL" and i[1] !="NA" and i[1] !="CL"):
            list_attypes_AA.append(i[-1])
    
    list_attypes_AA = list(dict.fromkeys(list_attypes_AA))        # Removing duplicates keeping the sorting order  
    str_attypes_AA = ' '.join(str(i) for i in list_attypes_AA)    # Transforming the entire list in str in order to be written in the file as a singole row 
    
    
    
    ## j- Writing other parameters: Creating CG-Protein group (CG) 
    list_attypes_CG = []       
    
    for i in list_survived_atoms:
        if(i[7] == "cg"):
            list_attypes_CG.append(i[-1])
    
    list_attypes_CG = list(dict.fromkeys(list_attypes_CG))      
    str_attypes_CG = ' '.join(str(i) for i in list_attypes_CG)

    ## k- Writing last parameters: groups, fix parameters, thermo_style, timestep, and so on...
    if(solvate_Flag == True): 
        f_input.write("group            water type {} {}\n\n".format(n_attype_OW, n_attype_HW))
    
    if(list_attypes_AA):   # namely the list is not empty, that is the list has an element, at least. 
        f_input.write("group            AA    type {}\n\n".format(str_attypes_AA))
    
    if(list_attypes_CG):
        f_input.write("group            CG    type {}\n\n".format(str_attypes_CG))
    
    
    f_input.write("neighbor         1.0 bin\n")
    f_input.write("neigh_modify     every 1 delay 10 check yes\n\n")
    
    f_input.write("thermo_style     custom step temp press density vol ke ebond eangle evdwl ecoul etotal\n\n")
    
    f_input.write("thermo_modify    flush yes\n")
    
    f_input.write("thermo           ${Nf}\n\n")
    
    f_input.write("dump             trj all custom ${Ndump} mult-res.lammpstrj id type x y z\n")
    f_input.write("dump_modify      trj sort id\n\n")
    
    f_input.write("minimize         1.0e-4 1.0e-6 100 1000\n\n")
    
    for i in bonds_mult_res_prot:
        if((i[2]==type_OW and i[3]==type_HW) or (i[2]==type_HW and i[3]==type_OW)):
            bondtype_water = i[-1]
            break 
    
    for i in angles_mult_res_prot:
        if(i[3]==type_HW and i[4]==type_OW and i[5]==type_HW):
            angletype_water = i[-1]
    
    if(solvate_Flag == True):  
        f_input.write("fix 1            water shake 0.0001 20 0 t {} {} b {} a {}\n".format(n_attype_OW, n_attype_HW, bondtype_water, angletype_water))
        f_input.write("fix 2            all npt temp ${Text} ${Text} 100.0 iso ${Pext} ${Pext} 1000.0\n\n")
    else:
        f_input.write("fix 1            all npt temp ${Text} ${Text} 100.0 iso ${Pext} ${Pext} 1000.0\n\n")

    f_input.write("#fix 3           all nve\n")
    f_input.write("#fix 4           all langevin ${Text} ${Text} 100 9892571\n\n")
    f_input.write("restart          ${Nrestart} restart.${root}\n\n")
    f_input.write("timestep         1.0\n")
    f_input.write("run              ${Nrun}")
