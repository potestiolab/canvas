import math 
import time 
import os 
import itertools
from operator import itemgetter

# Function that rewrite the solvated coordinate file after using "gmx solvate" command in GROMACS applied to only survived atoms 
# The purpose is to change the res_name, res_name, at_name, and at_number of the protein part if the corresponding survived atom  is CG. 
# The atomistic part, and the SOL one will be not affected by any variation.
def write_newsolvated_file(f, atoms, list_survived_atoms, n_solvated_atoms, Lx, Ly, Lz):  
    max = 0
    
    for j in atoms:
        if(j[1]!="SOL"):
            for i in list_survived_atoms: 
    
                if(i[9]==j[3]):# and j[1]!="SOL"):    # if at_number is the same in both lists: 
    
                    if(i[7] == "cg"):
                        j[0] = i[17]
                        j[1] = "MUL"
                        j[2] = i[16]
                        j[3] = i[9]
    
                    if(i[7] == "at"):
                        j[0] = i[17]
                        j[3] = i[9]
    
            last_res_number = i[17]
    
    
    f.write("Protein in Variable resolution, new indexes for Gromacs simulation\n")
    f.write("{:5d}\n".format(n_solvated_atoms))
        
    count = 1 
    
    for i in atoms:
        kk = math.floor(i[3]/100000)
        if(i[1]=="SOL"): 
            res = math.ceil(count/3)+last_res_number          # res increases by 3 every three for SOL 
            if(res<=99999):  
                f.write("{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(res,i[1],i[2],i[3]-kk*100000,i[4],i[5],i[6]))
            else:
                last_res_number = -1 
                count = 1 
                res = 0 
                f.write("{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(res,i[1],i[2],i[3]-kk*100000,i[4],i[5],i[6]))
    
            count = count + 1
    
        else:
            f.write("{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(i[0],i[1],i[2],i[3]-kk*100000,i[4],i[5],i[6]))
     
    f.write("{:10.5f}{:10.5f}{:10.5f}".format(Lx/10,Ly/10,Lz/10))     # divided by 10 because gromacs requires nm, not Amgstrom. 

    

# Function that writes only the coordinate file (.gro) of the atoms that survive ready for being solvated 
def write_outsurv_file(f, n_survived, list_survived_atoms, Lx, Ly, Lz):
    f.write("Protein in Variable Resolution\n")
    f.write("{:5d}\n".format(n_survived))

    for i in list_survived_atoms:
        kk = math.floor(i[3]/100000) 
        f.write("{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(i[9], i[1], i[2], i[3]-kk*100000, i[4], i[5], i[6]))

    f.write("{:10.5f}{:10.5f}{:10.5f}".format(Lx/10,Ly/10,Lz/10))    # divided by 10, because the module read_gromacs.py read coordinates in Angstrom
    							             # GROMACS, on the other hand, requires nm.



# Function that writes the file containing nonbonded parameters (epsilon and sigma) for each CG bead.
def write_ffnonbondednew_file(list_survived_atoms): 
    fnb_new = open("ffnonbonded_new.itp", "w")
    
    fnb_new.write("[ atomtypes ]\n")
    fnb_new.write("; name      at.num  mass     charge ptype  sigma      epsilon\n")
    
    charg = 0.0
    ptype = "A"
     
    for i in list_survived_atoms:
        if(i[7] == "cg"):        # Write the new atomtypes only for CG beads. 
            fnb_new.write("{:6s}{:>6d}{:10.4f}{:10.4f}{:>3s}   {:10.5e}{:15.5e}\n".format(i[16], i[15], i[10], charg, ptype, i[12], i[13]))
    


# Function that writes the coordinates of the initial frame (.gro), keeping only the the atoms that survive leaving position and name unchanged 
def write_initial_frame_file(f, n_survived, list_survived_atoms, Lx, Ly, Lz):  
    f.write("Initial frame for protein in Variable Resolution\n")
    f.write("{:5d}\n".format(n_survived))

    for i in list_survived_atoms: 
        if(i[7]=="cg"):
            f.write("{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(i[-1], "MUL", i[16], i[9], i[4], i[5], i[6])) 
        else:
            f.write("{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(i[-1], i[1], i[2], i[9], i[4], i[5], i[6]))   #i[2], NOT [16]

    f.write("{:10.5f}{:10.5f}{:10.5f}\n".format(Lx/10,Ly/10,Lz/10))



# Function that creates the new file of position restraints in the CANVAS model 
def write_posre_new_file(list_survived_atoms): 
    f_posre = open("posre_mul_res.itp", "w")

    f_posre.write("; In this topology include file, you will find position restraint\n")
    f_posre.write("; entries for all the heavy atoms in your original pdb file.\n")
    f_posre.write("; This means that all the protons which were added by pdb2gmx are\n")
    f_posre.write("; not restrained.\n\n")


    f_posre.write("[ position_restraints ]\n")
    f_posre.write("; atom  type      fx      fy      fz\n")

    typee = 1
    fx    = 1000
    fy    = 1000
    fz    = 1000

    for i in list_survived_atoms:
        if(i[10] > 2.0):                 # If mass is higher than 2.0, i.e. higher than the Hydrogen mass. 
            f_posre.write("{:>6d}{:>6d}{:6d}{:6d}{:6d}\n".format(i[9],typee, fx, fy, fz))



# Function that write the new file for bonded interaction parameters (bond- angle- dihedral-types) between atoms, where 
# at least one CG bead is involved (the original ffbonded.itp file reports the parameters between AT-AT atoms in the original forcefield
def write_ffbonded_new_file(f, bonds_mult_res_prot, angles_mult_res_prot, dihedrals_mult_res_prot):  

    # II.a- Bondtypes: in particular, here we use a different approach, at difference with angletypes 
    #       and dihedraltypes. The reason stems from the fact that a CG bead (e.g. A1) could be bonded, at the same time, at 
    #       two different Carbons (e.g CA and CB). These two carbons have the same atomtype (e.g. CT). There in bondtypes, 
    #       we would write two times "A1 CT b0 Kb". It is also evident that b0 will differ is both case and gromacs cannot 
    #       recognize which couple of atoms we are considering. 
    #       The solution is provided by making use of a string that identifies the specific bond in the topology file 
    #       (as we see later) that follows the bond itself. Example:   
    #
    #       1  2  Gb_$i where $i in an INTEGER.
    # 
    #       Then, in [bondtypes] (in this section) we write:
    #
    #       #define Gb_$i   b0  kb     ;func is not necessary since it is read in [bond] 
    #
    #       It is like for torsion part in dihedrals.
    
    
    k = 1
    
    for i in bonds_mult_res_prot:
        if(i[4] == "cg" or i[5] == "cg"):
            f.write(";\n")
            f.write(";new bondtypes with specific b0 and kb\n")
            f.write(";\n\n")
            break
    
    
    for i in bonds_mult_res_prot:
        if(i[4] == "cg" or i[5] == "cg"):
            f.write("#define Gb_{:d} {:11.4f}{:15.3f}\n".format(k, i[7], i[6]))
            k = k + 1
    
    # II.b- Angletypes: they will be written in standard format. 
    
    #func_angle = 1
    
    for i in angles_mult_res_prot:
        if(i[-3] == "cg" or i[-2] == "cg" or i[-1] == "cg"):   #If a CG bead is present in "angles_mult_res_prot", we write the angletypes section.
            f.write("\n[ angletypes ]\n")
            f.write(";  i    j    k  func       th0       cth\n")
            break
    
    
    for i in angles_mult_res_prot:
        if(i[-3] == "cg" or i[-2] == "cg" or i[-1] == "cg"):
            if(i[3]==5):  # If func=5:
                f.write("{:>6s}{:>6s}{:>6s}{:>6d}{:11.4f}{:15.3f}{:15.3f}{:15.3f}\n".format(i[4], i[5], i[6], i[3], i[8], i[7], i[9], i[10]))
            elif(i[3]==1 or i[3]==2):  # if func=1 or func=2    
                 f.write("{:>6s}{:>6s}{:>6s}{:>6d}{:11.4f}{:15.3f}\n".format(i[4], i[5], i[6], i[3], i[8], i[7]))
    
    # II.c- Dihedraltypes with funct = 4: they will be written in standard format. 
    
    func_dih4 = 4
    
    for i in dihedrals_mult_res_prot:
        if(i[4] == 4):
            if(i[12] == "cg" or i[13] == "cg" or i[14] == "cg" or i[15] == "cg"):
                f.write("\n[ dihedraltypes ]\n")
                f.write(";i  j   k  l     func      phase      kd      pn\n")
                break
    
    for i in dihedrals_mult_res_prot:
        if(i[4]==4):
            if(i[12] == "cg" or i[13] == "cg" or i[14] == "cg" or i[15] == "cg"):
                f.write("{:>6s}{:>6s}{:>6s}{:>6s}{:>6d}{:11.4f}{:15.5f}{:>6d}\n".format(i[5], i[6], i[7], i[8], func_dih4, i[9], i[10], i[11]))
    
    
    # II.d- Dihedraltypes with funct = 9: as for bondtypes, we use the same method, since some quartet of atoms 
    #       is repetead more than once. Therefore, in order to be safe, we use the string Dh_{$i} where $i is a number. 
    # 
    #       In case of torsion, we do not write anything, since they are completely atomistic are they are already treated 
    #       in ffbonded.itp file of Amber99sb-ildn forcefield.
    
    func_dih9 = 9 
    
    k = 1
    
    for i in dihedrals_mult_res_prot:
        if(i[4] == 9):
            if(not i[5].startswith("torsion")):
                if(i[12] == "cg" or i[13] == "cg" or i[14] == "cg" or i[15] == "cg"):
                    f.write("\n;\n")
                    f.write(";new dihedraltype9 with specific phase, kd and mult(n)\n")
                    f.write(";\n\n")
                    break
    
    for i in dihedrals_mult_res_prot:
        if(i[4] == 9):
            if(not i[5].startswith("torsion")):
                if(i[12] == "cg" or i[13] == "cg" or i[14] == "cg" or i[15] == "cg"):
                    f.write("#define Dh9_{:d} {:11.4f}{:15.5f}{:>6d}\n".format(k, i[9], i[10], i[11]))
                    k = k + 1


    # II.e- Dihedraltypes with funct = 2: they will be written in standard format, since there is no multiplicity, and thus 
    #       each quartet of atoms is written only once (as for dih4).  

    func_dih2 = 2

    for i in dihedrals_mult_res_prot:
        if(i[4] == 2):
            if(i[-4] == "cg" or i[-3] == "cg" or i[-2] == "cg" or i[-1] == "cg"):
                f.write("\n[ dihedraltypes ]\n")
                f.write(";i  j   k  l     func      phi0       kphi\n")
                break

    for i in dihedrals_mult_res_prot:
        if(i[4] == 2):
            if(i[-4] == "cg" or i[-3] == "cg" or i[-2] == "cg" or i[-1] == "cg"):
                f.write("{:>6s}{:>6s}{:>6s}{:>6s}{:>6d}{:11.4f}{:15.5f}\n".format(i[5], i[6], i[7], i[8], func_dih2, i[9], i[10]))
 

# Function that write the file named mult.itp that is essentially the topology file for our system that includes bonds, angles, dihedrals (and cmaps only 
# for charmm.ff) 

def write_mult_itp_file(f, list_survived_atoms, bonds_mult_res_prot, pairs_mult_res_prot, angles_mult_res_prot, dihedrals_mult_res_prot, \
                        k_strong, cmaps_mult_res_prot):

    # 1- Writing intro 
     
    f.write("[ moleculetype ]\n")
    f.write("; Name            nrexcl\n")
    f.write("Mult_Res_Prot       3\n\n")
    
    f.write("[ atoms ]\n")
    f.write(";   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB\n")
       
    for i in list_survived_atoms:
        if(i[7]=="cg"):
            f.write("{:6d}{:>11s}{:7d}{:>7s}{:>7s}{:7d}{:11.4f}{:11.3f}\n".format(i[9], i[16], i[17], "MUL", i[16], i[9], i[11], i[10]))
        else:
            f.write("{:6d}{:>11s}{:7d}{:>7s}{:>7s}{:7d}{:11.4f}{:11.3f}\n".format(i[9], i[16], i[17], i[1], i[2], i[9], i[11], i[10]))
    
    
    # 2- Writing bonds: we have just to pay attention at the fuction of bonds. 
    #    If an atom is AT, and the second one is CG, it means that we are on the interface, 
    #    therefore, in order to avoid problems of exclusions in between CG beads and AT part, 
    #    in this particular case, the bond function is 6, and not 1.
    #    However, in case of covalent bond, we prefer to keep the original chemical bond type, i.e. 1 
    # 
    #    If func = 6 means that there is an harmonic bond between CG and AT whose functional form 
    #    is always kb(b-b0)^2, but there no real spring, it is a "fake" spring whose purpose is to 
    #    avoid Coulomb and LJ exclusion between CG beads and some atomistic atoms. 
    
    
    func_bond  = 1
    func_harm  = 6 
    
    f.write("\n[ bonds ]\n")
    f.write(";  ai    aj funct            c0            c1            c2            c3\n")
    
    
    k = 1 
    
    for i in bonds_mult_res_prot:
    
        if(i[4] == "cg" and i[5] == "cg"):
            f.write("{:>5d}{:>6d}{:>6d}    Gb_{:d}\n".format(i[0], i[1], func_bond, k))
            k = k + 1
    
        elif(i[4] == "cg" and i[5] == "at"):
            if(i[6] <= k_strong):   # The weak bond has a variable value, thus a bond is covalent if knb > k_strong 
                f.write("{:>5d}{:>6d}{:>6d}    Gb_{:d}\n".format(i[0], i[1], func_harm, k))
                k = k + 1
            else: 
                f.write("{:>5d}{:>6d}{:>6d}    Gb_{:d}\n".format(i[0], i[1], func_bond, k))
                k = k + 1
    
        elif(i[4] == "at" and i[5] == "cg"):
            if(i[6] <= k_strong): 								 
                f.write("{:>5d}{:>6d}{:>6d}    Gb_{:d}\n".format(i[0], i[1], func_harm, k))
                k = k + 1
            else: 
                f.write("{:>5d}{:>6d}{:>6d}    Gb_{:d}\n".format(i[0], i[1], func_bond, k)) 
                k = k + 1
    
        elif(i[4] == "at" and i[5] == "at"):
            f.write("{:>5d}{:>6d}{:>6d}\n".format(i[0], i[1], func_bond))
    
    
    # 3- Writing pairs (we have already excluded pairs for one or both CG beads in section B.5) 
    
    func_pair = 1
    
    f.write("\n[ pairs ]\n")
    f.write(";  ai    aj funct            c0            c1            c2            c3\n")
    
    
    for i in pairs_mult_res_prot:
        f.write("{:>5d}{:>6d}{:>6d}\n".format(i[0], i[1], func_pair))
    
    
    
    # 4- Writing angles 
    
    #func_angle = 1
    
    f.write("\n[ angles ]\n")
    f.write(";  ai    aj    ak funct            c0            c1            c2            c3\n")
    
    for i in angles_mult_res_prot:
        f.write("{:>5d}{:>6d}{:>6d}{:>6d}\n".format(i[0], i[1], i[2], i[3]))    #i[3]==func_angle
    
    
    # 5- Writing dihedrals with funct = 9, included torsion 
    #    The charmm forcefield is employed if, at least, one of two conditions are fullfilled:
    #        o "dih2" are found in original topology file 
    #        o "cmap" are found in original topology file      
    #    If charmm.ff or Amber.ff are employed  in the CANVAS topology file (mult.itp) only once the quartet
    #    of atoms, even if the multiplicity is higher than one. The following command means that we remove duplicated dihedrals  
    #    # namely the first four elements ([x[0],x[1],x[2],x[3]]) of a longer list of lists... Look also A.5.5 in "CANVAS.py" script


    # - Splitting in dih_CG, dih_AT (for funct=2,4, and 9 excluded torsions), and dih_tors  

    dih2     = []
    dih4_AT  = []
    dih4_CG  = [] 
    dih9_AT  = []
    dih9_CG  = []
    dih_tors = []    


    for i in dihedrals_mult_res_prot:

        if(i[4]==2):
            dih2.append(i) 

        if(i[4]==4):
            if(i[12]=="cg" or i[13]=="cg" or i[14]=="cg" or i[15]=="cg"):
                dih4_CG.append(i)
            else:
                dih4_AT.append(i) 

        if(i[4]==9 and not i[5].startswith("torsion")):
            if(i[12]=="cg" or i[13]=="cg" or i[14]=="cg" or i[15]=="cg"):
                dih9_CG.append(i) 
            else:
                dih9_AT.append(i)

        if(i[4]==9 and i[5].startswith("torsion")):
            dih_tors.append(i) 


    # - Removing repeated dihedrals for the AT part only in Dih9 and Dih4, since in the CG part, afterwards, will be created a string Dh_{$i}    
    #   In the all_atom part, the dihedrals indexes (a1,a2,a3,a4) must be written only once even if the multiplicity is more than one. 
    #   Then according the multiplicity "n", gromacs reads automatically "n"-times from the ffbonded.itp file, section "dihedrals"

    dih4_AT =[list(my_iterator)[0] for g, my_iterator in itertools.groupby(dih4_AT, lambda x: [x[0],x[1],x[2],x[3]])]
    dih9_AT =[list(my_iterator)[0] for g, my_iterator in itertools.groupby(dih9_AT, lambda x: [x[0],x[1],x[2],x[3]])]
  

    # - Creating again "dihedrals_mult_res_prot" making the sum of dih2 + dih4_AT + dih4_CG + dih9_AT + dih9_CG + dih_tors
    #   and then sorting 

    dihedrals_mult_res_prot = dih2 + dih4_AT + dih4_CG + dih9_AT + dih9_CG + dih_tors

    dihedrals_mult_res_prot = sorted(dihedrals_mult_res_prot, key=itemgetter(0,1,2,3,-1))  # Sorting a1, a2, a3, a4, and finally the multiplicity n.  
                                                                                           # Note that in case of dih2, "-1" is kb, not multiplicity 
                                                                                           # but it is not a problem since we sort first in terms of 
                                                                                           # a1, a2, a3, a4; dih2 has no "n", therefore there is no 
                                                                                           # possibility to have identical a1,a2,a3,a4. Sorting based 
                                                                                           # on "-1" is done only for identical a1,a2,a3,a4 for dih4 and 
                                                                                           # dih9 having different multiplicity.  



    """
    Flag_charmm = False """

    func_dih9 = 9
    
    """                 
    ## 5.1 - Is charmm forcefield empoyed?
    for i in dihedrals_mult_res_prot:
        if(i[4] == 2):  
            Flag_charmm = True 
            dihedrals_mult_res_prot=[list(my_iterator)[0] for g, my_iterator in itertools.groupby(dihedrals_mult_res_prot, lambda x: [x[0],x[1],x[2],x[3]])]
            break 

    if(Flag_charmm == False):     # If dih2 are not found, let us give a look to cmap 
        if(Flag_cmaps == True):   # cmap are found, i.e. the charmm.ff is employed
            Flag_charmm = True
            dihedrals_mult_res_prot=[list(my_iterator)[0] for g, my_iterator in itertools.groupby(dihedrals_mult_res_prot, lambda x: [x[0],x[1],x[2],x[3]])]
    """

    

    ## 5.2 - Writing dihedral_9 
    for i in dihedrals_mult_res_prot:
        if(i[4] == 9):
            f.write("\n[ dihedrals ]\n")
            f.write(";  ai    aj    ak    al funct            c0            c1            c2            c3\n")
            break
    
    k = 1
    
    for i in dihedrals_mult_res_prot:
        if(i[4] == 9):                        # in case of  dihedrals 9 and torsion (special dihedrals with funct equals to 9)  
            if(i[5].startswith("torsion")):
                f.write("{:>5d}{:>6d}{:>6d}{:>6d}{:>6d}    {:s}\n".format(i[0], i[1], i[2], i[3], func_dih9, i[5]))
    
            else:  				  # If it is just dih9, no torsion and if it is involved, at least, one CG bead.
    
                if(i[12] =="cg" or i[13] =="cg" or i[14]=="cg" or i[15]=="cg"):   
                    f.write("{:>5d}{:>6d}{:>6d}{:>6d}{:>6d}    Dh9_{:d}\n".format(i[0], i[1], i[2], i[3], func_dih9, k))
                    k = k + 1
                else:
                    f.write("{:>5d}{:>6d}{:>6d}{:>6d}{:>6d}\n".format(i[0], i[1], i[2], i[3], func_dih9))
    
    
    # 6- Writing dihedrals with funct = 4 
    
    func_dih4 = 4

    for i in dihedrals_mult_res_prot:
        if(i[4] == 4):
            f.write("\n[ dihedrals ]\n")
            f.write(";  ai    aj    ak    al funct            c0            c1            c2            c3\n")
            break
    
    for i in dihedrals_mult_res_prot:
        if(i[4] == 4):
            f.write("{:>5d}{:>6d}{:>6d}{:>6d}{:>6d}\n".format(i[0], i[1], i[2], i[3], func_dih4))
    

    # 7- Writing dihedrals with funct = 2 

    func_dih2 = 2

    for i in dihedrals_mult_res_prot:
        if(i[4] == 2):
            f.write("\n[ dihedrals ]\n")
            f.write(";  ai    aj    ak    al funct            c0            c1            c2            c3\n")
            break


    for i in dihedrals_mult_res_prot:
        if(i[4] == 2):
            f.write("{:>5d}{:>6d}{:>6d}{:>6d}{:>6d}\n".format(i[0], i[1], i[2], i[3], func_dih2))


    # 8- Writing cmaps (Only for charmm.ff) 
    
    if cmaps_mult_res_prot:  # If CMAP list is NOT-empty (i.e. this list exists):
        f.write("\n[ cmap ]\n") 
        f.write(";  ai    aj    ak    al    am funct\n")

    if cmaps_mult_res_prot:
        for i in cmaps_mult_res_prot:
            f.write("{:>5d}{:>6d}{:>6d}{:>6d}{:>6d}{:>6d}\n".format(i[0], i[1], i[2], i[3], i[4], i[5]))
            

# Function that writes the new cmap file (Torsional Correction Map) if a charmm.ff is employed and if a CG bead is involved 
def write_new_cmap_file(f, cmaps_mult_res_prot): 
   
    f.write("[ cmaptypes ]\n\n")

    for i in cmaps_mult_res_prot:
        if(i[-5] =="cg" or i[-4] =="cg" or i[-3] =="cg" or i[-2]=="cg" or i[-1]=="cg"):
            f.write("{} {} {} {} {} {} {} {}\n".format(i[6], i[7], i[8], i[9], i[10], i[5], i[11], i[12]))

            for count in range(13,71): 
                f.write("{}".format(i[count])) 

            f.write("\n")
     

# Function that writes the new topology file for simulating system in CANVAS model 
def write_new_topology_file(f, path_ff_new, count_SOL, Flag_cmaps):

    timee  = time.strftime('%a %d %B %Y - %H:%M:%S')
    
    f.write(";\n")
    f.write(";       File 'topol_new.top' was generated\n")
    f.write(";       By user: USERNAME\n")#.format(os.getlogin()))  				# username = os.getlogin()
    f.write(";       On host: {}\n".format(os.uname()[1]))  				# hostname = os.uname()[1]
    f.write(";       At date: {}\n".format(timee))
    f.write(";\n")
    f.write(";       This is the topology file for Multiple resolution protein\n")
    f.write(";\n\n")
    
    f.write("; Include forcefield parameters\n")
    f.write('#include "{}forcefield.itp"\n'.format(path_ff_new))
    
    f.write('#include "ffnonbonded_new.itp"\n')
    f.write('#include "ffbonded_new.itp"\n')

    if(Flag_cmaps == True): 
        f.write('#include "cmap_new.itp"\n')

    f.write('#include "mult.itp"\n\n')
    
    f.write("; Include Position restraint file\n")
    f.write("#ifdef POSRES\n")
    f.write('#include "posre_mul_res.itp"\n')
    f.write("#endif\n\n")
    
    f.write("; Include water topology\n")
    f.write('#include "{}tip3p.itp"\n\n'.format(path_ff_new))
    
    f.write("#ifdef POSRES_WATER\n")
    f.write("; Position restraint for each water oxygen\n")
    f.write("[ position_restraints ]\n")
    f.write(";  i funct       fcx        fcy        fcz\n")
    f.write("   1    1       1000       1000       1000\n")
    f.write("#endif\n\n")
    
    f.write("; Include topology for ions\n")
    f.write('#include "{}ions.itp"\n\n'.format(path_ff_new))
    
    f.write("[ system ]\n")
    f.write("; Name\n")
    f.write("Multiple resolution protein in fully-at water\n\n")
    
    f.write("[ molecules ]\n")
    f.write("; Compound        #mols\n")
    f.write("Mult_Res_Prot         1\n")
    f.write("SOL         {}".format(count_SOL))          # count_SOL is obtained in Section 7.3 
     

# Function that writes the ions.mdp file for adding ions to a system for neutralizing charge
def write_ions_mdp(f): 
    f.write("; ions.mdp - used as input into grompp to generate ions.tpr\n\n")
    
    f.write("; Parameters describing what to do, when to stop and what to save\n")
    f.write("integrator      = steep         ; Algorithm (steep = steepest descent minimization)\n")
    f.write("emtol           = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm\n")
    f.write("emstep          = 0.01          ; Minimization step size\n")
    f.write("nsteps          = 50000         ; Maximum number of (minimization) steps to perform\n\n")
    
    f.write("; Parameters describing how to find the neighbors of each atom and how to calculate the interactions\n")
    f.write("nstlist         = 1             ; Frequency to update the neighbor list and long range forces\n")
    f.write("cutoff-scheme   = Verlet        ; Buffered neighbor searching\n")
    f.write("ns_type         = grid          ; Method to determine neighbor list (simple, grid)\n")
    f.write("coulombtype     = cutoff        ; Treatment of long range electrostatic interactions\n")
    f.write("rcoulomb        = 1.0           ; Short-range electrostatic cut-off (any value is fine)\n")
    f.write("rvdw            = 1.0           ; Short-range Van der Waals cut-off (any value is fine)\n")
    f.write("pbc             = xyz           ; Periodic Boundary Conditions in all 3 dimensions\n")


# Function that writes the em.mdp file for doing a steepest descent algorithm for minimation
def write_em_mdp(f, sig): 
    f.write("; minim.mdp - used as input into grompp to generate em.tpr\n\n")
    
    f.write("; Parameters describing what to do, when to stop and what to save\n")
    f.write("integrator      = steep                 ; Algorithm (steep = steepest descent minimization)\n")
    f.write("emtol           = 1000.0                ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm\n")
    f.write("emstep          = 0.01                  ; Minimization step size\n")
    f.write("nsteps          = 50000                 ; Maximum number of (minimization) steps to perform\n\n")
    
    f.write("; Parameters describing how to find the neighbors of each atom and how to calculate the interactions\n")
    f.write("nstlist         = 10                    ; Frequency to update the neighbor list and long range forces\n")
    f.write("cutoff-scheme   = Verlet                ; Buffered neighbor searching\n")
    f.write("ns_type         = grid                  ; Method to determine neighbor list (simple, grid)\n")
    f.write("coulombtype     = reaction-field        ; Treatment of long range electrostatic interactions\n")
    f.write("rcoulomb        = {:5.3f}                 ; Short-range electrostatic cut-off (nm) equals to 2.5*max_sigma\n".format(2.5*sig))
    f.write("rvdw            = {:5.3f}                 ; Short-range Van der Waals cut-off (nm) equals to 2.5*max_sigma\n".format(2.5*sig))
    f.write("pbc             = xyz                   ; Periodic Boundary Conditions in all 3 dimensions\n")
 

# Function that writes the npt.mdp file for doing a 50 ps equilibration in NVT ensamble 
def write_nvt_mdp(f, sig):
    f.write("title                   = NVT equilibration        ;\n\n")
    f.write("define                  = -DPOSRES                 ; position restrain the protein\n\n")
    
    f.write("; Run parameters\n")
    f.write("integrator              = md                       ; leap-frog integrator\n")
    f.write("nsteps                  = 25000                    ; 2 * 25000 = 50 ps\n")
    f.write("dt                      = 0.002                    ; 2 fs\n\n")
    
    f.write("; Output control\n")
    f.write("nstxout                 = 1000                     ; save coordinates every 1.0 ps\n")
    f.write("nstvout                 = 1000                     ; save velocities every 1.0 ps\n")
    f.write("nstenergy               = 1000                     ; save energies every 1.0 ps\n")
    f.write("nstlog                  = 1000                     ; update log file every 1.0 ps\n")
    f.write("xtc_grps                = Protein                  ;\n")
    f.write("nstxtcout               = 1000                     ;\n\n")
    
    f.write("; Bond parameters\n")
    f.write("continuation            = no                       ; first dynamics run\n")
    f.write("constraint_algorithm    = lincs                    ; holonomic constraints\n")
    f.write("constraints             = h-bonds                  ; bonds involving H are constrained\n")
    f.write("lincs_iter              = 1                        ; accuracy of LINCS\n")
    f.write("lincs_order             = 4                        ; also related to accuracy\n\n")
    
    f.write("; Nonbonded settings\n")
    f.write("cutoff-scheme           = Verlet                   ; Buffered neighbor searching\n")
    f.write("ns_type                 = grid                     ; search neighboring grid cells\n")
    f.write("nstlist                 = 10                       ; 20 fs, largely irrelevant with Verlet\n")
    f.write("rcoulomb                = {:5.3f}                    ; short-range electrostatic cutoff (nm) equals to 2.5*max_sigma\n".format(2.5*sig))
    f.write("rvdw                    = {:5.3f}                    ; short-range Van der Waals cutoff (nm) equals to 2.5*max_sigma\n".format(2.5*sig))
    f.write("DispCorr                = EnerPres                 ; account for cut-off vdW scheme\n\n")
    
    f.write("; Electrostatics\n")
    f.write("coulombtype             = reaction-field           ; Particle Mesh Ewald for long-range electrostatics\n")
    f.write("epsilon-rf              = 80                       ;\n")
    f.write(";pme_order              = 4                        ; cubic interpolation\n")
    f.write(";fourierspacing         = 0.16                     ; grid spacing for FFT\n\n")
    
    f.write("; Temperature coupling is on\n")
    f.write("tcoupl                  = V-rescale                ; modified Berendsen thermostat\n")
    f.write("tc-grps                 = Protein    Non-Protein   ; two coupling groups - more accurate\n")
    f.write("tau_t                   = 0.1        0.1           ; time constant, in ps\n")
    f.write("ref_t                   = 300        300           ; reference temperature, one for each group, in K\n\n")
    
    f.write("; Pressure coupling is off\n")
    f.write("pcoupl                  = no                       ; no pressure coupling in NVT\n\n")
    
    f.write("; Periodic boundary conditions\n")
    f.write("pbc                     = xyz                      ; 3-D PBC\n\n")
    
    f.write("; Velocity generation\n")
    f.write("gen_vel                 = yes                      ; assign velocities from Maxwell distribution\n")
    f.write("gen_temp                = 300                      ; temperature for Maxwell distribution\n")
    f.write("gen_seed                = -1                       ; generate a random seed\n")
    

# Function that writes the file run.mdp for doing a 50 ps equilibration in NPT ensemble  
def write_npt_mdp(f, sig):
    f.write("title                   = NPT equilibration       ;\n\n")
    f.write("define                  = -DPOSRES                ; position restrain the protein\n\n")
    
    f.write("; Run parameters\n")
    f.write("integrator              = md                      ; leap-frog integrator\n")
    f.write("nsteps                  = 25000                   ; 2 * 25000 = 50 ps\n")
    f.write("dt                      = 0.002                   ; 2 fs\n\n")
    
    f.write("; Output control\n")
    f.write("nstxout                 = 1000                    ; save coordinates every 1.0 ps\n")
    f.write("nstvout                 = 1000                    ; save velocities every 1.0 ps\n")
    f.write("nstenergy               = 1000                    ; save energies every 1.0 ps\n")
    f.write("nstlog                  = 1000                    ; update log file every 1.0 ps\n")
    f.write("xtc_grps                = Protein                 ;\n")
    f.write("nstxtcout               = 1000                    ;\n\n")
    
    f.write("; Bond parameters\n")
    f.write("continuation            = yes                     ; Restarting after NVT\n")
    f.write("constraint_algorithm    = lincs                   ; holonomic constraints\n")
    f.write("constraints             = h-bonds                 ; bonds involving H are constrained\n")
    f.write("lincs_iter              = 1                       ; accuracy of LINCS\n")
    f.write("lincs_order             = 4                       ; also related to accuracy\n\n")
    
    f.write("; Nonbonded settings\n")
    f.write("cutoff-scheme           = Verlet                  ; Buffered neighbor searching\n")
    f.write("ns_type                 = grid                    ; search neighboring grid cells\n")
    f.write("nstlist                 = 10                      ; 20 fs, largely irrelevant with Verlet scheme\n")
    f.write("rcoulomb                = {:5.3f}                   ; short-range electrostatic cutoff (nm) equals to 2.5*max_sigma\n".format(2.5*sig))
    f.write("rvdw                    = {:5.3f}                   ; short-range van der Waals cutoff (nm) equals to 2.5*max_sigma\n".format(2.5*sig))
    f.write("DispCorr                = EnerPres                ; account for cut-off vdW scheme\n\n")
    
    f.write("; Electrostatics\n")
    f.write("coulombtype             = reaction-field          ; Particle Mesh Ewald for long-range electrostatics\n")
    f.write("epsilon-rf              = 80                      ;\n")
    f.write(";pme_order              = 4                       ; cubic interpolation\n")
    f.write(";fourierspacing         = 0.16                    ; grid spacing for FFT\n\n")
    
    f.write("; Temperature coupling is on\n")
    f.write("tcoupl                  = V-rescale               ; modified Berendsen thermostat\n")
    f.write("tc-grps                 = Protein   Non-Protein   ; two coupling groups - more accurate\n")
    f.write("tau_t                   = 0.1       0.1           ; time constant, in ps\n")
    f.write("ref_t                   = 300       300           ; reference temperature, one for each group, in K\n\n")
    
    f.write("; Pressure coupling is on\n")
    f.write("pcoupl                  = Parrinello-Rahman       ; Pressure coupling on in NPT\n")
    f.write("pcoupltype              = isotropic               ; uniform scaling of box vectors\n")
    f.write("tau_p                   = 2.0                     ; time constant, in ps\n")
    f.write("ref_p                   = 1.0                     ; reference pressure, in bar\n")
    f.write("compressibility         = 4.5e-5                  ; isothermal compressibility of water, bar^-1\n")
    f.write("refcoord_scaling        = com\n\n")
    
    f.write("; Periodic boundary conditions\n")
    f.write("pbc                     = xyz                     ; 3-D PBC\n\n")
    
    f.write("; Velocity generation\n")
    f.write("gen_vel                 = no                      ; Velocity generation is off\n")


# Function that write the file run.mdp for doing a long simulating run
def write_run_mdp(f,sig): 
    f.write("title                   = run\n\n")

    f.write("; Run parameters\n")
    f.write("integrator              = md                      ; leap-frog integrator\n")
    f.write("nsteps                  = 250000000               ; 2 * 250000000 = 500'000 ps (500 ns)\n")
    f.write("dt                      = 0.002                   ; 2 fs\n\n")

    f.write("; Output control\n")
    f.write("nstxout                 = 0                       ; suppress bulky .trr file by specifying\n")
    f.write("nstvout                 = 0                       ; 0 for output frequency of nstxout\n")
    f.write("nstfout                 = 0                       ; nstvout, and nstfout\n")
    f.write("nstenergy               = 10000                   ; save energies every 10.0 ps\n")
    f.write("nstlog                  = 10000                   ; update log file every 10.0 ps\n")
    f.write("nstxout-compressed      = 10000                   ; save compressed coordinates every 10.0 ps\n")
    f.write("compressed-x-grps       = Protein                 ; save only Protein coordinate (incl. MUL res. --> residuetypes.dat file required)\n\n")

    f.write("; Bond parameters\n")
    f.write("continuation            = yes                     ; Restarting after NPT\n")
    f.write("constraint_algorithm    = lincs                   ; holonomic constraints\n")
    f.write("constraints             = h-bonds                 ; bonds involving H are constrained\n")
    f.write("lincs_iter              = 1                       ; accuracy of LINCS\n")
    f.write("lincs_order             = 4                       ; also related to accuracy\n\n")

    f.write("; Neighborsearching\n")
    f.write("cutoff-scheme           = Verlet                  ; Buffered neighbor searching\n")
    f.write("ns_type                 = grid                    ; search neighboring grid cells\n")
    f.write("nstlist                 = 10                      ; 20 fs, largely irrelevant with Verlet scheme\n")
    f.write("rcoulomb                = {:5.3f}                   ; short-range electrostatic cutoff (nm) equals to 2.5*max_sigma\n".format(2.5*sig))
    f.write("rvdw                    = {:5.3f}                   ; short-range Van der Waals cutoff (nm) equals to 2.5*max_sigma\n\n".format(2.5*sig))

    f.write("; Electrostatics\n")
    f.write("coulombtype              = reaction-field         ; Particle Mesh Ewald for long-range electrostatics\n")
    f.write("epsilon-rf               = 80                     ;\n")
    f.write(";pme_order               = 4                      ; cubic interpolation\n")
    f.write(";fourierspacing          = 0.16                   ; grid spacing for FFT\n\n")

    f.write("; Temperature coupling is on \n")
    f.write("tcoupl                  = V-rescale               ; modified Berendsen thermostat\n")
    f.write("tc-grps                 = Protein   non-Protein   ; two coupling groups - more accurate\n")
    f.write("tau_t                   = 0.1       0.1           ; time constant, in ps\n")
    f.write("ref_t                   = 300       300           ; reference temperature, one for each group, in K\n\n")

    f.write("; Pressure coupling is on\n")
    f.write("pcoupl                  = Parrinello-Rahman       ; Pressure coupling on in NPT\n")
    f.write("pcoupltype              = isotropic               ; uniform scaling of box vectors\n")
    f.write("tau_p                   = 2.0                     ; time constant, in ps\n")
    f.write("ref_p                   = 1.0                     ; reference pressure, in bar\n")
    f.write("compressibility         = 4.5e-5                  ; isothermal compressibility of water, bar^-1\n\n")

    f.write("; Periodic boundary conditions\n")
    f.write("pbc                     = xyz                     ; 3-D PBC\n\n")

    f.write("; Dispersion correction\n")
    f.write("DispCorr                = EnerPres                ; account for cut-off vdW scheme\n\n")

    f.write("; Velocity generation\n")
    f.write("gen_vel                 = no                      ; Velocity generation is off\n\n")
