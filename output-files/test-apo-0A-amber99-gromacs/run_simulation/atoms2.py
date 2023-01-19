"""
   This script reads the new topology file called "mult.itp" and it returns, in ascending order from 1 to N, 
   the CA index, and its nature, i.e. atomistic, medium-grained, or coarse-grained.  
   The script is very useful because in the calculation of RMSF the CA index of the x-axis 
   goes between 1 and the number N of CAs (o residues of course). 

   o In case of all-atom simulation, each CA index is the residue number 

   o In case of AT-CG simulation, each CA index is still the residue number since the residue number of each CA_cg
     is the same of the reference fullyAT   

   o In case of AT-MG-CG simulation, each bead in the MG region (N,CA,C,O) has its own residue number 
     even tough in the all-atom representation they belong to the same residue. Therefore the enumeration 
     undergoes variations. This script is useful, in particular, in this last case. Look the next example: 

         - The 1st, 2nd, and 8th residues are described in MG, the remainder as CG:  

   1 N
   2 CA
   3 C
   4 O
-----------
   5 N
   6 CA
   7 C
   8 O
-----------
   9 CA
-----------
  10 CA
-----------
  11 CA
-----------
  12 CA 
-----------
  13 CA 
-----------
  14 N
  15 CA
  16 C
  17 O
----------

   Taking only the CA index we have that 2,6,15 are CA_mg, while 9,10,11,12,13 are CA_cg. 
   We cannot take these value. Thus, sorting, and starting from 1, in the end, we have that:

   MG = [1,2,8]
   CG = [3,4,5,6,7]
"""

f  = open("mult.itp", "r") 
f2 = open("CA_regions.txt", "w")

AT        = []
MG        = []
MG_carbon = []
CG        = []

AT2 = []
MG2 = []
CG2 = []

for line in f:

    line =line.strip()

    if(line =="[ atoms ]"):
        break                        # We found the "atoms" part; we break and start again. 


for line in f:

    if(line =="\n"):                 # If the row is empty, then we exit the loop because we have read all the atoms in topology file...
        break

    if(line[0] !=";"):               # If the 1st letter of line is ";" then skip the line itself, otherwise we split it. 
        splt = line.split()

        at_number  = int(splt[0])
        at_type    = str(splt[1])    # atomistic type; very important 
        res_number = int(splt[2])
        res_name   = str(splt[3])
        at_name    = str(splt[4])
        cgnr       = int(splt[5])
        charge     = float(splt[6])
        mass       = float(splt[7])


        if(res_name != "MUL"):
            if(at_name == "CA"): 
                AT.append(at_number)

        if(res_name == "MUL"): 
            if(mass == 14.010): 
                MG.append(at_number+1)        	       # Because the order is always N-CA-C-O, thus after finding the Nitrogen N of MUL residue,
                MG_carbon.append(at_number+2)	       # at_number+1 and at_number+2 correspond at CA and C, respectively.   
 
            if(mass == 12.010):                        # We take all possible C (both CA and C of medium-grained regions and CA of coarse-grained  
                CG.append(at_number) 		       # region. In the next operation, only CA_cg will be kept   



CG = sorted(list(set(CG)-set(MG)-set(MG_carbon)))      # In this way CG region mantain only CA_cg, since the CA_mg (MG region) and C_mg (MG_carbon region)
						       # will be exlcuded after insiemistic operation: 
                                                       # CG (CA_cg) = CG (CA_cg + CA_mg + C_mg) - MG (CA_cg)  - MG_carbon (C_mg) 

max_regions =[max(AT), max(MG), max(CG)]

maximum = max(max_regions) 

i = 0 
for k in range(1, maximum+1): 

    if k in AT:
        i = i + 1
        AT2.append(i) 
        f2.write("{:s} {:d}\n".format("at", i))

    if k in MG:
        i = i + 1 
        MG2.append(i)
        f2.write("{:s} {:d}\n".format("mg", i))

    if k in CG:
        i = i + 1 
        CG2.append(i)
        f2.write("{:s} {:d}\n".format("cg", i))

print("\n")
print("AT = ",AT2 )
print("\n")
print("MG = ",MG2)
print("\n")
print("CG = ", CG2)
print("\n")
