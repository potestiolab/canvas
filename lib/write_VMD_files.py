# Function that writes a file in which, for each atom type (atomistic, medium-grained, and coarse-grained), 
# is indicated  the Van der Waals radius (sigma/2) annd the charge.  
def write_radius_charges(f_rc, list_survived_atoms): 

    f_rc.write("; index  at_type  VdW_radius (AA)    charge\n")

    for i in list_survived_atoms:
        f_rc.write("{:6d}    {:4s}  {:10.5f}        {:10.5f}\n".format(i[9], i[16], 10*i[12]/2, i[11]))     # 10*sigma/2 because sigma is a diameter, 
                                                                                                            # but VMD works with Angstrom and not nm. 

# Function that writes the .tkl script that reads "f_rc"  for visualizing radius and charge coloration in VMD   
def write_tcl_radius_charges(f_tcl_rc): 
 
    f_tcl_rc.write('set fp [open "radius_charge-allatom.txt" r]\n\n')

    f_tcl_rc.write('# Count the number of lines in a text file\n\n')

    f_tcl_rc.write("set number 0\n") 
    f_tcl_rc.write("while { [gets $fp line] >= 0 } {\n")
    f_tcl_rc.write("    incr number\n")
    f_tcl_rc.write("                               }\n")  							   
    f_tcl_rc.write("close $fp\n\n")

    f_tcl_rc.write('############################################\n\n')

    f_tcl_rc.write('set fp [open "radius_charge-allatom.txt" r]\n\n')

    f_tcl_rc.write('set file_data [read $fp]\n\n')

    f_tcl_rc.write('set data [split $file_data "\\n"]\n\n')

    f_tcl_rc.write('for {set i 1} {$i < $number} {incr i} {\n\n')

    f_tcl_rc.write('     set a [lindex $data $i 0];\n')
    f_tcl_rc.write('     set b [lindex $data $i 1];\n')
    f_tcl_rc.write('     set c [lindex $data $i 2];\n')
    f_tcl_rc.write('     set d [lindex $data $i 3];\n')

    f_tcl_rc.write('     set e [atomselect top "serial $a"];\n\n')

    f_tcl_rc.write('     $e set radius $c\n')
    f_tcl_rc.write('     $e set beta $d\n')
    f_tcl_rc.write('     		                      }\n')



# Function that writes the .tkl script that reads the CA list and creates a file for performing the RMSF in VMD of only CA atoms 
# (atomistic, medium-grained, and coarse-grained). 
def write_rmsf(f_tcl_rmsf, CA_list):

    f_tcl_rmsf.write("proc rmsf_all_ca {} {\n\n")
    f_tcl_rmsf.write('    set file_name "rmsf.dat"\n')
    f_tcl_rmsf.write("    set outfile [open $file_name w];\n")

    f_tcl_rmsf.write('    set sel [atomselect top "serial ')

    for i in CA_list:
        f_tcl_rmsf.write("{} ".format(i))

    f_tcl_rmsf.write('"]\n')

    f_tcl_rmsf.write("    set mol [$sel molindex]\n")
    f_tcl_rmsf.write("    for {set i 0} {$i < [$sel num]} {incr i} {\n")
    f_tcl_rmsf.write("        set rmsf [measure rmsf $sel]\n")
    f_tcl_rmsf.write('        puts $outfile "[expr {$i+1}] [lindex $rmsf $i]"\n')
    f_tcl_rmsf.write("        }\n")
    f_tcl_rmsf.write("     close $outfile\n\n")
    f_tcl_rmsf.write("     }\n")
