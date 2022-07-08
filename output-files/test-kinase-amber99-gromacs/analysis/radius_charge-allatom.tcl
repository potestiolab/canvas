set fp [open "radius_charge-allatom.txt" r]

# Count the number of lines in a text file

set number 0
while { [gets $fp line] >= 0 } {
    incr number
                               }
close $fp

############################################

set fp [open "radius_charge-allatom.txt" r]

set file_data [read $fp]

set data [split $file_data "\n"]

for {set i 1} {$i < $number} {incr i} {

     set a [lindex $data $i 0];
     set b [lindex $data $i 1];
     set c [lindex $data $i 2];
     set d [lindex $data $i 3];
     set e [atomselect top "serial $a"];

     $e set radius $c
     $e set beta $d
     		                      }
