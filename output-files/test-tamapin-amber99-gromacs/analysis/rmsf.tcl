proc rmsf_all_ca {} {

    set file_name "rmsf.dat"
    set outfile [open $file_name w];
    set sel [atomselect top "serial 2 6 10 14 18 22 26 30 34 38 42 46 50 53 55 59 63 67 71 76 98 108 127 134 149 164 174 195 199 203 207 "]
    set mol [$sel molindex]
    for {set i 0} {$i < [$sel num]} {incr i} {
        set rmsf [measure rmsf $sel]
        puts $outfile "[expr {$i+1}] [lindex $rmsf $i]"
        }
     close $outfile

     }
