# Load trajectory
mol new ../3_WT/3_production/seokWT_onlyProt.psf type psf
mol addfile ../3_WT/3_production/seokWT_50ns_onlyProt.dcd type dcd first 0 last -1 step 2 waitfor all
mol addfile ../3_WT/3_production/seokWT_50to100ns_onlyProt.dcd type dcd first 0 last -1 step 2 waitfor all

# Align to first frame
set reference [atomselect top "protein and backbone and noh" frame 0]
		# the frame being compared
set compare [atomselect top "protein and backbone and noh"]
set num_steps [molinfo top get numframes]
set num [expr {$num_steps - 1}]
for {set frame 0} {$frame < $num} {incr frame} {
		# get the correct frame
		$compare frame $frame
		# compute the transformation
		set trans_mat [measure fit $compare $reference]
		# do the alignment
		$compare move $trans_mat
}
# Backbone RMSD
set outfile [open ../3_WT/3_production/rmsd_seokWT_100ns.dat w]
for {set frame 0} {$frame < $num_steps} {incr frame} {
		# get the correct frame
   $compare frame $frame
		# compute the transformation
   set trans_mat [measure fit $compare $reference]
		# do the alignment
   $compare move $trans_mat
		# compute the RMSD
   set rmsd [measure rmsd $compare $reference]
		# print the RMSD
   puts $outfile "$frame $rmsd"
}
close $outfile

# RMSF calculation
set outfile [open ../3_WT/3_production/rmsf_seokWT_100ns.dat w]
set sel [atomselect top "name CA"]
set rmsf [measure rmsf $sel first 0 last $num step 1]
for {set i 0} {$i < [$sel num]} {incr i} {
  puts $outfile "[expr {$i+1}] [lindex $rmsf $i]"
} 
close $outfile

# #################################################

exit

