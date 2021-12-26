##################################################
## Project: 2021-M1AP
## Script purpose: RMSD calculation
## Date: September 11, 2021
## Author: Umut Gerlevik
##################################################

set case [lindex $argv 0]
set subgroup [lindex $argv 1]


# REPEAT 1 ---------------------------------------
# Load trajectory
mol new ../${subgroup}_${case}/3_production/seok${case}_onlyProt.psf type psf
mol addfile ../${subgroup}_${case}/3_production/seok${case}_250ns_onlyProt_each25k.dcd type dcd first 0 last -1 step 2 waitfor all
mol addfile ../${subgroup}_${case}/3_production/seok${case}_250to500ns_onlyProt_each25k.dcd type dcd first 0 last -1 step 2 waitfor all

# Backbone RMSD (and alignment to the first frame)
set nf [molinfo top get numframes]
set reference [atomselect top "backbone" frame 0]
set compare [atomselect top "backbone"]
set num [expr {$nf - 1}]
set outfile [open ../${subgroup}_${case}/3_production/rmsd_seok${case}.dat w]
for {set frame 0} {$frame < $nf} {incr frame} {
	# get the correct frame
	$compare frame $frame
	# compute the transformation
	set trans_mat [measure fit $compare $reference]
	# do the alignment
	$compare move $trans_mat
	# compute the RMSD
	set rmsd [measure rmsd $compare $reference]
	# print the RMSD
	puts $outfile "${frame},${rmsd}"
}
close $outfile

mol delete all

# REPEAT 2 ---------------------------------------
# Load trajectory
mol new ../${subgroup}_${case}/3_production/seok${case}_onlyProt.psf type psf
mol addfile ../${subgroup}_${case}/4_production_repeat2/seok${case}_250ns_r2_onlyProt_each25k.dcd type dcd first 0 last -1 step 2 waitfor all
mol addfile ../${subgroup}_${case}/4_production_repeat2/seok${case}_250to500ns_r2_onlyProt_each25k.dcd type dcd first 0 last -1 step 2 waitfor all

# Backbone RMSD (and alignment to the first frame)
set nf [molinfo top get numframes]
set reference [atomselect top "backbone" frame 0]
set compare [atomselect top "backbone"]
set num [expr {$nf - 1}]
set outfile [open ../${subgroup}_${case}/4_production_repeat2/rmsd_seok${case}.dat w]
for {set frame 0} {$frame < $nf} {incr frame} {
	# get the correct frame
	$compare frame $frame
	# compute the transformation
	set trans_mat [measure fit $compare $reference]
	# do the alignment
	$compare move $trans_mat
	# compute the RMSD
	set rmsd [measure rmsd $compare $reference]
	# print the RMSD
	puts $outfile "${frame},${rmsd}"
}
close $outfile



exit