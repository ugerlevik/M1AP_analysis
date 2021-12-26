##################################################
## Project: 2021-M1AP
## Script purpose: RMSF calculation
## Date: September 12, 2021
## Author: Umut Gerlevik
##################################################

set case [lindex $argv 0]
set subgroup [lindex $argv 1]


# REPEAT 1 ---------------------------------------
# Load trajectory
mol new ../${subgroup}_${case}/3_production/seok${case}_onlyProt.psf type psf
mol addfile ../${subgroup}_${case}/3_production/seok${case}_250to500ns_onlyProt_each25k.dcd type dcd first 3499 last -1 step 1 waitfor all

# Align
set nf [molinfo top get numframes]
set reference [atomselect top "backbone" frame 0]
set compare [atomselect top "backbone"]
set num [expr {$nf - 1}]
for {set frame 0} {$frame < $nf} {incr frame} {
	# get the correct frame
	$compare frame $frame
	# compute the transformation
	set trans_mat [measure fit $compare $reference]
	# do the alignment
	$compare move $trans_mat
}

# RMSF calculation
set sel [atomselect top "name CA"]
set rmsf [measure rmsf $sel first 0 last $num step 1]
set resid [${sel} get resid]
set outfile [open ../${subgroup}_${case}/3_production/rmsf_seok${case}.dat w]
for {set i 0} {$i < [$sel num]} {incr i} {
	puts $outfile "[lindex ${resid} ${i}],[lindex ${rmsf} ${i}]"
}
close $outfile

mol delete all

# REPEAT 2 ---------------------------------------
# Load trajectory
mol new ../${subgroup}_${case}/3_production/seok${case}_onlyProt.psf type psf
mol addfile ../${subgroup}_${case}/4_production_repeat2/seok${case}_250to500ns_r2_onlyProt_each25k.dcd type dcd first 3499 last -1 step 1 waitfor all

# Align
set nf [molinfo top get numframes]
set reference [atomselect top "backbone" frame 0]
set compare [atomselect top "backbone"]
set num [expr {$nf - 1}]
for {set frame 0} {$frame < $nf} {incr frame} {
	# get the correct frame
	$compare frame $frame
	# compute the transformation
	set trans_mat [measure fit $compare $reference]
	# do the alignment
	$compare move $trans_mat
}

# RMSF calculation
set sel [atomselect top "name CA"]
set rmsf [measure rmsf $sel first 0 last $num step 1]
set resid [${sel} get resid]
set outfile [open ../${subgroup}_${case}/4_production_repeat2/rmsf_seok${case}.dat w]
for {set i 0} {$i < [$sel num]} {incr i} {
	puts $outfile "[lindex ${resid} ${i}],[lindex ${rmsf} ${i}]"
}
close $outfile



exit