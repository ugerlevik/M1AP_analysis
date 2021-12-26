##################################################
## Project: 2021-M1AP
## Script purpose: Rg calculation
## Date: September 12, 2021
## Author: Umut Gerlevik
##################################################

set case [lindex $argv 0]
set subgroup [lindex $argv 1]


# REPEAT 1 ---------------------------------------
# Load trajectory
mol new ../${subgroup}_${case}/3_production/seok${case}_onlyProt.psf type psf
mol addfile ../${subgroup}_${case}/3_production/seok${case}_250ns_onlyProt_each25k.dcd type dcd first 0 last -1 step 2 waitfor all
mol addfile ../${subgroup}_${case}/3_production/seok${case}_250to500ns_onlyProt_each25k.dcd type dcd first 0 last -1 step 2 waitfor all

# SASA
set nf [molinfo top get numframes] 
set sel [atomselect top "protein"]
set protein [atomselect top "protein"]
set outfile [open ../${subgroup}_${case}/3_production/SASA_seok${case}.dat w]
for {set i 0} {$i < $nf} {incr i} {
	molinfo top set frame $i
	set sasa [measure sasa 1.4 $protein -restrict $sel]
	puts $outfile "${i},${sasa}"
}
close $outfile

mol delete all

# REPEAT 2 ---------------------------------------
# Load trajectory
mol new ../${subgroup}_${case}/3_production/seok${case}_onlyProt.psf type psf
mol addfile ../${subgroup}_${case}/4_production_repeat2/seok${case}_250ns_r2_onlyProt_each25k.dcd type dcd first 0 last -1 step 2 waitfor all
mol addfile ../${subgroup}_${case}/4_production_repeat2/seok${case}_250to500ns_r2_onlyProt_each25k.dcd type dcd first 0 last -1 step 2 waitfor all

# SASA
set nf [molinfo top get numframes] 
set sel [atomselect top "protein"]
set protein [atomselect top "protein"]
set outfile [open ../${subgroup}_${case}/4_production_repeat2/SASA_seok${case}.dat w]
for {set i 0} {$i < $nf} {incr i} {
	molinfo top set frame $i
	set sasa [measure sasa 1.4 $protein -restrict $sel]
	puts $outfile "${i},${sasa}"
}
close $outfile


exit