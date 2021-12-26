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

# Rg
proc center_of_mass {selection} {
        # some error checking
        if {[$selection num] <= 0} {
                error "center_of_mass: needs a selection with atoms"
        }
        # set the center of mass to 0
        set com [veczero]
        # set the total mass to 0
        set mass 0
        # [$selection get {x y z}] returns the coordinates {x y z} 
        # [$selection get {mass}] returns the masses
        # so the following says "for each pair of {coordinates} and masses,
	#  do the computation ..."
        foreach coord [$selection get {x y z}] m [$selection get mass] {
           # sum of the masses
           set mass [expr $mass + $m]
           # sum up the product of mass and coordinate
           set com [vecadd $com [vecscale $m $coord]]
        }
        # and scale by the inverse of the number of atoms
        if {$mass == 0} {
                error "center_of_mass: total mass is zero"
        }
        # The "1.0" can't be "1", since otherwise integer division is done
        return [vecscale [expr 1.0/$mass] $com]
}
proc gyr_radius {sel} {
  # make sure this is a proper selection and has atoms
  if {[$sel num] <= 0} {
    error "gyr_radius: must have at least one atom in selection"
  }
  # gyration is sqrt( sum((r(i) - r(center_of_mass))^2) / N)
  set com [center_of_mass $sel]
  set sum 0
  foreach coord [$sel get {x y z}] {
    set sum [vecadd $sum [veclength2 [vecsub $coord $com]]]
  }
  return [expr sqrt($sum / ([$sel num] + 0.0))]
}

set outfile [open ../${subgroup}_${case}/3_production/rog_seok${case}.dat w]
set nf [molinfo top get numframes] 
set i 0
while {$i < $nf} {
    set prot [atomselect top "protein" frame $i]
    set i [expr {$i + 1}]
    set rog [gyr_radius $prot]
    puts $outfile "${i},${rog}"
} 
close $outfile

mol delete all

# REPEAT 2 ---------------------------------------
# Load trajectory
mol new ../${subgroup}_${case}/3_production/seok${case}_onlyProt.psf type psf
mol addfile ../${subgroup}_${case}/4_production_repeat2/seok${case}_250ns_r2_onlyProt_each25k.dcd type dcd first 0 last -1 step 2 waitfor all
mol addfile ../${subgroup}_${case}/4_production_repeat2/seok${case}_250to500ns_r2_onlyProt_each25k.dcd type dcd first 0 last -1 step 2 waitfor all

set outfile [open ../${subgroup}_${case}/4_production_repeat2/rog_seok${case}.dat w]
set nf [molinfo top get numframes] 
set i 0
while {$i < $nf} {
    set prot [atomselect top "protein" frame $i]
    set i [expr {$i + 1}]
    set rog [gyr_radius $prot]
    puts $outfile "${i},${rog}"
} 
close $outfile



exit