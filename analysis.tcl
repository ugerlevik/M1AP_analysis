set case [lindex $argv 0]
set subgroup [lindex $argv 1]

# #########################################
# DEFINITIONS
# Radius of gyration
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
#########################################

# #########################################
# REPEAT 1
# Load trajectory
mol new ../${subgroup}_${case}/3_production/seok${case}_onlyProt.psf type psf
mol addfile ../${subgroup}_${case}/3_production/seok${case}_250ns_onlyProt_each25k.dcd type dcd first 0 last -1 step 1 waitfor all

# Unwrap
package require pbctools
pbc unwrap -sel "protein"

# Number of H-bonds
package require hbonds
hbonds -sel1 [atomselect top "protein"] -writefile yes -plot no -outfile ../${subgroup}_${case}/3_production/hbonds_seok${case}.dat

# Salt bridges
package require saltbr
file mkdir ../${subgroup}_${case}/3_production/saltbridges_seok${case}
saltbr -sel [atomselect top protein] -log saltbridges_seok${case}.log -outdir ../${subgroup}_${case}/3_production/saltbridges_seok${case}

# Rg calculation
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

# SASA calculation
set sel [atomselect top "protein"]
set protein [atomselect top "protein"]
set outfile [open ../${subgroup}_${case}/3_production/SASA_seok${case}.dat w]
for {set i 0} {$i < $nf} {incr i} {
	molinfo top set frame $i
	set sasa [measure sasa 1.4 $protein -restrict $sel]
	puts $outfile "${i},${sasa}"
}
close $outfile

# FoldX Stability
mkdir foldx_stability_${case}
cd foldx_stability_${case}
for {set i 0} {$i <= $nf} {incr i 100} {
    [atomselect top all frame ${i}] writepdb ${i}.pdb
	set fd [open "${i}.pdb" r]
	set newfd [open "${i}.pdb.tmp" w]
	while {[gets $fd line] >= 0} {
		set newline [string map {HSD HIS HSE HIS HSP HIS} $line]
		puts $newfd $newline
	}
	close $fd
	close $newfd
	file rename -force "${i}.pdb.tmp" "${i}.pdb"
	exec foldx --command=Stability --pdb=${i}.pdb
}
set i [expr {${i} - 1}]
[atomselect top all frame ${i}] writepdb ${i}.pdb
set fd [open "${i}.pdb" r]
set newfd [open "${i}.pdb.tmp" w]
while {[gets $fd line] >= 0} {
	set newline [string map {HSD HIS HSE HIS HSP HIS} $line]
	puts $newfd $newline
}
close $fd
close $newfd
file rename -force "${i}.pdb.tmp" "${i}.pdb"
exec foldx --command=Stability --pdb=${i}.pdb

set outfile [open stability_seok${case}.dat w]
set i 0
while {$i < $nf} {
	puts $outfile "$i,[exec gawk {FNR==1 {print $2}} ${i}_0_ST.fxout]"
	set i [expr {${i} + 100}]
}
set i [expr {${i} - 1}]
puts $outfile "${i},[exec gawk {FNR==1 {print $2}} ${i}_0_ST.fxout]"
close $outfile
file rename stability_seok${case}.dat ../../${subgroup}_${case}/3_production/.
cd ..
file delete -force ./foldx_stability_${case}

# Backbone RMSD (and alignment to the first frame)
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

# RMSF calculation
set sel [atomselect top "name CA"]
set rmsf [measure rmsf $sel first 0 last $num step 1]
set resid [${sel} get resid]
set outfile [open ../${subgroup}_${case}/3_production/rmsf_seok${case}.dat w]
for {set i 0} {$i < [$sel num]} {incr i} {
	puts $outfile "[lindex ${resid} ${i}],[lindex ${rmsf} ${i}]"
}
close $outfile

# #################################################

# #################################################
# Clean Memory
mol delete all
# #################################################

# #################################################
# REPEAT 2
# Load trajectory
mol new ../${subgroup}_${case}/3_production/seok${case}_onlyProt.psf type psf
mol addfile ../${subgroup}_${case}/4_production_repeat2/seok${case}_250ns_r2_onlyProt_each25k.dcd type dcd first 0 last -1 step 1 waitfor all

# Unwrap
package require pbctools
pbc unwrap -sel "protein"

set nf [molinfo top get numframes]

# FoldX Stability
mkdir foldx_stability_${case}
cd foldx_stability_${case}
for {set i 0} {$i <= $nf} {incr i 100} {
    [atomselect top all frame ${i}] writepdb ${i}.pdb
	set fd [open "${i}.pdb" r]
	set newfd [open "${i}.pdb.tmp" w]
	while {[gets $fd line] >= 0} {
		set newline [string map {HSD HIS HSE HIS HSP HIS} $line]
		puts $newfd $newline
	}
	close $fd
	close $newfd
	file rename -force "${i}.pdb.tmp" "${i}.pdb"
	exec foldx --command=Stability --pdb=${i}.pdb
}
set i [expr {${i} - 1}]
[atomselect top all frame ${i}] writepdb ${i}.pdb
set fd [open "${i}.pdb" r]
set newfd [open "${i}.pdb.tmp" w]
while {[gets $fd line] >= 0} {
	set newline [string map {HSD HIS HSE HIS HSP HIS} $line]
	puts $newfd $newline
}
close $fd
close $newfd
file rename -force "${i}.pdb.tmp" "${i}.pdb"
exec foldx --command=Stability --pdb=${i}.pdb

set outfile [open stability_seok${case}.dat w]
set i 0
while {$i < $nf} {
	puts $outfile "$i,[exec gawk {FNR==1 {print $2}} ${i}_0_ST.fxout]"
	set i [expr {${i} + 100}]
}
set i [expr {${i} - 1}]
puts $outfile "${i},[exec gawk {FNR==1 {print $2}} ${i}_0_ST.fxout]"
close $outfile
file rename stability_seok${case}.dat ../../${subgroup}_${case}/4_production_repeat2/.
cd ..
file delete -force ./foldx_stability_${case}

# Backbone RMSD (and alignment to the first frame)
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

# RMSF calculation
set sel [atomselect top "name CA"]
set rmsf [measure rmsf $sel first 0 last $num step 1]
set resid [${sel} get resid]
set outfile [open ../${subgroup}_${case}/4_production_repeat2/rmsf_seok${case}.dat w]
for {set i 0} {$i < [$sel num]} {incr i} {
	puts $outfile "[lindex ${resid} ${i}],[lindex ${rmsf} ${i}]"
}
close $outfile

# #################################################

exit

