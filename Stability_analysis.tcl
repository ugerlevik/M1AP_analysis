##################################################
## Project: 2021-M1AP
## Script purpose: Rg calculation
## Date: September 15, 2021
## Author: Umut Gerlevik
##################################################

set case [lindex $argv 0]
set subgroup [lindex $argv 1]


# REPEAT 1 ---------------------------------------
# Load trajectory
mol new ../${subgroup}_${case}/3_production/seok${case}_onlyProt.psf type psf
mol addfile ../${subgroup}_${case}/3_production/seok${case}_250ns_onlyProt_each25k.dcd type dcd first 0 last -1 step 2 waitfor all
mol addfile ../${subgroup}_${case}/3_production/seok${case}_250to500ns_onlyProt_each25k.dcd type dcd first 0 last -1 step 2 waitfor all

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
file rename stability_seok${case}.dat ../../${subgroup}_${case}/3_production/.
cd ..
file delete -force ./foldx_stability_${case}

mol delete all

# REPEAT 2 ---------------------------------------
# Load trajectory
mol new ../${subgroup}_${case}/3_production/seok${case}_onlyProt.psf type psf
mol addfile ../${subgroup}_${case}/4_production_repeat2/seok${case}_250ns_r2_onlyProt_each25k.dcd type dcd first 0 last -1 step 2 waitfor all
mol addfile ../${subgroup}_${case}/4_production_repeat2/seok${case}_250to500ns_r2_onlyProt_each25k.dcd type dcd first 0 last -1 step 2 waitfor all

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



exit