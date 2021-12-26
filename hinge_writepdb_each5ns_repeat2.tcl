set case [lindex $argv 0]
set subgroup [lindex $argv 1]

mol new ../${subgroup}_${case}/3_production/seok${case}_onlyProt.psf type psf
mol addfile ../${subgroup}_${case}/4_production_repeat2/seok${case}_last75nsOf500ns_r2_onlyProt_each25k.dcd type dcd first 0 last -1 step 1 waitfor all

set frameNumber [molinfo top get frame]
set i 0

set sel [atomselect 0 "protein" frame $i]
$sel writepdb ../_outputs/Hinge_repeat2/${case}_prot_${i}ns.pdb

set i 99
while {$i <= $frameNumber} {
	set sel [atomselect 0 "protein" frame $i]
	$sel writepdb ../_outputs/Hinge_repeat2/${case}_prot_[expr {$i + 1}]ns.pdb
	set i [expr {$i + 100}]
}

exit