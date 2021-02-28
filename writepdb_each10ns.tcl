set case [lindex $argv 0]
set subgroup [lindex $argv 1]

mol new ../${subgroup}_${case}/3_production/seok${case}_onlyProt.psf type psf
mol addfile ../${subgroup}_${case}/3_production/seok${case}_250ns_onlyProt_each25k.dcd type dcd first 0 last -1 step 1 waitfor all

set frameNumber [molinfo top get frame]
set i 0

set sel [atomselect 0 "protein" frame $i]
$sel writepdb ../_outputs/Hinge/${case}_prot_${i}ns.pdb

set i 499
while {$i <= $frameNumber} {
	set sel [atomselect 0 "protein" frame $i]
	$sel writepdb ../_outputs/Hinge/${case}_prot_[expr {$i + 1}]ns.pdb
	set i [expr {$i + 500}]
}

exit