set case [lindex $argv 0]
set subgroup [lindex $argv 1]

mol new ../${subgroup}_${case}/3_production/seok${case}.psf type psf
mol addfile ../${subgroup}_${case}/3_production/seok${case}.pdb

set sel [atomselect top "protein"]

$sel writepdb ../${subgroup}_${case}/3_production/seok${case}_onlyProt.pdb
#$sel writepsf pca_${dcdname}_cMD100ns_CA.psf

exit