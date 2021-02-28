set mutname [lindex $argv 0]
set subgroup [lindex $argv 1]
set res [lindex $argv 2]

set wtname "WT"
set case "resid ${res}"

# Repeat 1 -------------------------------------------
# Load trajectory
mol new ../3_WT/3_production/seokWT_onlyProt.psf type psf
mol addfile ../3_WT/3_production/seokWT_250ns_onlyProt_each25k.dcd type dcd first 0 last -1 step 1 waitfor all

# Number of H-bonds
package require hbonds
hbonds -sel1 [atomselect top "same residue as (within 10 of ${case})"] -writefile yes -plot no -outfile ../_outputs/hbonds/hbonds10of${case}_wt.dat

# Load trajectory
mol delete all
mol new ../${subgroup}_${mutname}/3_production/seok${mutname}_onlyProt.psf type psf
mol addfile ../${subgroup}_${mutname}/3_production/seok${mutname}_250ns_onlyProt_each25k.dcd type dcd first 0 last -1 step 1 waitfor all

# Number of H-bonds
package require hbonds
hbonds -sel1 [atomselect top "same residue as (within 10 of ${case})"] -writefile yes -plot no -outfile ../_outputs/hbonds/hbonds10of${case}_${mutname}.dat


# Repeat 2 -------------------------------------------
# Load trajectory
mol delete all
mol new ../3_WT/3_production/seokWT_onlyProt.psf type psf
mol addfile ../3_WT/4_production_repeat2/seokWT_250ns_r2_onlyProt_each25k.dcd type dcd first 0 last -1 step 1 waitfor all

# Number of H-bonds
package require hbonds
hbonds -sel1 [atomselect top "same residue as (within 10 of ${case})"] -writefile yes -plot no -outfile ../_outputs/hbonds/hbonds10of${case}_wt_r2.dat

# Load trajectory
mol delete all
mol new ../${subgroup}_${mutname}/3_production/seok${mutname}_onlyProt.psf type psf
mol addfile ../${subgroup}_${mutname}/4_production_repeat2/seok${mutname}_250ns_r2_onlyProt_each25k.dcd type dcd first 0 last -1 step 1 waitfor all

# Number of H-bonds
package require hbonds
hbonds -sel1 [atomselect top "same residue as (within 10 of ${case})"] -writefile yes -plot no -outfile ../_outputs/hbonds/hbonds10of${case}_${mutname}_r2.dat


exit

