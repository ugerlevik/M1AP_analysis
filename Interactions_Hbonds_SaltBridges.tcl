##################################################
## Project: 2021-M1AP
## Script purpose: H-bonds and salt bridges
## Date: September 17, 2021
## Author: Umut Gerlevik
##################################################

set case [lindex $argv 0]
set subgroup [lindex $argv 1]


# REPEAT 1 ---------------------------------------
# Load trajectory
mol new ../${subgroup}_${case}/3_production/seok${case}_onlyProt.psf type psf
mol addfile ../${subgroup}_${case}/3_production/seok${case}_250ns_onlyProt_each25k.dcd type dcd first 0 last -1 step 2 waitfor all
mol addfile ../${subgroup}_${case}/3_production/seok${case}_250to500ns_onlyProt_each25k.dcd type dcd first 0 last -1 step 2 waitfor all

# Number of H-bonds
package require hbonds
hbonds -sel1 [atomselect top "protein"] -writefile yes -plot no -outfile ../${subgroup}_${case}/3_production/hbonds_seok${case}.dat

# Salt bridges
package require saltbr
file mkdir ../${subgroup}_${case}/3_production/saltbridges_seok${case}
saltbr -sel [atomselect top protein] -log saltbridges_seok${case}.log -outdir ../${subgroup}_${case}/3_production/saltbridges_seok${case}


mol delete all

# REPEAT 2 ---------------------------------------
# Load trajectory
mol new ../${subgroup}_${case}/3_production/seok${case}_onlyProt.psf type psf
mol addfile ../${subgroup}_${case}/4_production_repeat2/seok${case}_250ns_r2_onlyProt_each25k.dcd type dcd first 0 last -1 step 2 waitfor all
mol addfile ../${subgroup}_${case}/4_production_repeat2/seok${case}_250to500ns_r2_onlyProt_each25k.dcd type dcd first 0 last -1 step 2 waitfor all

# Number of H-bonds
package require hbonds
hbonds -sel1 [atomselect top "protein"] -writefile yes -plot no -outfile ../${subgroup}_${case}/4_production_repeat2/hbonds_seok${case}.dat

# Salt bridges
package require saltbr
file mkdir ../${subgroup}_${case}/4_production_repeat2/saltbridges_seok${case}
saltbr -sel [atomselect top protein] -log saltbridges_seok${case}.log -outdir ../${subgroup}_${case}/4_production_repeat2/saltbridges_seok${case}


exit