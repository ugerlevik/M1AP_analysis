set mutname [lindex $argv 0]
set subgroup [lindex $argv 1]
set res [lindex $argv 2]

# Load trajectory
mol new ../3_WT/3_production/seokWT_onlyProt.psf type psf
mol addfile ../3_WT/3_production/seokWT_250ns_onlyProt_each25k.dcd type dcd first 0 last -1 step 1 waitfor all

proc uniqueList {list} {
  set new {}
  foreach item $list {
    if {[lsearch $new $item] < 0} {
      lappend new $item
    }
  }
  return $new
}

set a [atomselect top "within 10 of resid ${res} frame 1999"]

set outfile [open ../_outputs/saltbridges/localResidueList_${res}.txt w]
#puts $outfile "[uniqueList [$a get {resname resid chain}]]"
puts $outfile "[uniqueList [$a get {resname resid}]]"
close $outfile

exit