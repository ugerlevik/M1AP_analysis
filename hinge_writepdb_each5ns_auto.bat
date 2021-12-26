@echo off
for /F "tokens=* delims=\t" %%a in (case_subgroup.txt) do (
echo %%a
vmd -dispdev text -e hinge_writepdb_each5ns.tcl -args %%a
vmd -dispdev text -e hinge_writepdb_each5ns_repeat2.tcl -args %%a
)