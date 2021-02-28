@echo off
for /F "tokens=* delims=\t" %%a in (localHbonds_resid_case.txt) do (
	echo %%a
	vmd -dispdev text -e LocalResidueID.tcl -args %%a
)