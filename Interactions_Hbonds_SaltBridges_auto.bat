@echo off
for /F "tokens=* delims=\t" %%a in (case_subgroup.txt) do (
	echo %%a
	vmd -dispdev text -e Interactions_Hbonds_SaltBridges.tcl -args %%a
)