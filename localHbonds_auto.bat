@echo off
for /F "tokens=* delims=\t" %%a in (localHbonds_resid_case.txt) do (
	echo %%a
	vmd -dispdev text -e localHbonds_analyze.tcl -args %%a
)