@echo off
for /F "tokens=* delims=\t" %%a in (Interactions_localHbonds_resid_case.txt) do (
	echo %%a
	vmd -dispdev text -e Interactions_localHbonds_analyze.tcl -args %%a
)