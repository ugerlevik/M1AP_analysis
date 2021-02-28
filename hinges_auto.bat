@echo off
for /F "tokens=* delims=\t" %%a in (case_subgroup.txt) do (
echo %%a
py hinges.py %%a
)