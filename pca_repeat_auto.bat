@echo off
for /F "tokens=* delims=\t" %%a in (pca_cases.txt) do (
	echo %%a
	py pca_repeat.py %%a
)