#!/bin/bash
filename='case_subgroup.txt'
exec 4<$filename
echo Start
while read -u4 p ; do
    vmd -dispdev text -e analysis.tcl -args $p
done
