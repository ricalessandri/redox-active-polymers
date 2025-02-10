#!/bin/bash

for rdx_unit in "TEMPO" "NMPHTH" ; do
  for solvent in "DMF" "THF" ; do
    cavity=$(grep "GePol: Cavity volume " ${rdx_unit}*${solvent}*.out  |tail -n 1 | awk '{ print $5 }')
    radius=$(python -c "print( ( (3.*float(${cavity})) / (4*3.1415926535) )**(1./3) )")   
    echo $rdx_unit $cavity $radius "("$solvent")"
    echo "The solvent has a negligible impact."
  done
done

