#!/bin/bash


polymer="PMAP"
phth_atoms="C08 N09 C10 O11 C12 C13 C14 C15 C16 C17 C18 O19 H10 H11 H12 H13"
solvent="DME"; cation="TBA";anion="PF6"
N="30" # degree of polymerization
TEMP="300"
replica="D"
SAMPLE="1-relax-100ns-whole.gro"

baseFOLDER=${PWD}

## for percent in "05" "10" "20" ; do
##   for charge in "000" "020" "060" ; do
for percent in "05" ; do
  for charge in "000" ; do

    FOLDER="../../configurations/${polymer}${charge}charge_${solvent}_${cation}${anion}_${percent}percent/relax-${N}mer-${TEMP}K-${replica}"
    
    cd ${FOLDER}
    
    # Generate index for all the following groups
    string_0="a ${phth_atoms} & r ${polymer}" # N-methyl-phthalimide selection string (all the atoms describing the N-methyl-phthalimide) 
    name_0="name 0 NMPHTHs"                   # N-methyl-phthalimide group name
    string_1="r TBA"                          # TBA selection string
    name_1="name 1 TBAs"                      # TBA group name
    
    for cmd in "del 0-20" "${string_0}" "${string_1}" "${name_0}" "${name_1}" "q"; do echo ${cmd}; done | gmx make_ndx -f ${SAMPLE} -o index_${charge}charge_${percent}percent.ndx

    rm \#index*
    
    cd ${baseFOLDER}

  done
done

