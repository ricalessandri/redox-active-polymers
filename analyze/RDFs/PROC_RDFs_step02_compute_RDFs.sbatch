#!/bin/bash


polymer="PMAP"
#polymer="PEPP"
#polymer="PVBP"
solvent="DME"; cation="TBA";anion="PF6"
N="30" # degree of polymerization
TEMP="300"
replica="D"
SAMPLE="1-relax-100ns-whole.gro"
TOP_FILE="system_melt.top"
TPR_FILE="../local_topol.tpr"

dt=250    # default setting for gmx rdf
bin=02    # default setting for gmx rdf

currentFOLDER=${PWD}

#for percent in "05" "10" "20" ; do
#  for charge in "000" "020" "060" ; do
for percent in "05" ; do
  for charge in "000" ; do

    FOLDER="../../configurations/${polymer}${charge}charge_${solvent}_${cation}${anion}_${percent}percent/relax-${N}mer-${TEMP}K-${replica}"
    INDEX="index_${charge}charge_${percent}percent.ndx"

    echo "${TEMP}K"
    echo "${FOLDER}"
    cd    ${FOLDER}

    # Generate tpr file
    gmx grompp -p ${TOP_FILE} -c ${SAMPLE} -f ../../../MD_settings/eq_step1_min.mdp -o ${TPR_FILE}
    rm mdout.mdp

    # RDF between the centers of geometry of N-methyl-phthalimides
    for cmd in "NMPHTHs" "NMPHTHs"; do echo $cmd; done | gmx rdf -f ${SAMPLE} -n ${INDEX} -dt  ${dt} \
                                                         -selrpos part_res_cog -seltype part_res_cog -s ${TPR_FILE} -bin 0.${bin} -o rdf_NMPHTHs_dt0${dt}bin0${bin}${MODE}.xvg

    # TBA cation - N-methyl phthalimide center of geometry RDF
    for cmd in "NMPHTHs" "TBAs"; do echo $cmd; done | gmx rdf -f ${SAMPLE} -n ${INDEX} -dt  ${dt} \
                                                     -selrpos part_res_cog -seltype part_res_cog -s ${TPR_FILE} -bin 0.${bin} -o rdf_NMPHTHs-TBAs_dt0${dt}bin0${bin}${MODE}.xvg

    rm \#rdf*xvg*

    cd ${currentFOLDER}

  done
done

