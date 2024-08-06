#!/bin/bash


# ** MONO007 - "monomers - 15x sets (~145k)" **
MODELlabel="MODELMONO007-monomers-2x600K595K585K575K475K465K455K445K435K345K335K325K315K305K"


polymer="PMAP"
#polymer="PEPP"
#polymer="PVBP"

solvent="DME";cation="TBA";anion="PF6"
N="30"

TEMP=300
SNAP="100ns"; TRAJSTEP="1"
CUTOFF="9"; cutoffLBL=$(printf "%02d" ${CUTOFF})

replica="D"

SoC="000"
percent="05"


LABEL="${TEMP}K${replica}${SNAP}"

FOLDER_LABEL="${solvent}_${cation}${anion}_${percent}percent"
#FOLDERin="../configurations/${polymer}${SoC}charge_${FOLDER_LABEL}/relax-${N}mer-${TEMP}K-${replica}"

POLYMERfolder="${polymer}${SoC}charge_${FOLDER_LABEL}_${N}mer_${LABEL}"

SLOPE="0.8091639975519355"      # from [https://doi.org/10.1021/jacsau.4c00276]
INTERCEPT="-1.5699211355668465" # from [https://doi.org/10.1021/jacsau.4c00276]

cd predictions/${MODELlabel}/${POLYMERfolder}/pair-predictions-${TEMP}K-${replica}-${SNAP}

python ../../../../bin/log10overlap_to_ECeV.py --overlaps overlaps_predicted_cutoff09A.csv --verbose False --slope ${SLOPE} --intercept ${INTERCEPT}

FILE_1='log_couplings.csv'
NAME="PMAP"

python ../../../../bin/just_plot_couplings_vs_state_of_charge.py --what coupling      --csvfile01 ${FILE_1} --polymer ${NAME} --percent ${percent}
mv 0_data_couplings.pdf 0_data_couplings_cutoff09A.pdf
mv 0_data_couplings.txt 0_data_couplings_cutoff09A.txt

cd ../../../../

