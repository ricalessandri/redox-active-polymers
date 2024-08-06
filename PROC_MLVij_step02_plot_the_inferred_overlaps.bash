#!/bin/bash



# ** MONO007 - "monomers - 15x sets (~145k)" **
MODELlabel="MODELMONO007-monomers-2x600K595K585K575K475K465K455K445K435K345K335K325K315K305K"


polymer="PMAP"
#polymer="PEPP"
#polymer="PVBP"

solvent="DME";cation="TBA";anion="PF6"

SoC="000"
#SoC="020"
#SoC="060"
percent=05
#percent=10
#percent=20
N="30"
resname="${polymer}"

TEMP=300
SNAP="100ns"; TRAJSTEP="1"
CUTOFF="9"; cutoffLBL=$(printf "%02d" ${CUTOFF})



mkdir -p predictions
mkdir -p predictions/${MODELlabel}
cd       predictions/${MODELlabel}


# Just plot the predicted overlap distribution
for replica in "D"; do

  LABEL="${TEMP}K${replica}${SNAP}"

  cd ${polymer}${SoC}charge_${solvent}_${cation}${anion}_${percent}percent_${N}mer_${LABEL}/pair-predictions-${TEMP}K-${replica}-${SNAP}
  echo $pwd

  python3 ../../../../bin/just_plot_data_and_fit.py --what overlap --csvfile overlaps_predicted_cutoff${cutoffLBL}A.csv --v-type log --ylim 1.4
  mv 0_data_overlaps.pdf 0_data_overlaps_cutoff${cutoffLBL}A.pdf 
  
  cd ../../

done

