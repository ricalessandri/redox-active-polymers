#!/bin/bash


if [ "${HOSTNAME}" == "macbkpro2022ric.local" ] ; then
  echo "Don't forget to: "
  echo "conda activate ML"
fi

# ** MONO007 - "monomers - 15x sets (~145k)" **
BESTmodel="../NN/overlaps-MONOMERS007-monomers-2x600K595K585K575K475K465K455K445K435K345K335K325K315K305K/distlog100lr0.00100bs0512ne01000nn0400k5/trained_model_AA"
MODELlabel="MODELMONO007-monomers-2x600K595K585K575K475K465K455K445K435K345K335K325K315K305K"


polymer="PMAP"
#polymer="PEPP"
#polymer="PVBP"

solvent="DME";cation="TBA";anion="PF6"
N="30"
resname="${polymer}"

TEMP=300
SNAP="100ns"; TRAJSTEP="1"
CUTOFF="9"; cutoffLBL=$(printf "%02d" ${CUTOFF})

replica="D"

#for SoC in "000" "020" "060" ; do
for SoC in "000" ; do
#for percent in "05" "10" "20"; do
for percent in "05" ; do

LABEL="${TEMP}K${replica}${SNAP}"

FOLDER_LABEL="${solvent}_${cation}${anion}_${percent}percent"
FOLDERin="../configurations/${polymer}${SoC}charge_${FOLDER_LABEL}/relax-${N}mer-${TEMP}K-${replica}"

GRO_FILE="${FOLDERin}/1-relax-${SNAP}-whole.gro"
TOP_FILE="${FOLDERin}/system_melt.top"
TPR_FILE="${FOLDERin}/../local_topol.tpr"

POLYMERfolder="${polymer}${SoC}charge_${FOLDER_LABEL}_${N}mer_${LABEL}"

mkdir -p predictions
mkdir -p predictions/${MODELlabel}
mkdir -p predictions/${MODELlabel}/${POLYMERfolder}
cd       predictions/

# Generate tpr file
gmx grompp -p ${TOP_FILE} -c ${GRO_FILE} -f ../MD_settings/eq_step1_min.mdp -o ${TPR_FILE}
rm mdout.mdp

# Predict
python3 ../bin/predict_overlaps--pairs.py \
        -s ${TPR_FILE} -f ${GRO_FILE} -t ${TEMP} -l ${replica} -n ${SNAP} \
        --mapping AA --resname ${polymer} \
        --pdb-to-align-n ../bin/lib/NMPHTH-${polymer}-opt-neutral-wB97X.pdb \
        --pdb-to-align-r ../bin/lib/NMPHTH-${polymer}-opt-radical_anion-wB97X.pdb \
        --ML-model ${BESTmodel} \
        --map-file  ../mappings_${polymer}.json \
        --trajstep ${TRAJSTEP} \
        --cutoff ${CUTOFF}
        ##--test True \
        ##--verbose True \

mkdir -p                                                                    ${MODELlabel}/${POLYMERfolder}/pair_DATs
mv pair_DATs/pairs_info_cutoff${cutoffLBL}A_${TEMP}K_${replica}_${SNAP}.dat ${MODELlabel}/${POLYMERfolder}/pair_DATs/.
mkdir -p                                                                    ${MODELlabel}/${POLYMERfolder}/pair-predictions-${TEMP}K-${replica}-${SNAP}
mv pair-predictions-${TEMP}K-${replica}-${SNAP}/*                           ${MODELlabel}/${POLYMERfolder}/pair-predictions-${TEMP}K-${replica}-${SNAP}/.

# Clean-up
rm -r pair_DATs/ pair-predictions-${TEMP}K-${replica}-${SNAP}/ ${TPR_FILE}

cd ../

done
done

