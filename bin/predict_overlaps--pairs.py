#!/usr/bin/env python3
# coding: utf-8
"""
Prediction of orbital overlaps based on feature vectors of *pairs* of monomers.
"""

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os, sys
import datetime
from MDAnalysis.analysis import distances as mda_dist
import time
from pathlib import Path
import argparse
import json
import src.functions as srcfunctions 
from tensorflow import keras
from sklearn import preprocessing
import csv
import logging
from contextlib import redirect_stdout
import joblib

# Check that I'm using python3
if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")


# Just plot formatting
plt.style.use('seaborn-talk')
plt.rcParams['font.family'] = 'sans'
plt.rcParams['font.size'] = 18
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.labelweight'] = 'normal'
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['legend.fontsize'] = 18   # For CONCISE legends (paper)
plt.rcParams['legend.fontsize'] = 14   # For VERBOSE legends
plt.rcParams['figure.titlesize'] = 18 



# ================================================
# =               INPUT DATA                     =
# ================================================
# Parse the arguments
parser = argparse.ArgumentParser(description='Generation of feature vectors and Gaussian inputs from an atomistic morphology for *pairs* of monomers.')
parser.add_argument('-f', '--gro-file', required=True, type=str  , help='name of GRO file of the morphology')
parser.add_argument('-s', '--tpr-file', required=True, type=str  , help='name of TPR file of the morphology')
parser.add_argument('-t', '--temp'    , required=True, type=int  , help='temperature at which the sample was taken (e.g., "300")') 
parser.add_argument('-l', '--label'   , required=True, type=str  , help='label identifying the run (e.g., "A"')
parser.add_argument('-n', '--snap'    , required=True, type=str  , help='time at which the snapshot was taken (e.g., "49ns")')
parser.add_argument('-c', '--cutoff'  , required=True, type=int  , help='cutoff for pair selection')
parser.add_argument('-r', '--resname' , required=True, type=str  , help='name of the residue to be analyzed')
parser.add_argument('--map-file'      , required=True, type=str  , help='name of JSON file containing mappings') # OLDER default='all_mappings_PTMA_39atoms.json'
parser.add_argument('--pdb-to-align-n', required=True, type=str  , help='name of PDB file to be aligned')
parser.add_argument('--pdb-to-align-r', required=True, type=str  , help='name of PDB file to be aligned')
parser.add_argument('--verbose'       , default=False,             help='True if you want verbose output')
parser.add_argument('--also-cm'       , default=False,             help='True if you want *also* Coulomb Matrices') 
parser.add_argument('--also-1d'       , default=False,             help='True if you want *also* 1D flattened matrices containing the ecli/stag identities')
parser.add_argument('--test'          , default=False,             help='True if you want to run a test')
parser.add_argument('--trajstep'      , default=1    , type=int  , help='Step size for trajectory; default = 1 = read all frames')
parser.add_argument('--v-type'        , default='log', type=str  , help='type of transformation to apply to couplings; *abs*, *log*, or *signed*')
parser.add_argument('--feature'       , default='distmat', type=str, help='type of input featurization: distmat OR coulmat')
parser.add_argument('--make-plot'     , default=False , type=bool , help='toggles on/off plotting')
parser.add_argument('--mapping'       , default=False, type=str  , help='tested for "M3COG" and "GBNO2"')
parser.add_argument('--ML-model'      , default=False, type=str  , help='e.g., "ML/NN/overlaps-NEWB-.../distlog090.../trained_model_M3COG"')

args = parser.parse_args()

GRO_FILE        = args.gro_file       # e.g., '../2A_sample_conformations/sampling-at-{TEMP}K_{LABEL}_oplswB97XD/run.tpr'
TPR_FILE        = args.tpr_file       # e.g., '../2A_sample_conformations/sampling-at-{TEMP}K_{LABEL}_oplswB97XD/run-snap{SNAP}-whole.gro'
TEMP            = args.temp           # e.g., '300'
LABEL           = args.label          # e.g., 'G'
SNAP            = args.snap           # e.g., '49ns'
CUTOFF          = args.cutoff         # in AA systems, it was 10 (hence the default)
RESNAME         = args.resname        # e.g., "PEOPH", "NMPH"
MAPjsonFILE     = args.map_file       # default = 'all_mappings_PTMA_39atoms.json' 
PDBtoALIGNn     = args.pdb_to_align_n # PDB file of the structure to be aligned
PDBtoALIGNr     = args.pdb_to_align_r # PDB file of the structure to be aligned
ALSO_CM         = args.also_cm        # 'True' if you want *also* Coulomb Matrices 
ALSO_1D         = args.also_1d        # 'True' if you want *also* 1D flattened matrices containing the ecli/stag identities 
VERBOSE         = args.verbose
TEST            = args.test
EVERY_NTH_FRAME = args.trajstep
OVERLAP_TYPE    = args.v_type
FEATURE         = args.feature
MAKE_PLOT       = args.make_plot
MAPPING         = args.mapping
TRAINED_ML_MODEL_PATH = args.ML_model
print(f"\nReading data from {GRO_FILE} and {TPR_FILE} ({TEMP}K, {LABEL}, {SNAP}ns); VERBOSE is set to {VERBOSE}.")

# Set up some paths and folders
REPOBASE   = os.path.realpath(os.path.join(os.path.dirname(__file__), '.')) # script now resides in the "REPOBASE"
WORKDIR    = os.getcwd()
print(f"REPOBASE is {REPOBASE}")
print(f"WORKDIR  is {WORKDIR}" )
OUTPUTDIR  = f"pair-predictions-{TEMP}K-{LABEL}-{SNAP}"
OUTPUTDATs = "pair_DATs"
PDBtoALIGNn = mda.Universe( os.path.join( WORKDIR, PDBtoALIGNn ) )
PDBtoALIGNr = mda.Universe( os.path.join( WORKDIR, PDBtoALIGNr ) )
Path(f"{OUTPUTDATs}").mkdir(parents=True, exist_ok=True)
Path(f"{OUTPUTDIR}").mkdir(parents=True, exist_ok=True)

# Set up logging to file
logfilepath = os.path.join(OUTPUTDIR,'predict_overlaps_cutoff{0:02d}A.log'.format(int(CUTOFF)))
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename=logfilepath,
                    filemode='w')
console = logging.StreamHandler()                                          # define a Handler which writes INFO messages or higher to the sys.stderr
console.setLevel(logging.INFO)                                             
formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')  # set a format which is simpler for console use
console.setFormatter(formatter)                                            # tell the handler to use this format
logging.getLogger().addHandler(console)                                    # add the handler to the root logger
mda.stop_logging()                                                         # stop MDAnalysis from filling the logging file
logging.info('Predicting electronic couplings (orbital overlaps) w/ ECG.') # we can log to the root logger, or any other logger. First the root.
logger1 = logging.getLogger('myapp.area1')                                 # define another logger (different loggers might represent areas in the code)


# Load the morphology
u = mda.Universe(TPR_FILE, GRO_FILE) 
logger1.info(f"System (resname={RESNAME}) with {len(u.atoms.fragments)} chains, {len(u.atoms)} atoms, and {len(u.trajectory)} frames.")
logger1.info(f'- NOTE! - resname = {RESNAME} --> will reduce the distance matrix from 16x16 to 12x12 (i.e., discard the Hs)')



#=================================================#
# 1) Load the trained NN                          #
#=================================================#

# Loading of the trained NN and of the *fitted* `X_scaler` and `Y_scaler`
model = keras.models.load_model( os.path.join( WORKDIR, TRAINED_ML_MODEL_PATH ) ) 
X_scaler = joblib.load( os.path.join(WORKDIR, TRAINED_ML_MODEL_PATH, 'X_scaler.joblib') )
Y_scaler = joblib.load( os.path.join(WORKDIR, TRAINED_ML_MODEL_PATH, 'Y_scaler.joblib') )

logger1.info(f"\nINFO - printing a summary of the model (might be at the bottom of the log file): {model.summary()}")
with open( logfilepath, 'a') as log: # needed to print `model.summary()` to file
    with redirect_stdout(log):
        model.summary()



#=================================================#
# 2) Go through the pairs and predict couplings   #
#=================================================#

# We use a MONOMER-MONOMER COM distance of {CUTOFF} ang as cutoff. 
# 
# Pseudo-code to obtain list of residue pairs between which to compute :
# 1. compute distance matrix betweem COMs of all MONOMERS
#     1. compute and store in an array all MONOMER COMs 
#     2. compute the distance matrix of those MONOMER COMs
# 2. for each MONOMER-MONOMER pair:
#     1. If COM-COM distance is equal or less than the {CUTOFF}:
#         1. write DAT file with info
#             - INFO: `pair_index`; `COM-COM distance`; `resid1`; `resid2`;
#         2. infer orbital overlap by using the loaded ML surrogate model


# 2A. Compute distance matrix between COMs

# Make sure chains are whole (take 5 secs for 100 chains with N=30)
for fragment in u.atoms.fragments:
    mda.lib.mdamath.make_whole(fragment)

MONOMERs = u.select_atoms(f"resname {RESNAME}")
logger1.info(f"- INFO - There are {len(MONOMERs.residues)} MONOMER radical sites.")

# read-in sites of N-methyl-phthalimide
if os.path.exists( os.path.join(WORKDIR, MAPjsonFILE) ):
    with open( os.path.join(WORKDIR, MAPjsonFILE) ) as json_mappings:
        mappings = json.load(json_mappings)

GROUP_selection = mappings["AA"]
logger1.info(f'- INFO - GROUP_selection = {GROUP_selection}')


# 2B. Iterate over all pairs and predict couplings (and save them to file) for the pairs within the cutoff 

pair_index =    0  # initialize pair_index
predicted_overlaps = []
pair_indices       = []

logger1.info(f"\nINFO - prediction - Starting with the prediction.")
if TEST:
    N_max =  2
else:
    N_max = -1

begin = time.time()

with open( os.path.join(OUTPUTDATs,'pairs_info_cutoff{0:02d}A_{1}K_{2}_{3}.dat'.format(int(CUTOFF),TEMP,LABEL,SNAP)), 'w') as pairs_info:

    for frame in u.trajectory[::EVERY_NTH_FRAME]: # Iterate over the frames
        if VERBOSE:
            print(f'frame = {frame}; time = {u.trajectory.time}')

        # Get the coordinates for the reference atoms and store them into an array with dimensions (COM_MONOMERs, 3)
        COM_MONOMER = []
        for residue in MONOMERs.residues:
            COM_MONOMER.append(residue.atoms.select_atoms(f'{GROUP_selection}').center_of_mass()) # all the atoms of N-methyl-phthalimide
        COM_MONOMER = np.row_stack(COM_MONOMER).astype('float32')
        # Distances matrix between N09s of the MONOMER
        COM_COM = mda_dist.distance_array(COM_MONOMER, COM_MONOMER, box=u.dimensions)
        logger1.info(f"- INFO - COM-COM distance matrix with dimensions {COM_COM.shape} computed.")

        # Iterate over the upper triangle of the N09-N09 distance matrix (diagonal excluded!)
        for i_MONOMER, MONOMERi_residue in enumerate(MONOMERs.residues[:N_max], start=0):
            for j_MONOMER, MONOMERj_residue in enumerate(MONOMERs.residues[i_MONOMER+1:], start=i_MONOMER+1):

                if COM_COM[i_MONOMER, j_MONOMER] <= CUTOFF:
                    
                    pair_index+=1 # pair_index starts from 1, effectively!

                    pairs_info.write("{0:10d} {1:12.8f} {2:10d} {3:10d} {4:15.3f} ".format(
                                     pair_index, COM_COM[i_MONOMER, j_MONOMER], MONOMERi_residue.resid, MONOMERj_residue.resid, u.trajectory.time)
                                     + " # pair_index  COM-COM_dist  MONOMERi_resID  MONOMERj_resID  timestamp_in_ps\n"
                                    )

                    MONOMERi_residue_aligned, MONOMERi_ref_residue = srcfunctions.align_MONOMER(u, MONOMERi_residue, RESNAME, "neutral"      , GROUP_selection, PDBtoALIGNn, VERBOSE)
                                                                                                                                                               
                    MONOMERj_residue_aligned, MONOMERj_ref_residue = srcfunctions.align_MONOMER(u, MONOMERj_residue, RESNAME, "radical_anion", GROUP_selection, PDBtoALIGNr, VERBOSE)

                    # Retrieve the distance matrix of this pair
                    reciprocal_distmatrix = srcfunctions.fast_compute_and_store_reciprocal_distmat(
                                                         "pair",
                                                         pair_index, 
                                                         MONOMERi_residue_aligned, 
                                                         MONOMERj_residue_aligned, 
                                                         MAPjsonFILE,
                                                         0, 0,  # --> only needed for distinguishing multiple conformations in the 1D-matrix case (see function docs) 
                                                         False, # --> do NOT save matrices to file
                                                         OUTPUTDIR, 
                                                         ALSO_CM, ALSO_1D, 
                                                         u, VERBOSE)

                    # AD-HOC FIX for PMAP/PEPP/PVBP
                    if RESNAME in ["PMAP","PEPP","PVBP"] and reciprocal_distmatrix.shape[0] == 16:
                        # Reduce matrix from 16x16 to 12x12 (i.e., I want to disregard the hydrogens, if present)
                        reciprocal_distmatrix = reciprocal_distmatrix[:12,:12]

                    # Pre-process the distance matrix
                    reciprocal_distmatrix = reciprocal_distmatrix.flatten()     # --> (36,0)
                    reciprocal_distmatrix = reciprocal_distmatrix.reshape(-1,1) # --> (36,1)
                    reciprocal_distmatrix = np.transpose(reciprocal_distmatrix) # --> (1,36)

                    # Use loaded model to make prediction using the reciprocal_distmatrix
                    y_predicted = model.predict( X_scaler.transform( reciprocal_distmatrix) ) 
                    y_predicted = Y_scaler.inverse_transform( y_predicted )
                    predicted_overlaps.append(y_predicted[0,0]) 
                    pair_indices.append(f'{pair_index:06}')

end = time.time()
logger1.info(f"\nINFO - prediction - Prediction completed in {round(end-begin,3)} seconds.")

logger1.info(f"\n**DONE** {pair_index} COM-COM distances are within the CUTOFF.")

with open(os.path.join(OUTPUTDIR, "overlaps_predicted_cutoff{0:02d}A.csv".format(int(CUTOFF))), 'w') as f:
    writer = csv.writer(f)
    writer.writerows(zip(pair_indices, predicted_overlaps))

logger1.info("**DONE** predicted overlaps written to 'overlaps_predicted_cutoff{0:02d}A.csv'.\n".format(int(CUTOFF)))


