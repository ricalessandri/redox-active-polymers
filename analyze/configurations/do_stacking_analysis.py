#!/usr/bin/env python3
# coding: utf-8
"""
USAGE:

    python do_stacking_analysis.py -n PMAPH -a PMAPH.json
"""


# TODO:
# > change this hardcoded lines (#@RA@#) if wanting to generalize beyond phthalimide !!! 
# > merge with the other `run/do_stacking_analysis.py`? 
# > "Error - atom mismatch" - either implement something smarter or perhaps not needed?

import MDAnalysis as mda
from MDAnalysis.analysis import distances as mda_dist
import numpy as np
import matplotlib.pyplot as plt
import time
import argparse
import sys
import json
import os
from pathlib import Path

PDBs=True
PDBs=False

# Parse the arguments
parser = argparse.ArgumentParser(description='Run with:\n\n  python run_stacking_analysis.py -n PMAPH')
parser.add_argument('-n', '--name'     , required=True, type=str, help='name of the molecule')
parser.add_argument('-a', '--atomnames', required=True, type=str, help='name of JSON file containing atomnames of molecular unit (e.g., phthalimide)')
parser.add_argument('-c', '--cutoff'   , required=True, type=int, help='cutoff for the dimers') 
parser.add_argument('-f', '--gro-file' , required=True, type=str, help='name of GRO file to analyze')
parser.add_argument('-s', '--tpr-file' , required=True, type=str, help='name of TPR file to analyze')
parser.add_argument('-l', '--label'    , required=True, type=str, help='label')
parser.add_argument('--test'           , default=False,             help='True if you want to run a test')

args = parser.parse_args()

NAME          = args.name # e.g., 'PMAPH' 
jsonATOMNAMES = args.atomnames
CUTOFF        = int(args.cutoff)  # e.g., 8 A
GRO_FILE      = args.gro_file
TPR_FILE      = args.tpr_file
WORKDIR       = os.getcwd()
LABEL         = args.label
TEST          = args.test


INPUTTPR= os.path.join(WORKDIR, TPR_FILE) 
INPUTGRO= os.path.join(WORKDIR, GRO_FILE) 
print(f"Reading files '{TPR_FILE}' and '{GRO_FILE}' ...")
u = mda.Universe(INPUTTPR, INPUTGRO)

print(f"Settings: -n {NAME}; -a {jsonATOMNAMES}")
print(f"System with {len(u.atoms.fragments)} chains, {len(u.atoms)} atoms, and {len(u.trajectory)} frames.")
print(f"System with {len(u.atoms.fragments)} chains, {len(u.atoms)} atoms, and {len(u.residues)} residues.")

MONOMERs = u.select_atoms(f"resname {NAME}")
print(f"- INFO - There are {len(MONOMERs.residues)} MONOMER radical sites.")


# Let's produce a `{NAME}_dimers_cutoff10ang.dat` file that contains:
#
# distCOMiCOMj  angle_between_plane_normals  dist_NiNj  dih_NiCOMiCOMjNj  resID_i  resID_j # r_COMiCOMj ang_planes r_NN dih_NiCOMiCOMjNj resID_i resid_j

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    """ Returns the angle in degrees between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            90.0
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            180.0
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.degrees( np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)) )


def normal_to_plane(p1,p2,p3):
    """
    Return the vector normal to the plane defined by p1, p2, and p3.
    """
    # These two vectors are in the plane
    vec1 = p2 - p1
    vec2 = p3 - p1
    return np.cross(vec1, vec2)


def angle_between_plane_normals(PHTH_i, PHTH_j):
    """ 
    Returns the angle between the vector normal to two molecular units.
    """
    ## if PHTH_i.atoms[0].name != 'N08' or PHTH_i.atoms[4].name != 'C13' or PHTH_i.atoms[5].name != 'C14' or \
    ##    PHTH_j.atoms[0].name != 'N08' or PHTH_j.atoms[4].name != 'C13' or PHTH_j.atoms[5].name != 'C14' :
    ##     sys.exit('Error - atom mismatch: should be N08, C13, C14 (x2) but they are: {0} {1} {2} and {3} {4} {5}'
    ##              .format(PHTH_i.atoms[0].name, PHTH_i.atoms[4].name, PHTH_i.atoms[5].name,
    ##                      PHTH_j.atoms[0].name, PHTH_j.atoms[4].name, PHTH_j.atoms[5].name) )
        
    N_i  = PHTH_i.atoms[1].position # #@RA@# hard-coded for phthalamide !! ; was "0"
    Ca_i = PHTH_i.atoms[6].position # #@RA@# hard-coded for phthalamide !! ; was "4"
    Cb_i = PHTH_i.atoms[7].position # #@RA@# hard-coded for phthalamide !! ; was "5"
    N_j  = PHTH_j.atoms[1].position # #@RA@# hard-coded for phthalamide !! ; was "0"
    Ca_j = PHTH_j.atoms[6].position # #@RA@# hard-coded for phthalamide !! ; was "4"
    Cb_j = PHTH_j.atoms[7].position # #@RA@# hard-coded for phthalamide !! ; was "5"
    
    # Compute normal vectors to PHTH_i & PH_TH_j planes
    PHTH_i_normal = normal_to_plane(N_i, Ca_i, Cb_i)
    PHTH_j_normal = normal_to_plane(N_j, Ca_j, Cb_j)
    
    return angle_between(PHTH_i_normal,PHTH_j_normal)


# Read in atomnames of molecular unit for the polymer {NAME}
if os.path.exists( os.path.join(WORKDIR, jsonATOMNAMES) ):
    with open( os.path.join(WORKDIR, jsonATOMNAMES) ) as json_atomnames:
        atomnames = json.load(json_atomnames)
else:
    sys.exit(f"File {jsonATOMNAMES} not found.")

NMPHTH_selection = atomnames["AA"]
print(f'- INFO - NMPHTH_selection = {NMPHTH_selection}') 

if TEST:
    N_max = 10
else:
    N_max = -1   # all of them
pair_ID = 0  # initialize pair count

tic = time.time()
with open( os.path.join(WORKDIR, '{0}_dimers_{1}_cutoff{2:02d}A.dat'.format(NAME, LABEL, CUTOFF)), 'w') as outdat:

    for frame in u.trajectory: # Iterate over the frames
        ## print(f'frame = {frame}; time = {u.trajectory.time}') ## VERBOSE

        ## for idx_i, residue_i in MONOMERs.residues[:]:
        ##     for residue_j in MONOMERs.residues[residue_i.resid:]: # do only upper triangle *and* skip the diagonal
        for i_MONOMER, MONOMERi_residue in enumerate(MONOMERs.residues[:N_max], start=0):
            for j_MONOMER, MONOMERj_residue in enumerate(MONOMERs.residues[i_MONOMER+1:], start=i_MONOMER+1): # do only upper triangle *and* skip the diagonal

                # COG-COG distance
                plane_i = MONOMERi_residue.atoms.select_atoms(f'{atomnames["AA"]}')
                plane_j = MONOMERj_residue.atoms.select_atoms(f'{atomnames["AA"]}')

                dist_COM = mda_dist.distance_array(plane_i.atoms.center_of_mass(), plane_j.atoms.center_of_mass(), box=u.dimensions)[0][0] # use the minimum image convention!

                if dist_COM < CUTOFF:
                    
                    pair_ID = pair_ID + 1

                    # N-N distance
                    N_i = MONOMERi_residue.atoms.select_atoms(f'name {atomnames["AA"].split()[2]}')
                    N_j = MONOMERj_residue.atoms.select_atoms(f'name {atomnames["AA"].split()[2]}')
                    dist_NN = mda_dist.distance_array(N_i.atoms[0].position, N_j.atoms[0].position, box=u.dimensions)[0][0] # use the minimum image convention!

                    # angle between plane normals 
                    angle = angle_between_plane_normals(plane_i,plane_j)

                    # dih
                    dihedral = np.rad2deg( mda.lib.distances.calc_dihedrals(N_i.positions[0], plane_i.atoms.center_of_mass(),
                                                                            plane_j.atoms.center_of_mass(), N_j.positions[0],
                                                                            box=u.dimensions) )

                    # Check if intra or inter
                    pair = MONOMERi_residue + MONOMERj_residue
                    if len(pair.atoms.fragments) > 1:
                        pair_type = 1 # 'inter'
                    else:
                        pair_type = 0 # 'intra'

                    # MAKE PDBs
                    if PDBs:
                        plane_i.write('plane_i_pair_{0:06d}_resid{1:06d}.pdb'.format(pair_ID, MONOMERi_residue.resid) )
                        plane_j.write('plane_j_pair_{0:06d}_resid{1:06d}.pdb'.format(pair_ID, MONOMERj_residue.resid) )
                        
                    outdat.write('{0:8d} {1:12.6f} {2:12.6f} {3:12.6f} {4:12.6f} {5} {6:6d} {7:6d} # pair_ID dist_COM ang dist_NN dih 0(=intra)/1(=inter) resID_i resID_j \n'
                                .format(pair_ID, dist_COM, angle, dist_NN, dihedral, pair_type, MONOMERi_residue.resid, MONOMERj_residue.resid) )

toc = time.time()
print(f"It took {round(toc-tic,3)} seconds for scanning through {N_max*(N_max-1)/2} resid dimers and {len(u.trajectory)} frames.")

