#!/usr/bin/env python3

import os, sys
from MDAnalysis.analysis import align
from MDAnalysis.analysis import distances as mda_dist
import numpy as np
import json

# Set up some paths and folders
REPOBASE   = os.path.realpath(os.path.join(os.path.dirname(__file__), '..')) # define path of "REPOBASE"
WORKDIR    = os.getcwd()
LIBDIR     = os.path.join( REPOBASE, "lib")
print(f"@RA@ REPOBASE={REPOBASE}")
print(f"@RA@ WORKDIR={WORKDIR}")
print(f"@RA@ LIBDIR={LIBDIR}")
 

def align_MONOMER(universe, residue, RESNAME, neutral_or_cation,  selection_for_alignment, pdb_to_align, VERBOSE=False):
    """
    ...

    The `align.alignto` uses the fast QCP algorithm to calculate the root mean square distance (RMSD) 
    between two coordinate sets [Theobald2005] and the rotation matrix R that minimizes the RMSD [Liu2010].
    
    [Theobald2005] Douglas L. Theobald (2005), Acta Crystallographica A 61(4):478-480.
    [Liu2010]      Pu Liu, Dmitris K. Agrafiotis, and Douglas L. Theobald (2010), J. Comput. Chem. 31, 1561-1563.

    NOTE: It takes only 36 seconds for a morphology with 100 chains with $N=30$.
    """
    
    QC_opt = pdb_to_align

    if VERBOSE:
        print("Loading the QC-optimized *{0}* conformation; contains {1} atoms "
              .format( neutral_or_cation, len(QC_opt.atoms) )) 

    mobile = QC_opt.atoms # reset the mobile structure to the QC-optimized one (may not be needed)

    if residue.resname == RESNAME: # skip the termini, which are called something else 

        ref_MD_frame = universe.select_atoms(f"resid {residue.resid}")
        
        align_rmsd = align.alignto(mobile, ref_MD_frame, 
                                   select=f"{selection_for_alignment}",
                                   weights="mass")
        
        if (align_rmsd[0] - align_rmsd[1]) <= -0.1 :
            sys.exit(f"ERROR!? RMSD before ({round(align_rmsd[0],4)}) < after ({round(align_rmsd[1],4)}) alignment!? Please check what's up.")
        else:
            if VERBOSE:
                print(f"- INFO - RMSD before ({round(align_rmsd[0],4)}) and after ({round(align_rmsd[1],4)}) alignment.")
    
    else:
        sys.exit("ERROR! Residue is not *PTMA*? Exiting...")
    
    return mobile, ref_MD_frame


def fast_compute_and_store_reciprocal_distmat(mon_or_pair, index, MOLi_residue, MOLj_residue, MAPjsonFILE, i_ID, j_ID, 
                                              STORE, OUTDIR, ALSO_CM, ALSO_1D, universe, VERBOSE):
    """
    Computes and (optionally) stores to `.npy` files a reciprocal distance matrix for each mapping resolution.
    
    Parameters
    ----------
    mon_or_pair: string
        Either "mon" or "pair" are available; based on that, 
        differet labels and mappings are used.
    index: int
        Index of the monomer/pair.
    MOLi_residue: MDAnalysis residue
        Residue of molecule I in the IJ pair.
    MOLj_residue: MDAnalysis residue
        Residue of molecule J in the IJ pair.
    MAPjsonFILE: string
        Name of .json file containing the mapping(s) for which the reciprocal distance matrix/matrices should be computed.
    i_ID: int
        If i_ID == 2, MOL_i is a staggered conformation; if i_ID == 1, MOL_i is an eclipsed conformation
        Only relevant to the 1D featurization.
    j_ID: int
        If j_ID == 2, MOL_i is a staggered conformation; if j_ID == 1, MOL_i is an eclipsed conformation
        Only relevant to the 1D featurization.
    STORE: bool
        If True, stores the computed matrices to file.
    OUTDIR: string
        Name of output directory. 
    ALSO_CM: bool
        If True, produce also a Coulomb matrix. 
    ALSO_1D: bool
        If True, produce also 1D flattened versions of the matrices. 
    universe: MDAnalysis.universe
        MDAnalysis universe.
    VERBOSE: bool
        If True, prints more information. 

    Returns
    --------
    reciprocal_distmatrix: ndarray
        The reciprocal distance matrix. 
    """
    
    if mon_or_pair in ["pair", "mon"]:
        
        # Read mappings from file
        if os.path.exists( os.path.join(WORKDIR, MAPjsonFILE) ):
            with open( os.path.join(WORKDIR, MAPjsonFILE) ) as json_mappings:
                mappings = json.load(json_mappings)
        else:
            with open( os.path.join(LIBDIR, MAPjsonFILE) ) as json_mappings:
                mappings = json.load(json_mappings)
    else:
        sys.exit("I know what to do only with the keywords 'mon' or 'pair'. Check that.")
    
    for resolution in mappings:
        
        if "GB" in resolution or "COG" in resolution:
            
            positions_i = np.zeros(( len(mappings[resolution].split(',')) ,3))
            positions_j = np.zeros(( len(mappings[resolution].split(',')) ,3))
   
            for bead_idx, atoms in enumerate(mappings[resolution].split(','), start=0):
                if not len(MOLi_residue.atoms.select_atoms(f"name {atoms}")) == 0:
                    positions_i[bead_idx] = MOLi_residue.atoms.select_atoms(f"name {atoms}").center_of_geometry()
                else:
                    sys.exit(f'MOLi_residue.atoms.select_atoms(f"name {atoms}") is EMPTY. Exiting...')
                if not len(MOLj_residue.atoms.select_atoms(f"name {atoms}")) == 0:
                    positions_j[bead_idx] = MOLj_residue.atoms.select_atoms(f"name {atoms}").center_of_geometry()
                else:
                    sys.exit(f'MOLj_residue.atoms.select_atoms(f"name {atoms}") is EMPTY. Exiting...')
   
            dist_matrix = mda_dist.distance_array(positions_i,
                                                  positions_j,
                                                  box=universe.dimensions)

            # Take the reciprocal of it
            reciprocal_distmatrix = 1. / dist_matrix
            
            # Fix the diagonal by replacing 'inf' with '0' (necessary for the monomer case)
            reciprocal_distmatrix[reciprocal_distmatrix == np.inf] = 0

            if STORE:
                np.save( os.path.join(OUTDIR,'{0}{1:06d}{2}.npy'.format(mon_or_pair, index, resolution)),
                        reciprocal_distmatrix)
   
            if VERBOSE:
                print(f"- INFO - residues with # atoms: {len(positions_i)} and {len(positions_j)}")
             
        else:
            MOLi_residue_mapped = MOLi_residue.atoms.select_atoms(f"{mappings[resolution]}")
            MOLj_residue_mapped = MOLj_residue.atoms.select_atoms(f"{mappings[resolution]}")

            reciprocal_distmatrix = compute_recip_distmat(MOLi_residue_mapped, # MOLi
                                                          MOLj_residue_mapped, # MOLj
                                                          universe, VERBOSE)

            if STORE:
                np.save( os.path.join(OUTDIR,'{0}{1:06d}{2}.npy'.format(mon_or_pair, index, resolution)),
                        reciprocal_distmatrix)
    
            if ALSO_CM and STORE:
                coulomb_matrix = compute_coulmat_inter(MOLi_residue_mapped, # MOLi
                                                       MOLj_residue_mapped, # MOLj
                                                       resolution,
                                                       universe, VERBOSE)
     
                np.save( os.path.join(OUTDIR,'{0}{1:06d}{2}_CM.npy'.format(mon_or_pair, index, resolution)),
                        coulomb_matrix)
         
        if ALSO_1D and STORE:
            # Save arrays with the stag/ecli info
            reciprocal_distmatrix_1D = np.append(reciprocal_distmatrix.flatten(), [i_ID, j_ID])
            np.save( os.path.join(OUTDIR,'{0}{1:06d}{2}_1D.npy'.format(mon_or_pair, index, resolution)),
                    reciprocal_distmatrix_1D)

            if ALSO_CM and STORE:
                coulomb_matrix_1D        = np.append(coulomb_matrix.flatten(),        [i_ID, j_ID])
                np.save( os.path.join(OUTDIR,'{0}{1:06d}{2}_CM_1D.npy'.format(mon_or_pair, index, resolution)),
                        coulomb_matrix_1D)

        if VERBOSE:
            print(f"- INFO - reciprocal distance matrix with dimensions {reciprocal_distmatrix.shape} computed.")
            print(f"{reciprocal_distmatrix}\n")
        if VERBOSE and ALSO_CM:
            print(f"- INFO - Coulomb Matrix *inter* with dimensions {coulomb_matrix.shape} computed.")
            print(f"{coulomb_matrix}\n")

    return reciprocal_distmatrix


def compute_recip_distmat(MOLi_residue, MOLj_residue, universe, VERBOSE=False):
    """
    Computes the *reciprocal* distance Matrix between residues *i* and *j.
    """    
    # Compute distance matrix
    dist_matrix = mda_dist.distance_array(MOLi_residue.atoms.positions, # MOLi
                                          MOLj_residue.atoms.positions, # MOLj
                                          box=universe.dimensions)
    # Take the reciprocal of it
    reciprocal_distmatrix = 1. / dist_matrix

    # Fix the diagonal by replacing 'inf' with '0' (necessary for the monomer case)
    reciprocal_distmatrix[reciprocal_distmatrix == np.inf] = 0

    if VERBOSE:
        print(f"- INFO - residues with # atoms: {len(MOLi_residue.atoms)} and {len(MOLj_residue.atoms)}")
    
    return reciprocal_distmatrix


def compute_coulmat_inter(MOLi_residue, MOLj_residue, resolution, universe, VERBOSE=False):
    """
    Computes the *intermolecular* Coulomb Matrix between residues *i* and *j.
    According to [Rupp, et al., PRL 108, 058301 (2012)], the i-j terms are defined as:

        C_ij = Z_i*Z_j / ||R_i - R_j||

    Note that we do NOT have the C_ii = 0.5 * Z_i**2.4 terms.
    """
    # Compute distance matrix to be used later
    dist_matrix = mda_dist.distance_array(MOLi_residue.atoms.positions, # MOLi
                                          MOLj_residue.atoms.positions, # MOLj
                                          box=universe.dimensions)

    # Source the charges for the given resolution
    charges = np.load(os.path.join('lib',f'charges_MOL_{resolution}.npy'))

    # Initialize Coulomb Matrix "inter"
    coulomb_matrix = np.zeros((len(MOLi_residue.atoms),len(MOLj_residue.atoms)))

    # Coulomb Matrix "inter"
    for idx_i, atom_i in enumerate(MOLi_residue.atoms, start=0):
        for idx_j, atom_j in enumerate(MOLj_residue.atoms, start=0):
            # Z_i*Z_j / ||R_i - R_j||
            coulomb_matrix[idx_i, idx_j] = charges[idx_i]*charges[idx_j] / dist_matrix[idx_i, idx_j]

    if VERBOSE:
        print(f"- INFO - residues with # atoms: {len(MOLi_residue.atoms)} and {len(MOLj_residue.atoms)}")

    return coulomb_matrix


