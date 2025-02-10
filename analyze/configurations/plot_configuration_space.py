#!/usr/bin/env python3
# coding: utf-8
"""

"""

import MDAnalysis as mda
from MDAnalysis import selections
from MDAnalysis.analysis import distances as mda_dist
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import math
import argparse
import sys
import os
import pandas as pd
import seaborn as sns


plt.style.use('seaborn-talk')
plt.rcParams['font.family'] = 'sans'
plt.rcParams['font.size'] = 18
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.labelweight'] = 'normal'
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['legend.fontsize'] = 18
plt.rcParams['figure.titlesize'] = 18 



# Parse the arguments
parser = argparse.ArgumentParser(description='Run with:\n\n  python run_stacking_analysis.py -n PMAPH')
parser.add_argument('-n', '--name'     , required=True, type=str, help='name of the molecule')
parser.add_argument('-f', '--folder'   , required=True, type=str, help='name of the folder')
parser.add_argument('-c', '--cutoff'   , required=True, type=int, help='cutoff for the dimers (ang)') 
parser.add_argument('-l', '--label'    , required=True, type=str, help='label')
parser.add_argument('-r', '--resname'  , default=None,  type=str, help='RESNAME, if different from MOLNAME (-name)')
args = parser.parse_args()
MOL         = args.name        # e.g., 'PMMPH' 
FOLDER      = args.folder      # e.g., PMMPH000charge_DME_TBAPF6_05percent
CUTOFF      = int(args.cutoff) # e.g., 10 (ang)
LBL         = args.label       # e.g. "30mer300KD100ns"
if args.resname == None:
    RESNAME = MOL
else:
    RESNAME = args.resname

# Setup
LABEL  = 'dimers_{0}_cutoff{1:02d}A'.format(LBL, CUTOFF)
WORKDIR = os.getcwd()
FOLDER  = os.path.join( WORKDIR, FOLDER ) 


# ### 2. EXPLORING THE CONFORMATIONAL LANDSCAPE
# We load the coordinates of the donor-acceptor (DA) complexes generated with step 1.
# We then prepare the *aggregate* arrays, *i.e.*, arrays which contain the $r_{DA}$ and $angle_{DA}$ coordinates 
# for the 1st, 2nd, and 3rd closest thiophene to the *i-th* fullerene.

# Load the coordinates file
file_cutoff      = [ (os.path.join( FOLDER, '{0}_{1}.dat'.format(RESNAME, LABEL)), 'cutoff') ]
datalist_cutoff      = [ ( np.loadtxt(filename), label ) for filename, label in file_cutoff ]

#-------------------------------#
#    CUTOFF-BASED ARRAYS        #
#-------------------------------#
# Retrieve the r_DA and angle_DA based on the cutoff.
x_COM      = np.column_stack([data[:, 1] for data, label in datalist_cutoff]).flatten() # COM-COM distance
ang_PLANES = np.column_stack([data[:, 2] for data, label in datalist_cutoff]).flatten()
x_NN       = np.column_stack([data[:, 3] for data, label in datalist_cutoff]).flatten() # N-N distance
dih_NCoCoN = np.column_stack([data[:, 4] for data, label in datalist_cutoff]).flatten() # DIHEDRAL N-COM-COM-N
intra_inter = np.column_stack([data[:, 5] for data, label in datalist_cutoff]).flatten() # intra=0; inter=1 

intra_count = np.count_nonzero(intra_inter==0)
inter_count = np.count_nonzero(intra_inter==1)

print(f'% of intra-molecular pairs = {round(intra_count/len(intra_inter)*100,2)}')
print(f'% of intra-molecular pairs = {round(inter_count/len(intra_inter)*100,2)}')
with open( os.path.join( FOLDER, f'intrainter_stats_{LABEL}.txt'), 'w') as stats:
    stats.write(f'Total number of dimers     = {len(intra_inter)}\n')
    stats.write(f'% of intra-molecular pairs = {round(intra_count/len(intra_inter)*100,1)}\n')
    stats.write(f'% of inter-molecular pairs = {round(inter_count/len(intra_inter)*100,1)}\n')

## print(riccardo)

print("\n Based on the cutoff ({0} ang), there are {1} DA pairs -- r_DA = {1}; angle_DA = {2}"
      .format(CUTOFF, len(x_COM), len(ang_PLANES)))

print( x_COM.min()    , x_COM.max()     )
print( ang_PLANES.min(), ang_PLANES.max() )


# 2a. Make cutoff plot:
def make_cutoff_plot(x_COM, ang_PLANES, x_min, x_max, y_min, y_max, xbin_width, ybin_width):
    """
    """
    # modify x_min and x_max so to have centered bins already
    x_min      = x_min - xbin_width/2
    x_max      = x_max + xbin_width/2
    y_min      = y_min - ybin_width/2
    y_max      = y_max + ybin_width/2
    #
    no_xbins = int((x_max-x_min) / xbin_width)
    no_ybins = int((y_max-y_min) / ybin_width)
    
    print("Bin widths for x and y :: {0} ang; {1} deg".format(xbin_width,ybin_width))
    print("Number of x and y bins :: {0}     ; {1}    ".format(no_xbins, no_ybins))
    
    # Do the binning.
    H, xedges, yedges = np.histogram2d(x_COM, ang_PLANES, 
                                       range=[[x_min, x_max], [y_min,y_max]], 
                                       bins=[no_xbins, no_ybins])
    
    print("xedges, yedges and histogram shapes:")
    print(xedges.shape, yedges.shape, H.shape)
    
    # contours are *point* based plots, so convert our bound into point centers
    X_for_contourf = xedges[:-1]+xbin_width/2
    Y_for_contourf = yedges[:-1]+ybin_width/2
    print("xedges, yedges  vs  X_for_contourf, Y_for_contourf: ")
    print(xedges.shape, yedges.shape, " vs ", X_for_contourf.shape, Y_for_contourf.shape)
    
    # Prepare X and Y
    X, Y = np.meshgrid(X_for_contourf,Y_for_contourf)
    
    print("\nThe following 3 must have the same shape:")
    print(X.shape, Y.shape, H.T.shape)

    # **Note**: `xedges` and `yedges` are basically the following:
    # ``` 
    # xedges = np.arange(x_min, x_max+0.01, xbin_width)
    # yedges = np.arange(y_min, y_max+0.01, ybin_width)
    # ```
    
    # The radial distribution function is usually determined by calculating the distance between 
    # all particle pairs and binning them into a histogram. The histogram is then normalized with respect to an ideal gas, 
    # where particle histograms are completely uncorrelated. For three dimensions, this normalization is the number density 
    # of the system ($\rho$) multiplied by the volume of the spherical shell, which symbolically can be expressed as $\rho \,4\pi r^{2}dr$.
    
    density = 22.5*0.001 # ang^{-3} = 22.5 # nm^{-3}
    density = 1 # Normalize only for the volume
    
    def volume_normalization_factor(density, r, dr):
        return density * 4 * np.pi * r**2 * dr
    
    HT_volu = np.zeros(H.T.shape)
    HT_norm = np.zeros(H.T.shape)
    
    for index, r in enumerate(X_for_contourf): 
        HT_volu[:,index] = H.T[:,index] / volume_normalization_factor(density,r,xbin_width)
    
    # Normalize by dividing by the maximum of `HT_norm`
    HT_norm = HT_volu / np.max(HT_volu)
    # HT_norm = HT_volu / (22.5*0.001)

    return X, Y, H, HT_volu, HT_norm, X_for_contourf, Y_for_contourf


def mirror_angles(angles):
    """
    """
    new_ang = []
    for ang in angles:
        if ang > 90:
            alpha = ang - 90
            ang_prime = 90 - alpha
            new_ang.append(ang_prime)
        else:
            new_ang.append(ang)
    return np.array(new_ang) 

def mirror_dihedrals(dihedrals):
    """
    """
    return np.abs(dihedrals)

ang_PLANES_mirrored = mirror_angles(ang_PLANES)
dih_NCoCoN_mirrored = mirror_dihedrals(dih_NCoCoN)

print(f"INFO - ang_PLANES_mirrored - {np.max(ang_PLANES_mirrored)}; {np.min(ang_PLANES_mirrored)}") 
print(f"INFO - ang_PLANES          - {np.max(ang_PLANES)}; {np.min(ang_PLANES)}") 

X, Y, H, HT_volu, HT_norm, X_for_contourf, Y_for_contourf = make_cutoff_plot(x_COM         , ang_PLANES         ,
                                                                             x_min=3.0     , x_max=float(CUTOFF),
                                                                             y_min=0.0     , y_max=180.0        ,
                                                                             xbin_width=0.5, ybin_width=10.0    )

X3, Y3, H3, HT_volu3, HT_norm3, X_for_contourf3, Y_for_contourf3 = make_cutoff_plot(x_COM         , ang_PLANES_mirrored,
                                                                                    x_min=3.0     , x_max=float(CUTOFF),
                                                                                    y_min=0.0     , y_max= 90.0        ,
                                                                                    xbin_width=0.5, ybin_width=10.0    )

## X2, Y2, H2, HT_volu2, HT_norm2, X_for_contourf2, Y_for_contourf2 = make_cutoff_plot(x_NN          , ang_PLANES         ,
##                                                                                     x_min=3.0     , x_max=float(CUTOFF),
##                                                                                     y_min=0.0     , y_max=180.0        ,
##                                                                                     xbin_width=0.5, ybin_width=10.0    )

X2, Y2, H2, HT_volu2, HT_norm2, X_for_contourf2, Y_for_contourf2 = make_cutoff_plot(x_COM         , dih_NCoCoN         ,
                                                                                    x_min=3.0     , x_max=float(CUTOFF),
                                                                                    y_min=-180.0  , y_max=180.0        ,
                                                                                    xbin_width=0.5, ybin_width=20.0    )

X4, Y4, H4, HT_volu4, HT_norm4, X_for_contourf4, Y_for_contourf4 = make_cutoff_plot(x_COM         , dih_NCoCoN_mirrored,
                                                                                    x_min=3.0     , x_max=float(CUTOFF),
                                                                                    y_min=0.0     , y_max=180.0        ,
                                                                                    xbin_width=0.5, ybin_width=20.0)

#-------------------------------#
#           PLOTS               #
#-------------------------------#
# Plot 1 - 
fig, ax = plt.subplots(2,2,figsize=(5*2,4.5*2))

ax[0,0].set_xlabel('r$_{ij}$ ($\AA$)')
ax[0,0].set_ylabel(r'$\Theta_{ij}$ (deg)')
cf1 = ax[0,0].contourf(X, Y, HT_norm,cmap ="Blues", levels = np.arange(0,1.1,0.1), antialiased=True)#, levels=levels)
ax[0,0].set_xticks(np.arange(X_for_contourf[0],X_for_contourf[-1]+1, step=1))
ax[0,0].set_yticks(np.arange(Y_for_contourf[0],Y_for_contourf[-1]+1, step=30)) # TODO: remove this hardcoded "20"
ax[0,0].set_xlim(X_for_contourf[0],X_for_contourf[-1])
ax[0,0].set_ylim(Y_for_contourf[0],Y_for_contourf[-1])

fig.colorbar(cf1, ax=ax[0,0])

ax[0,1].set_xlabel('r$_{ij}$ ($\AA$)')
ax[0,1].set_ylabel(r'$\phi_{N_i-COG_i-COG_j-N_j}$ (deg)')
cf2 = ax[0,1].contourf(X2, Y2, HT_norm2,cmap ="Blues", levels = np.arange(0,1.1,0.1), antialiased=True)#, levels=levels)
ax[0,1].set_xticks(np.arange(X_for_contourf2[0],X_for_contourf2[-1]+1, step=1))
ax[0,1].set_yticks(np.arange(Y_for_contourf2[0],Y_for_contourf2[-1]+1, step=60)) # TODO: remove this hardcoded "40"
ax[0,1].set_xlim(X_for_contourf2[0],X_for_contourf2[-1])
ax[0,1].set_ylim(Y_for_contourf2[0],Y_for_contourf2[-1])

fig.colorbar(cf2, ax=ax[0,1])

ax[1,0].set_xlabel('r$_{ij}$ ($\AA$)')
ax[1,0].set_ylabel(r'$\Theta_{ij}$ (deg)')
cf3 = ax[1,0].contourf(X3, Y3, HT_norm3,cmap ="Blues", levels = np.arange(0,1.1,0.1), antialiased=True)#, levels=levels)
ax[1,0].set_xticks(np.arange(X_for_contourf3[0],X_for_contourf3[-1]+1, step=1))
ax[1,0].set_yticks(np.arange(Y_for_contourf3[0],Y_for_contourf3[-1]+1, step=15)) # TODO: remove this hardcoded "20"
ax[1,0].set_xlim(X_for_contourf3[0],X_for_contourf3[-1])
ax[1,0].set_ylim(Y_for_contourf3[0],Y_for_contourf3[-1])

fig.colorbar(cf3, ax=ax[1,0])

ax[1,1].set_xlabel('r$_{ij}$ ($\AA$)')
ax[1,1].set_ylabel(r'$\phi_{N_i-COG_i-COG_j-N_j}$ (deg)')
cf4 = ax[1,1].contourf(X4, Y4, HT_norm4,cmap ="Blues", levels = np.arange(0,1.1,0.1), antialiased=True)#, levels=levels)
ax[1,1].set_xticks(np.arange(X_for_contourf4[0],X_for_contourf4[-1]+1, step=1))
ax[1,1].set_yticks(np.arange(Y_for_contourf4[0],Y_for_contourf4[-1]+1, step=30)) # TODO: remove this hardcoded "40"
ax[1,1].set_xlim(X_for_contourf4[0],X_for_contourf4[-1])
ax[1,1].set_ylim(Y_for_contourf4[0],Y_for_contourf4[-1])

fig.colorbar(cf4, ax=ax[1,1])



plt.tight_layout()

plt.savefig( os.path.join( FOLDER, f'conf_stack_space_{LABEL}.pdf'), bbox_inches='tight')


# Plot #2 (separate ang_PLANES (mirrored))
plt.figure(figsize=(6*0.8,5*0.8))
plt.xlabel('r$_{ij}$ ($\AA$)')
plt.ylabel(r'$\Theta_{ij}$ (deg)')
cf = plt.contourf(X3, Y3, HT_norm3,cmap ="Blues", levels = np.arange(0,1.1,0.1), antialiased=True)#, levels=levels)
plt.colorbar(cf)
plt.xticks(np.arange(X_for_contourf3[0],X_for_contourf3[-1]+1, step=1))
plt.yticks(np.arange(Y_for_contourf3[0],Y_for_contourf3[-1]+1, step=15)) # TODO: remove this hardcoded "20"
plt.xlim(X_for_contourf3[0],X_for_contourf3[-1])
plt.ylim(Y_for_contourf3[0],Y_for_contourf3[-1])
plt.tight_layout()
plt.savefig( os.path.join( FOLDER, f'conf_stack_space_{LABEL}_angPLANES.pdf'), bbox_inches='tight')


# Plot #3 (separate dih_NCoCoN (mirrored))
plt.figure(figsize=(6*0.8,5*0.8))
plt.xlabel('r$_{ij}$ ($\AA$)')
plt.ylabel(r'$\phi_{N-COG-COG-N}$ (deg)')
cf = plt.contourf(X4, Y4, HT_norm4,cmap ="Blues", levels = np.arange(0,1.1,0.1), antialiased=True)#, levels=levels)
plt.colorbar(cf)
plt.xticks(np.arange(X_for_contourf4[0],X_for_contourf4[-1]+1, step=1))
plt.yticks(np.arange(Y_for_contourf4[0],Y_for_contourf4[-1]+1, step=30)) # TODO: remove this hardcoded "20"
plt.xlim(X_for_contourf4[0],X_for_contourf4[-1])
plt.ylim(Y_for_contourf4[0],Y_for_contourf4[-1])
plt.tight_layout()
plt.savefig( os.path.join( FOLDER, f'conf_stack_space_{LABEL}_dihNCoCoN.pdf'), bbox_inches='tight')







# Plot #X (separate ang_PLANES (mirrored))
# ### 2b. Make cutoff plots: raw *vs* normalized by volume *vs* normalized by volume and # of DA pairs

# smaller font sizes!
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['figure.titlesize'] = 14 


fig, ax = plt.subplots(2,2,figsize=(10,9))


ax[0,0].set_title("Raw")
ax[0,0].set_ylabel(r'$\Theta_{ij}$ (deg)')
ax[0,0].set_xlabel('r$_{ij}$ ($\AA$)')
cf1 = ax[0,0].contourf(X3, Y3, H3.T, cmap ="Blues")
ax[0,0].set_ylim(X_for_contourf3[0],X_for_contourf3[-1])
ax[0,0].set_yticks(np.arange(0, 91, step=10))
ax[0,0].set_ylim(Y_for_contourf3[0],Y_for_contourf3[-1])

fig.colorbar(cf1, ax=ax[0,0])


ax[0,1].set_title("Normalized by volume")
ax[0,1].set_ylabel(r'$\Theta_{ij}$ (deg)')
ax[0,1].set_xlabel('r$_{ij}$ ($\AA$)')
cf2 = ax[0,1].contourf(X3, Y3, HT_volu3, cmap ="Blues")
ax[0,1].set_ylim(X_for_contourf3[0],X_for_contourf3[-1])
ax[0,1].set_yticks(np.arange(0, 91, step=10))
ax[0,1].set_ylim(Y_for_contourf3[0],Y_for_contourf3[-1])

fig.colorbar(cf2, ax=ax[0,1])


ax[1,0].set_title(r"Normalized by volume and # of $ij$ dimers")
ax[1,0].set_ylabel(r'$\Theta_{ij}$ (deg)')
ax[1,0].set_xlabel('r$_{ij}$ ($\AA$)')
cf3 = ax[1,0].contourf(X3, Y3, HT_norm3,cmap ="Blues", levels = np.arange(0,1.1,0.1), antialiased=True)#, levels=levels)
ax[1,0].set_ylim(X_for_contourf3[0],X_for_contourf3[-1])
ax[1,0].set_yticks(np.arange(0, 91, step=10))
ax[1,0].set_ylim(Y_for_contourf3[0],Y_for_contourf3[-1])

fig.colorbar(cf3, ax=ax[1,0])


ax[1,1].set_title("")
cf4 = ax[1,1].contourf(X3, Y3, HT_norm3,cmap ="Blues", levels = np.arange(0,100,10))
fig.colorbar(cf3, ax=ax[1,1])
ax[1,1].set_yticks(np.arange(0, 90, step=90))
ax[1,1].set_xticks(np.arange(0, 16, step=16))



plt.tight_layout()
plt.savefig( os.path.join( FOLDER, f'conf_stack_space_{LABEL}-raw-vs-volu-vs-norm.pdf'), bbox_inches='tight')
    

