#!/usr/bin/env python
# coding: utf-8


# USAGE:
# python just_plot_data_and_fit.py --what overlap

# =============================================== #
# 0) Import modules & matplotlib settings         #
# =============================================== #
import numpy as np
import matplotlib.pylab as plt
from sklearn.metrics import mean_squared_error, r2_score
import time
import matplotlib as mpl
import argparse
import sys
from cycler import cycler
from scipy.stats import skewnorm
from scipy.signal import find_peaks


mpl.use('Agg') # set a non-interactive matplotlib backend


plt.rcParams['font.family'] = 'sans'
plt.rcParams['font.size'] = 13  #15  # 18
plt.rcParams['axes.labelsize'] = 15  # 20
plt.rcParams['axes.labelweight'] = 'normal'
plt.rcParams['xtick.labelsize'] = 15  # 18
plt.rcParams['ytick.labelsize'] = 15  # 18
plt.rcParams['legend.fontsize'] = 15  # 18
mpl.rcParams['axes.linewidth'] = 1.5  # set the thickness of the graph border globally


Tol_bright = ('#4477AA', '#EE6677', '#AA3377', '#BBBBBB') 
plt.rcParams["axes.prop_cycle"] = cycler(color=Tol_bright)


# Parse the arguments
parser = argparse.ArgumentParser(description='Just plot the data')
parser.add_argument('--what'      , required=True      , type=str  , help='what it is that we are plotting?')
parser.add_argument('--v-type'    , default='raw'      , type=str  , help='type of coupling data provided: *raw* or *log')
parser.add_argument('--verbose'   , default=False      , type=bool , help='more printing...')
parser.add_argument('--csvfile'   , default='default'  , type=str  , help='name of CSV file if not default')
parser.add_argument('--units'     , default='au'       , type=str  , help='units; only relevant to MO energies; either "au" or "eV" accepted')
parser.add_argument('--label'     , default=''         , type=str  , help='*optional* label to be printed at the top of the graph')
parser.add_argument('--ylim'      , default=None       , type=float, help='y-axis upper limit')

args = parser.parse_args()

WHAT             = args.what
OVERLAP_TYPE     = args.v_type
VERBOSE          = args.verbose
CSVFILE          = args.csvfile
UNITS            = args.units
LABEL            = args.label
YLIM             = args.ylim
if CSVFILE == 'default':
    CSVFILE = f'{WHAT}s.csv'
print(f"\nParameters (taking into account user input): what={WHAT}, verbose={VERBOSE}; csvfile={CSVFILE}.")
if not LABEL:
    print('We have no label')





# =============================================== #
# 1) Data loading                                 #
# =============================================== #

def load_overlaps(CSVFILE, OVERLAP_TYPE):
    """
    """
    data = np.genfromtxt(CSVFILE, delimiter=',')
    IDs      = np.transpose(data)[0].astype(int) # IDs are integers
    if "predicted" in CSVFILE:
        overlaps = np.transpose(data)[1]             # this is a 'overlaps_predicted.csv' so there's only 2 columns 
    else:
        overlaps = np.transpose(data)[1]             # load the "raw" overlap IJ 
        overlaps = np.transpose(data)[2]             # load the "raw" overlap JI
        overlaps = np.transpose(data)[3]             # load the "raw" overlap (IJ+JI)/2
    print(f'Sizes of the IDs ({len(IDs)}) and overlap ({len(overlaps)}) vectors ("as recevied").')
    
    if OVERLAP_TYPE == 'raw':
        indices_of_zero_overlap = np.argwhere(overlaps == 0.0) # get indices of zero overlap
        IDs  = np.delete(IDs,  indices_of_zero_overlap)
        overlaps = np.delete(overlaps, indices_of_zero_overlap)
        print(f'Sizes of the IDs ({len(IDs)}) and overlap ({len(overlaps)}) vectors (after removing zeros).')
        indices_of_nan_overlap = np.argwhere(np.isnan(overlaps)) # get indices of NaN overlap
        IDs  = np.delete(IDs,  indices_of_nan_overlap)
        overlaps = np.delete(overlaps, indices_of_nan_overlap)
        print(f'Sizes of the IDs ({len(IDs)}) and overlap ({len(overlaps)}) vectors (after removing NaNs).')

        # Use the log of the overlaps
        dataY = np.log10(np.abs(overlaps))

    elif OVERLAP_TYPE == 'log':
        dataY = overlaps

    else:
        sys.exit('In what format is the overlap data? *raw* or already in *log* form? I only know how to deal with these two cases...')

    return IDs, dataY


def remove_zeros(IDs, dataY):
    """
    """
    indices_of_zero_dataY = np.argwhere(dataY == 0.0) # get indices of zeros in the dataY
    IDs  = np.delete(IDs,  indices_of_zero_dataY)
    dataY = np.delete(dataY, indices_of_zero_dataY)
    print(f'Sizes of the IDs ({len(IDs)}) and dataY ({len(dataY)}) vectors (after removing zeros).')
    return IDs, dataY


def fit_skewnorm_to_the_data(dataY, initial_guess):
    """
    """
    params = skewnorm.fit(dataY, initial_guess[0], loc=initial_guess[1], scale=initial_guess[2]) # a, log, scale = skewness, mean, st_dev
    return params


def extract_mean_as_peak_position(x, params, dataY, bin_edges):
    """
    """
    pdf_values = skewnorm.pdf(x, *params)*len(dataY) * (bin_edges[1] - bin_edges[0])
    # Find peaks in the PDF using scipy.signal.find_peaks
    peaks, _ = find_peaks(pdf_values)
    # Select the peak with the highest PDF value
    peak_position = x[peaks[np.argmax(pdf_values[peaks])]]
    return peak_position


if WHAT=="overlap":

    IDs , dataY = load_overlaps(CSVFILE, OVERLAP_TYPE)

    print(f'Sizes of the IDs ({len(IDs)}) and log10(overlap) ({len(dataY)}) vectors (after removing zeros).')
else:
    sys.exit("WHAT did you pass to the script? I only understand *overlap*.")


IDs, dataY = remove_zeros(IDs, dataY)


if WHAT=="overlap":
    f = plt.figure(figsize=(5,5), dpi=300)
    ax = plt.axes()
    ax.tick_params(direction='in',width=1.5)
    if LABEL:
        plt.title(f"{LABEL}")
        #plt.title(r"Distribution of $\log_{10}\langle\phi_{SOMO}|\phi_{LUMO}\rangle$")

    bin_edges = np.linspace(-8,-1,100)
    hist, binss, _ = plt.hist(dataY, bins=bin_edges,density=False,alpha=0.0,color='black')

    # Fit a skewed gaussian
    initial_guess = [100, np.mean(dataY), np.std(dataY)]
    params = fit_skewnorm_to_the_data(dataY, initial_guess)
    x = (binss[:-1] + binss[1:]) / 2  # Use the bin centers
    peak_position = extract_mean_as_peak_position(x, params, dataY, bin_edges)

    # Actual plot
    hist, binss, _ = plt.hist(dataY, bins=bin_edges, label=r'$\mu=$'+f'{round(peak_position,3)} (pred.)',density=False,alpha=0.7)
    ax.plot(x, skewnorm.pdf(x, *params)*len(dataY) * (bin_edges[1] - bin_edges[0]), linewidth=3, color='darkblue')

    plt.xlabel(r"$\log_{10}\langle\phi_{SOMO}|\phi_{LUMO}\rangle$")
    plt.xlim(-7.5,-1)
    plt.legend(loc='upper left',fontsize=12)
    plt.tight_layout()
    #f.savefig(f'0_data_{WHAT}s.png', bbox_inches='tight')
    f.savefig(f'0_data_{WHAT}s.pdf', bbox_inches='tight')

