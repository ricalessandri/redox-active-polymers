#!/usr/bin/env python
# coding: utf-8
"""
DESCRIPTION:
    TBD
USAGE:
    TBD
"""

import numpy as np
import matplotlib.pylab as plt
from sklearn.metrics import mean_squared_error, r2_score
import time
import matplotlib as mpl
import argparse
import sys
import seaborn as sns
from cycler import cycler
from scipy.stats import skewnorm
from scipy.signal import find_peaks

mpl.rcParams.update({'font.size': 18})  # You can adjust the value as needed
sns.set_palette("colorblind") # Set the colorblind-friendly palette
mpl.use('Agg') # set a non-interactive matplotlib backend
plt.rcParams['font.family'] = 'sans'
plt.rcParams['font.size'] = 15 # 16
plt.rcParams['axes.labelsize'] = 17 # 16
plt.rcParams['axes.labelweight'] = 'normal'
plt.rcParams['xtick.labelsize'] = 17 # 18
plt.rcParams['ytick.labelsize'] = 17 # 18
plt.rcParams['legend.fontsize'] = 17 # 18 
mpl.rcParams['axes.linewidth'] = 1.5  # set the thickness of the graph border globally
Tol_bright = ('#4477AA', '#EE6677', '#AA3377', '#BBBBBB') 
plt.rcParams["axes.prop_cycle"] = cycler(color=Tol_bright)


# Parse the arguments
parser = argparse.ArgumentParser(description='Just plot the data')
parser.add_argument('--what'      , required=True      , type=str  , help='what it is that we are plotting?')
parser.add_argument('--verbose'   , default=False      , type=bool , help='more printing...')
parser.add_argument('--csvfile01' , required=True      , type=str  , help='name of CSV file if not default')
parser.add_argument('--csvfile02' , default='unknown'  , type=str  , help='name of CSV file if not default')
parser.add_argument('--csvfile03' , default='unknown'  , type=str  , help='name of CSV file if not default')
parser.add_argument('--reference' , default=False      , type=str  , help='plot against reference data?')
parser.add_argument('--csvfileREF', default=False      , type=str  , help='name of CSV file containing the REFerence distribution')
parser.add_argument('--boltzmann' , default=False      , type=float, help='Boltzmann-averaged value provided? (for couplings only)')
parser.add_argument('--units'     , default='au'       , type=str  , help='units; only relevant to MO energies; either "au" or "eV" accepted')
parser.add_argument('--label'     , default=''         , type=str  , help='*optional* label to be printed at the top of the graph')
parser.add_argument('--polymer'   , default=''         , type=str  , help='*optional* label to be printed at the top of the graph')
parser.add_argument('--percent'   , default=''         , type=str  , help='*optional* percent of solvent')

args = parser.parse_args()

WHAT             = args.what
VERBOSE          = args.verbose
CSVFILE_1        = args.csvfile01
CSVFILE_2        = args.csvfile02
CSVFILE_3        = args.csvfile03
REF              = args.reference
REFCSVFILE       = args.csvfileREF
BOLTZMANN        = args.boltzmann
UNITS            = args.units
LABEL            = args.label
POLYMER          = args.polymer
PERCENT          = args.percent
if REF:
    if not REFCSVFILE:
        sys.exit(f"Error: requested to plot against reference data (REF = {REF}) but no CSV file which the reference data was provided. Exiting.'")
print(f"\nParameters (taking into account user input): what={WHAT}, verbose={VERBOSE}; csvfile(s)={CSVFILE_1},{CSVFILE_2},{CSVFILE_3}; Boltzmann={BOLTZMANN}.")
if not LABEL:
    print('We have no label')



# =============================================== #

def load_couplings(CSVFILE):
    """
    """
    data = np.genfromtxt(CSVFILE, delimiter=',')
    IDs       = np.transpose(data)[0].astype(int) # IDs are integers
    couplings = np.transpose(data)[1]
    print(f'Sizes of the IDs ({len(IDs)}) and couplings ({len(couplings)}) vectors ("as recevied").')
    
    indices_of_zero_coupling = np.argwhere(couplings == 0.0) # get indices of zeros in the overlap
    IDs       = np.delete(IDs,  indices_of_zero_coupling)
    couplings = np.delete(couplings, indices_of_zero_coupling)
    print(f'Sizes of the IDs ({len(IDs)}) and couplings ({len(couplings)}) vectors (after removing zeros).')

    dataY = couplings

    return IDs, dataY


def remove_zeros(IDs, dataY):
    """
    """
    indices_of_zero_dataY = np.argwhere(dataY == 0.0) # get indices of zeros in the dataY
    IDs  = np.delete(IDs,  indices_of_zero_dataY)
    dataY = np.delete(dataY, indices_of_zero_dataY)
    print(f'Sizes of the IDs ({len(IDs)}) and dataY ({len(dataY)}) vectors (after removing zeros).')
    return IDs, dataY


def load_couplings_and_remove_zeros(CSVFILE):
    """
    """
    IDs, dataY = load_couplings(CSVFILE)
    IDs, dataY = remove_zeros(IDs, dataY)
    print(f'np.min(dataY) = {np.min(dataY)} log10[eV]; np.max(dataY) = {np.max(dataY)} log10[eV]')
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


# Before starting, CHECK
if WHAT != "coupling":
    sys.exit("WHAT did you pass to the script? I only understand *coupling*.")


IDs_1, dataY_1 = load_couplings_and_remove_zeros(CSVFILE_1)
initial_guess_1 = [10, np.mean(dataY_1)+0.5, np.std(dataY_1)]
params_1 = fit_skewnorm_to_the_data(dataY_1, initial_guess_1)
if CSVFILE_2 != "unknown":
    IDs_2, dataY_2 = load_couplings_and_remove_zeros(CSVFILE_2)
    initial_guess_2 = [10, np.mean(dataY_2)+0.5, np.std(dataY_2)]
    params_2 = fit_skewnorm_to_the_data(dataY_2, initial_guess_2)
if CSVFILE_3 != "unknown":
    IDs_3, dataY_3 = load_couplings_and_remove_zeros(CSVFILE_3)
    initial_guess_3 = [10, np.mean(dataY_3)+0.8, np.std(dataY_3)]
    params_3 = fit_skewnorm_to_the_data(dataY_3, initial_guess_3)


# Actual plot
f = plt.figure(figsize=(4.8, 4.8), dpi=300)
ax = plt.axes()
ax.tick_params(direction='in',width=1.5)
if LABEL:
    plt.title(f"{LABEL}")

if CSVFILE_2 != "unknown" and CSVFILE_3 != "unknown":
    mu_1 = np.mean(dataY_1); mu_2 = np.mean(dataY_2); mu_3 = np.mean(dataY_3)
    print(f' OLD mu_1 = {mu_1}; mu_2 = {mu_2}; mu_3 = {mu_3}')
    print(f' OLD mu_1 = {params_1[1]}; mu_2 = {params_2[1]}; mu_3 = {params_3[1]}')
    print(f' OLD mu_1 = {10**mu_1}; mu_2 = {10**mu_2}; mu_3 = {10**mu_3}')
    print(f' OLD mu_1 = {10**params_1[1]}; mu_2 = {10**params_2[1]}; mu_3 = {10**params_3[1]}')
    sigma_1 = np.std(dataY_1); sigma_2 = np.std(dataY_2); sigma_3 = np.std(dataY_3)
    print(f' NP.STD    sigma_1 = {sigma_1}; sigma_2 = {sigma_2}; sigma_3 = {sigma_3}')
    print(f' PARAMS[2] sigma_1 = {params_1[2]}; sigma_2 = {params_2[2]}; sigma_3 = {params_3[2]}')

bin_edges       = np.linspace(-6,1,280) # better (sharper) and more sensitivity for the fit

# The following 3 lines are just for obtaining bins_1, bins_2, and bins_3, but they'll plot nothing (because of alpha=0.0)
hist_1      , bins_1      , _ = plt.hist(dataY_1, bins=bin_edges      , density=False,alpha=0.0, color='black')
x_1 = (bins_1[:-1] + bins_1[1:]) / 2  # Use the bin centers
peak_position_1 = extract_mean_as_peak_position(x_1, params_1, dataY_1, bin_edges)
average_Vij_from_peak_1  = 10**peak_position_1 # eV
print(f"The range_Vij_incl_std_1 is from {10**(peak_position_1+params_1[2])} to {10**(peak_position_1-params_1[2])}")
print(f"The peak is at {peak_position_1}, which means {10**(peak_position_1)} eV")
if CSVFILE_2 != "unknown":
    hist_2, bins_2, _ = plt.hist(dataY_2, bins=bin_edges, density=False,alpha=0.0, color='black')
    x_2 = (bins_2[:-1] + bins_2[1:]) / 2  # Use the bin centers
    peak_position_2 = extract_mean_as_peak_position(x_2, params_1, dataY_2, bin_edges)
    average_Vij_from_peak_2  = 10**peak_position_2 # eV
    print(f"The range_Vij_incl_std_2 is from {10**(peak_position_2+params_2[2])} to {10**(peak_position_2-params_2[2])}")
    print(f"The peak is at {peak_position_2}, which means {10**(peak_position_2)} eV")
if CSVFILE_3 != "unknown":
    hist_3, bins_3, _ = plt.hist(dataY_3, bins=bin_edges, density=False,alpha=0.0, color='#F4B942')
    x_3 = (bins_3[:-1] + bins_3[1:]) / 2  # Use the bin centers
    peak_position_3 = extract_mean_as_peak_position(x_3, params_3, dataY_3, bin_edges)
    average_Vij_from_peak_3  = 10**peak_position_3 # eV
    print(f"The range_Vij_incl_std_3 is from {10**(peak_position_3+params_3[2])} to {10**(peak_position_3-params_3[2])}")
    print(f"The peak is at {peak_position_3}, which means {10**(peak_position_3)} eV")


if CSVFILE_2 != "unknown" and CSVFILE_3 != "unknown":
    print("Position of the peak:", peak_position_1        , peak_position_2        , peak_position_3    )
    print("Means (in eV)       :", average_Vij_from_peak_1, average_Vij_from_peak_2, average_Vij_from_peak_3)

with open (f'0_data_{WHAT}s.txt', 'w') as txtout:
    txtout.write('polymer SoC percent mean_Vij (eV)\n')
    txtout.write('---------------------------------------\n')
    txtout.write(f'{POLYMER}    000    {PERCENT}    {average_Vij_from_peak_1}\n') 
    if CSVFILE_2 != "unknown":
        txtout.write(f'{POLYMER}    020    {PERCENT}    {average_Vij_from_peak_2}\n')
    if CSVFILE_3 != "unknown":
        txtout.write(f'{POLYMER}    060    {PERCENT}    {average_Vij_from_peak_3}\n')

if POLYMER == "PTMA":
    color_1 = "#ad4059"
else:
    color_1 = '#4059AD'

hist_1, bins_1, _ = plt.hist(dataY_1, bins=bin_edges, label=r'  0%; $\langle V_{ij}\rangle=$'+f'{int(round(average_Vij_from_peak_1*1000,0))} meV',
                             density=False,alpha=0.4, color=color_1)
if CSVFILE_2 != "unknown":
    hist_2, bins_2, _ = plt.hist(dataY_2, bins=bin_edges, label=r'20%; $\langle V_{ij}\rangle=$' +f'{int(round(average_Vij_from_peak_2*1000,0))} meV',
                                 density=False,alpha=0.4, color='#9A8978')
if CSVFILE_3 != "unknown":
    hist_3, bins_3, _ = plt.hist(dataY_3, bins=bin_edges, label=r'60%; $\langle V_{ij}\rangle=$' +f'{int(round(average_Vij_from_peak_3*1000,0))} meV',
                                 density=False,alpha=0.4, color='#F4B942')

plt.annotate(f"{POLYMER}", (-5.8, 10 ), horizontalalignment='left', color='gray', weight='bold') # if bin_edges = np.linspace(-6,1,280)

# draw the pdf of the fitted skewnorm
ax.plot(x_1, skewnorm.pdf(x_1, *params_1)*len(dataY_1) * (bin_edges[1] - bin_edges[0]), color=color_1, linewidth=3)
if CSVFILE_2 != "unknown":
    ax.plot(x_2, skewnorm.pdf(x_2, *params_1)*len(dataY_2) * (bin_edges[1] - bin_edges[0]), color='#9A8978', linewidth=3)
if CSVFILE_3 != "unknown":
    ax.plot(x_3, skewnorm.pdf(x_3, *params_3)*len(dataY_3) * (bin_edges[1] - bin_edges[0]), color='#F4B942', linewidth=3)

plt.xlabel(r"$\log_{10}V_{ij}$ [$\log_{10}$(eV)]")
plt.ylabel(r"counts")
plt.xlim(-6,0.5)
plt.ylim(0,225) # if bin_edges = np.linspace(-6,1,280)
plt.legend(loc='upper left',fontsize=10)
plt.tight_layout()
#f.savefig(f'0_data_{WHAT}s.png', bbox_inches='tight')
f.savefig(f'0_data_{WHAT}s.pdf', bbox_inches='tight')

