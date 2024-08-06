#!/usr/bin/env python
# coding: utf-8
"""
USAGE: 
  python log10overlap_to_ALMOeV.py --overlaps_file overlaps.csv
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse


def from_log_overlap_to_ec_in_eV(slope, intercept, log10_of_overlap):
    """
    """
    return 10**((log10_of_overlap - intercept)/slope)

def from_log_overlap_to_log_ec(slope, intercept, log10_of_overlap):
    """
    """
    return (log10_of_overlap - intercept)/slope

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert the log10 overlaps to ALMO electronic couplings (eV).")
    parser.add_argument("--overlaps_file", type=str  , required=True, help="Path to a file containing coupling values")
    parser.add_argument('--verbose'      ,             default=True , help='True if you want verbose output (default behavior)')
    parser.add_argument('--slope'        , type=float, required=True, help='slope     of linear relationship between log(overlaps) and log(couplings)')
    parser.add_argument('--intercept'    , type=float, required=True, help='intercept of linear relationship between log(overlaps) and log(couplings)')

    args = parser.parse_args()
    VERBOSE = args.verbose
    SLOPE     = args.slope     #  0.8091639975519355 for NMPHTH [https://doi.org/10.1021/jacsau.4c00276]
    INTERCEPT = args.intercept # -1.5699211355668465 for NMPHTH [https://doi.org/10.1021/jacsau.4c00276]

    # Read in the overlaps
    if args.overlaps_file:
        data = np.loadtxt(args.overlaps_file, delimiter=',', dtype=str)
        IDs, overlaps = data[:, 0], data[:, 1].astype(float)

    if VERBOSE == True:
        print(f'Linear fit: y_fit = {SLOPE}*x_fit + {INTERCEPT}')

    couplings = []
    for overlap in overlaps:
        coupling_in_eV = from_log_overlap_to_ec_in_eV(SLOPE, INTERCEPT, overlap) 
        coupling = coupling_in_eV
        if VERBOSE == True:
            print(f"Coupling : {coupling:.2e} eV (for overlap {overlap:.2e})")
        couplings.append(coupling)

    couplings = np.array(couplings)
    IDs = np.array(IDs)
    if VERBOSE == True:
        print(f' len of couplings is {len(couplings)}; len of IDs is {len(IDs)}')

    # Save IDs and couplings to output file
    with open('couplings.csv', 'w') as output:
        for ID, coupling in zip(IDs, couplings):
            output.write(f"{ID},{coupling:.7f}\n")

    log_couplings = []
    for overlap in overlaps:
        log_coupling_in_eV = from_log_overlap_to_log_ec(SLOPE, INTERCEPT, overlap) 
        if VERBOSE == True:
            print(f"np.log10(coupling) : {coupling:.2e} meV (for overlap {overlap:.2e})")
        log_couplings.append(log_coupling_in_eV)

    log_couplings = np.array(log_couplings)
    IDs = np.array(IDs)
    if VERBOSE == True:
        print(f' len of log_couplings is {len(log_couplings)}; len of IDs is {len(IDs)}')

    # Save IDs and couplings to output file
    with open('log_couplings.csv', 'w') as output:
        for ID, log_coupling in zip(IDs, log_couplings):
            output.write(f"{ID},{log_coupling:.7f}\n")

