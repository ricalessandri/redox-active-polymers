#!/usr/bin/env python
# coding: utf-8
"""
We first determine the dielectric constant of the mixture.
We compute lambda_out in atomi units first (hence, we can ignore the 4*pi*epsilon_0 = 1 and we need to convert ang to bohr).
Then, we convert the outer reorg. energy from Hartree back to eV.
"""

import numpy as np
import csv


def Landau_Lifshitz_Looyenga(eps_nonpolar, eps_polar, frac_polar):
    """
    Eq. (4) of [10.1002/bbpc.19910950801].
    """
    eps_mixture = eps_nonpolar + eps_nonpolar**(1./3) + (frac_polar * (eps_polar**(1./3) - eps_nonpolar**(1./3)))
    return eps_mixture 


def Debye(eps_nonpolar, eps_polar, frac_polar):
    """
    Eq. (1) of [10.1002/bbpc.19910950801].
    """
    xx = frac_polar * ((eps_polar-1)/(eps_polar+2)) + (1-frac_polar) * ((eps_nonpolar-1)/(eps_nonpolar+2))
    eps_mixture = (2*xx + 1) / (1. - xx)
    return eps_mixture 


def lambda_outer(delta_e, r_D, r_A, r_DA, n, eps):
    """
    Eq. XX of ... 
    In *atomic units*.

    We compute lambda_out in atomi units first. Hence, we can ignore the 4*pi*epsilon_0 and convert ang to bohr.
    Then, we convert the outer reorg. energy from Hartree back to eV.
    """
    return  delta_e**2 * ( 1./(2*r_D) + 1./(2*r_A) - 1/r_DA) * (1/n**2 - 1/eps) # atomic units


# Input the system in the below dictionary
r_D_eps_refndx_dict = {'NMPHTH_DME': (3.74,  7.2, 1.379), # cavity radius of NMPHTH (ang), relative permettivity of DME, refractive index of DME 
                       'TEMPO_H2O' : (3.94, 80.1, 1.334)} # cavity radius of TEMPO  (ang), relative permettivity of H2O, refractive index of H2O

# Parameters
delta_e = 1 # e
eps_nitroxide_polymers =  3.0  # [10.1021/jacs.1c02571]


for label, (r_D, eps_solvent, refractive_ndx) in r_D_eps_refndx_dict.items():

    # Epsilons 
    epsilon_with_05percent_DME = Debye(eps_nitroxide_polymers, eps_solvent, 0.05) # using model from Debye 
    epsilon_with_10percent_DME = Debye(eps_nitroxide_polymers, eps_solvent, 0.10) # using  
    epsilon_with_20percent_DME = Debye(eps_nitroxide_polymers, eps_solvent, 0.20) # using  

    starting_eps   = [eps_nitroxide_polymers, epsilon_with_05percent_DME, epsilon_with_10percent_DME, epsilon_with_20percent_DME, eps_solvent]   # dimethoxyethane (dimensionless)

    r_D = r_D * 1.8897259886  # ang-to-bohr
    r_A = r_D
    r_DA = r_D + r_A

    epsilon_values = starting_eps + [10**i for i in range(3, 5)]  # Values of epsilon from pure polymer to 10^5

    csv_filename = f"lambda_outer_values_{label}.csv"

    with open(csv_filename, 'w', newline='') as csvfile:
        fieldnames = ['epsilon', 'lambda_outer']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()  # Write the header row

        for epsilon in epsilon_values:
            lambda_out = lambda_outer(delta_e, r_D, r_A, r_DA, refractive_ndx, epsilon)
            lambda_out = lambda_out * 27.2114  # hartree to eV
            writer.writerow({'epsilon': epsilon, 'lambda_outer': lambda_out})

    print(f"CSV file '{csv_filename}' has been created.")

