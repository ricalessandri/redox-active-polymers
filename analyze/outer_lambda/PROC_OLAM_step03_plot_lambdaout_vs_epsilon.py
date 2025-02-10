#!/usr/bin/env python
# coding: utf-8


import matplotlib.pyplot as plt
import matplotlib
import matplotlib as mpl
import numpy as np
import seaborn as sns

matplotlib.rcParams.update({'font.size': 18})  # You can adjust the value as needed
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


plt.figure(figsize=(6.4, 4.8))

data_file_1 = f'lambda_outer_values_NMPHTH_DME.csv'
data_file_2 = f'lambda_outer_values_TEMPO_H2O.csv'
dataX1 = np.loadtxt(data_file_1, delimiter=',', skiprows=1, usecols=(0))
dataY1 = np.loadtxt(data_file_1, delimiter=',', skiprows=1, usecols=(1))
dataX2 = np.loadtxt(data_file_2, delimiter=',', skiprows=1, usecols=(0))
dataY2 = np.loadtxt(data_file_2, delimiter=',', skiprows=1, usecols=(1))

plt.scatter(dataX1,dataY1, marker='o', color=sns.color_palette()[0]) #, color='C0')
plt.scatter(dataX2,dataY2, marker='o', color=sns.color_palette()[1]) #, color='C1')


plt.axhline(y=dataY1[2], color='C0', linestyle='--'    , label=r'$\lambda_{\mathrm{out, NMePh}}=$'+ f'{round(dataY1[2],3)} eV (10%DME)')
plt.axhline(y=dataY2[2], color='C1', linestyle='--'    , label=r'$\lambda_{\mathrm{out, TEMPO}}=$'+ f'{round(dataY2[2],3)} eV (10%H2O)')

plt.xlabel(r'$\epsilon_r$')
plt.ylabel(r'$\lambda_{\mathrm{out}}$ (eV)')

plt.xscale('log')
plt.ylim(0,1.1)

plt.legend(fontsize=12,loc='lower right')
plt.grid(True)
plt.tight_layout()
plt.savefig(f'outer_lambda_vs_epsilon.pdf')

