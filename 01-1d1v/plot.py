#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 15:57:12 2018

@author: luca_pezzini
"""

import os
import sys
import numpy as np
from textwrap import wrap
import matplotlib.pyplot as plt

#v = [3]
#Time = [[None] for k in range(len(v))]
#Energy = [[None] for k in range(len(v))]
Output_folder = os.getcwd() + '/Output/'


plt.xlabel(r'Time [$\omega_p^{-1}$]', fontsize = 20) 
plt.ylabel(r'$E_{Elec}$ []', fontsize = 20) 
 
#Input_folder=os.getcwd()
            
#Time[i] = np.loadtxt(os.path.normpath(Input_folder + 'energy.dat), skiprows=0, unpack=True)[0]
#Energy[i] = np.loadtxt(os.path.normpath(Input_folder + 'energy.dat'), skiprows=0, unpack=True)[1] 

f, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3, sharex='col', sharey='row')
fig = plt.figure(figsize=(10,7)) 
ax1.Time = np.loadtxt(os.path.normpath('energy.dat'), skiprows=0, unpack=True)[0]
ax1.Energy = np.loadtxt(os.path.normpath('energy.dat'), skiprows=0, unpack=True)[1]

plt.semilogy(Time, Energy, label=r'$v_b$: ')

        
plt.title(r"Variation of $v_b$ with" + '\n' + r"$N_{grid}$=64, nppc=50 and $\Delta$t=$\Delta$x/4",  fontsize=25)        
plt.legend(loc='lower right')
fig.tight_layout()
plt.savefig('V_zero_variation.png')      
plt.show()        
plt.close(fig)