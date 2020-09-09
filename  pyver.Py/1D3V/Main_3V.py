#import os
import sys
#import numpy as np
#import matplotlib.pyplot as plt

import Species_3V as S
#import PIC_1D as P1D
#import PIC_1D_rel as P1D
import PIC_1D3V as P1D

#"""
###############################################################################
#                                                                             #
#                                two Stream                                   #
#                                                                             #
###############################################################################

#input parameters for the 1D grid
N_grid = 64.  # input: Number of grid cells
L = 10.0  # input: length of the domain


charge_to_mass = 1.  # q_p/m_p
thermal_velocity = 0.01 # input
rho_zero = 1. # background charge density
nppc = 400. # input: Number of particles per cell

v_zero = 0.1
thermal_velocity = v_zero/10. # input

#define the different particle species
Electrons = S.Species_3V(-1., thermal_velocity, rho_zero, nppc, color='red')
Positrons = S.Species_3V(1., thermal_velocity, rho_zero, nppc, color='k')


delta_t = L/(4.*N_grid)

PIC = P1D.PIC_1D3V(N_grid, L, delta_t, Species = [Electrons, Positrons], time_steps = 8000, folder='3_vector_test_2.0')
PIC.EM_two_Stream(v_zero)
#PIC.Bx_grid_center = np.ones(10)
#PIC.Bx_grid_center[5] = 2.
PIC.Initiate_Plasma(snapshot = 400, Output_xv = 10)

#"""
sys.exit()
#"""     
###############################################################################
#                                                                             #
#                         two Stream: Vary Parameters                         #
#                                                                             #
###############################################################################
# NG for 64 128 256
# nppc: 10 50 100
# delta t: delta_X/2 delta_X/4 and delta_x/8 
# v_0 = 0.01, 0.1, 0.001
# take L = 10

L = 10.

rho_zero = 1.
charge_to_mass = 1.

N_grid = [128.] #[64., 128., 256]
nppc =  [50.] #[10., 50., 100.]
v_zero = 0.01  # [0.01, 0.05, 0.1]
thermal_velocity = v_zero/10

#delta_x = Domain_Length/N_grid
t_factor = [1./2., 1./4., 1./8.]

for i in N_grid:
    for j in nppc:
        for k in t_factor:
            delta_x = L/i
            delta_t = k*delta_x
            print delta_t, i, j, int(10000/2*(1./k))

            #define the different particle species
            Electrons = S.Species_3V(-charge_to_mass, thermal_velocity, rho_zero, j, color='red')
            Positrons = S.Species_3V(charge_to_mass, thermal_velocity, rho_zero, j, color='k')
        
            Output_folder='N_grid_' + str(i) + '_nppc_' + str(j) + '_delta_t_' + str(delta_t)
            PIC = P1D.PIC_1D3V(i, L, delta_t, Species = [Electrons, Positrons], time_steps = int(10000/2*(1./k)), folder=Output_folder)
            PIC.EM_two_Stream(v_zero)
            PIC.Initiate_Plasma(Output_xv= 10)     #snapshot = 20
            
#"""
#"""     
###############################################################################
#                                                                             #
#                         two Stream: Vary Parameters                         #
#                                                                             #
###############################################################################
# NG for 64 128 256
# nppc: 10 50 100
# delta t: delta_X/2 delta_X/4 and delta_x/8 
# v_0 = 0.01, 0.1, 0.001
# take L = 10

L = 10.

rho_zero = 1.
charge_to_mass = 1.

N_grid = [128.] #[64., 128., 256]
nppc =  [10., 50., 100.]
v_zero = 0.01  # [0.01, 0.05, 0.1]
thermal_velocity = v_zero/10

#delta_x = Domain_Length/N_grid
t_factor = [1./4.]

for i in N_grid:
    for j in nppc:
        for k in t_factor:
            delta_x = L/i
            delta_t = k*delta_x
            print delta_t, i, j, int(10000/2*(1./k))

            #define the different particle species
            Electrons = S.Species_3V(-charge_to_mass, thermal_velocity, rho_zero, j, color='red')
            Positrons = S.Species_3V(charge_to_mass, thermal_velocity, rho_zero, j, color='k')
        
            Output_folder='N_grid_' + str(i) + '_nppc_' + str(j) + '_delta_t_' + str(delta_t)
            PIC = P1D.PIC_1D3V(i, L, delta_t, Species = [Electrons, Positrons], time_steps = int(10000/2*(1./k)), folder=Output_folder)
            PIC.EM_two_Stream(v_zero)
            PIC.Initiate_Plasma(Output_xv= 10)     #snapshot = 20       
#"""
#"""            
###############################################################################
#                                                                             #
#                         two Stream: Vary Parameters                         #
#                                                                             #
###############################################################################
# NG for 64 128 256
# nppc: 10 50 100
# delta t: delta_X/2 delta_X/4 and delta_x/8 
# v_0 = 0.01, 0.1, 0.001
# take L = 10

L = 10.

rho_zero = 1.
charge_to_mass = 1.

N_grid = [64., 128., 256.] #[64., 128., 256]
nppc =  [50.] #[10., 50., 100.]
v_zero = 0.01  # [0.01, 0.05, 0.1]
thermal_velocity = v_zero/10

#delta_x = Domain_Length/N_grid
t_factor = [1./4.]

for i in N_grid:
    for j in nppc:
        for k in t_factor:
            delta_x = L/128.
            delta_t = k*delta_x
            print delta_t, i, j, int(10000/2*(1./k))

            #define the different particle species
            Electrons = S.Species_3V(-charge_to_mass, thermal_velocity, rho_zero, j, color='red')
            Positrons = S.Species_3V(charge_to_mass, thermal_velocity, rho_zero, j, color='k')
        
            Output_folder='N_grid_' + str(i) + '_nppc_' + str(j) + '_delta_t_' + str(delta_t)
            PIC = P1D.PIC_1D3V(i, L, delta_t, Species = [Electrons, Positrons], time_steps = int(10000/2*(1./k)), folder=Output_folder)
            PIC.EM_two_Stream(v_zero)
            PIC.Initiate_Plasma(Output_xv= 10)     #snapshot = 20 
#"""
#"""            
###############################################################################
#                                                                             #
#                       two Stream: Variations of v_zero                      #
#                                                                             #
###############################################################################
# v_0 = 0.01, 0.1, 0.001

L = 10.
rho_zero = 1.
charge_to_mass = 1.

N_grid = 128. #[32., 64., 128.]
nppc =  50. #[10., 50., 100.]
v_zero = [0.01, 0.05, 0.1]

delta_t = L/(N_grid*4)

for i in v_zero:
    thermal_velocity = i/10

    #define the different particle species
    Electrons = S.Species_3V(-charge_to_mass, thermal_velocity, rho_zero, nppc, color='red')
    Positrons = S.Species_3V(charge_to_mass, thermal_velocity, rho_zero, nppc, color='k')
        
    Output_folder='v_zero_' + str(i)
    PIC = P1D.PIC_1D3V(N_grid, L, delta_t, Species = [Electrons, Positrons], time_steps = 40000, folder=Output_folder)
    PIC.EM_two_Stream(i)
    PIC.Initiate_Plasma()
#"""