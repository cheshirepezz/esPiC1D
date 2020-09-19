#import os
import sys
import numpy as np
import matplotlib.pyplot as plt

import Species as S
#import PIC_1D as P1D
import PIC_1D_rel as P1D
#"""

###############################################################################
#                                                                             #
#                               Relativistic                                  #
#                                                                             #
###############################################################################

#input parameters for the 1D grid
N_grid = 64.  # input: Number of grid cells
L = 100.  # input: length of the domain


charge_to_mass = 1.  # q_p/m_p
thermal_velocity = 0.2 # input
rho_zero = 1. # background charge density
nppc = 200. # input: Number of particles per cell

gamma_zero = 4

#define the different particle species
Electrons = S.Species(-1., thermal_velocity, rho_zero, nppc, color='red')
Positrons = S.Species(1., thermal_velocity, rho_zero, nppc, color='k')

delta_t = L/(8.*N_grid)

PIC = P1D.PIC_1D_rel(N_grid, L, delta_t, Species = [Electrons, Positrons], time_steps = 1000, folder='Relativistic_gif' )
PIC.relativistic_two_Stream(gamma_zero)
PIC.Initiate_Relativistic_Plasma(snapshot = 1000)

#"""
sys.exit()
#"""
###############################################################################
#                                                                             #
#                                two Stream                                   #
#                                                                             #
###############################################################################

#input parameters for the 1D grid
N_grid = 64.  # input: Number of grid cells
k_max = np.sqrt(3)*1.0/(2*0.1)
L =  4.15*2*np.pi/k_max  # input: length of the domain


charge_to_mass = 1.  # q_p/m_p
thermal_velocity = 0.01 # input
rho_zero = 1. # background charge density
nppc = 500. # input: Number of particles per cell

v_zero = 0.1

#define the different particle species
Electrons = S.Species(-1., thermal_velocity, rho_zero, nppc, color='red')
Positrons = S.Species(1., thermal_velocity, rho_zero, nppc, color='k')


delta_t = L/(8.*N_grid)

PIC = P1D.PIC_1D_rel(N_grid, L, delta_t, Species = [Electrons, Positrons], time_steps = 8000, folder='Two_Stream_test_gif_3.0')
PIC.two_Stream(v_zero)
PIC.Initiate_Plasma(snapshot = 800)

#"""
sys.exit()
#"""
###############################################################################
#                                                                             #
#                                Electrostatic                                #
#                                                                             #
###############################################################################

#input parameters for the 1D grid
N_grid = 32.  # input: Number of grid cells
L = 1.  # input: length of the domain

charge_to_mass = 1.  # q_p/m_p
thermal_velocity = 0.01 # input
rho_zero = 1. # background charge density
nppc = 10. # input: Number of particles per cell

#define the different particle species
Electrons = S.Species(charge_to_mass, thermal_velocity, rho_zero, nppc, name='Electrons')

delta_t = L/(2.*N_grid)


PIC = P1D.PIC_1D(N_grid, L, delta_t, Species = [Electrons], time_steps = 5000)

# velocity distribution at the start
plt.hist(Electrons.v_pos)
plt.show()

PIC.Initiate_Plasma(snapshot = 10)

# velocity distribution at the End
plt.hist(Electrons.v_pos)
plt.show()

#"""
#"""
###############################################################################
#                                                                             #
#                                two Stream                                   #
#                                                                             #
###############################################################################

#input parameters for the 1D grid
N_grid = 32.  # input: Number of grid cells
L = 10.0  # input: length of the domain


charge_to_mass = 1.  # q_p/m_p
thermal_velocity = 0.01 # input
rho_zero = 1. # background charge density
nppc = 100. # input: Number of particles per cell

v_zero = 0.1

#define the different particle species
Electrons = S.Species(-1., thermal_velocity, rho_zero, nppc, color='red')
Positrons = S.Species(1., thermal_velocity, rho_zero, nppc, color='k')


delta_t = L/(2.*N_grid)

PIC = P1D.PIC_1D_rel(N_grid, L, delta_t, Species = [Electrons, Positrons], time_steps = 300, folder='Two_Stream_test')
PIC.two_Stream(v_zero)
PIC.Initiate_Plasma(snapshot = 10)

#"""
#"""
###############################################################################
#                                                                             #
#                         two Stream: Vary Parameters                         #
#                                                                             #
###############################################################################
# NG for 32 64 128
# nppc: 10 50 100
# delta t: delta_X/2 delta_X/4 and delta_x/8 
# v_0 = 0.01, 0.1, 0.001
# take L = 10

L = 10.

rho_zero = 1.
charge_to_mass = 1.

N_grid = [32., 64., 128.]
nppc =  [10., 50., 100.]
v_zero = 0.01  # [0.001, 0.01, 0.1]
thermal_velocity = v_zero/10

#delta_x = Domain_Length/N_grid
t_factor = [1./2., 1./4., 1./8.]

for i in N_grid:
    for j in nppc:
        for k in t_factor:
            delta_x = L/i
            delta_t = k*delta_x
#            print delta_t, i, int(1500/2*(1./k))

            #define the different particle species
            Electrons = S.Species(-charge_to_mass, thermal_velocity, rho_zero, j, color='red')
            Positrons = S.Species(charge_to_mass, thermal_velocity, rho_zero, j, color='k')
        
            Output_folder='N_grid_' + str(i) + '_nppc_' + str(j) + '_delta_t_' + str(delta_t)
            PIC = P1D.PIC_1D_rel(i, L, delta_t, Species = [Electrons, Positrons], time_steps = int(1500/2*(1./k)), folder=Output_folder)
            PIC.two_Stream(v_zero)
            PIC.Initiate_Plasma(snapshot = 20)
            
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

N_grid = 512. #[32., 64., 128.]
nppc =  100. #[10., 50., 100.]
v_zero = [0.01, 0.05, 0.1]

delta_t = L/(N_grid*4)

for i in v_zero:
    thermal_velocity = i/10

    #define the different particle species
    Electrons = S.Species(-charge_to_mass, thermal_velocity, rho_zero, nppc, color='red')
    Positrons = S.Species(charge_to_mass, thermal_velocity, rho_zero, nppc, color='k')
        
    Output_folder='v_zero_' + str(i)
    PIC = P1D.PIC_1D_rel(N_grid, L, delta_t, Species = [Electrons, Positrons], time_steps = 2000, folder=Output_folder)
    PIC.two_Stream(i)
    PIC.Initiate_Plasma(snapshot = 20)

#"""
#"""

###############################################################################
#                                                                             #
#                               Relativistic                                  #
#                                                                             #
###############################################################################

#input parameters for the 1D grid
N_grid = 64.  # input: Number of grid cells
L = 100.  # input: length of the domain


charge_to_mass = 1.  # q_p/m_p
thermal_velocity = 0.2 # input
rho_zero = 1. # background charge density
nppc = 200. # input: Number of particles per cell

gamma_zero = 6

#define the different particle species
Electrons = S.Species(-1., thermal_velocity, rho_zero, nppc, color='red')
Positrons = S.Species(1., thermal_velocity, rho_zero, nppc, color='k')

delta_t = L/(2.*N_grid)

PIC = P1D.PIC_1D_rel(N_grid, L, delta_t, Species = [Electrons, Positrons], time_steps = 3000, folder='Relativistic_test' )
PIC.relativistic_two_Stream(gamma_zero)
PIC.Initiate_Relativistic_Plasma(snapshot = 100)

#"""
#"""
###############################################################################
#                                                                             #
#                       Relativistic variation of gamma                       #
#                                                                             #
###############################################################################

#input parameters for the 1D grid
N_grid = 64.  # input: Number of grid cells

# for instabilities to occur  L > np.sqrt(2)*np.pi*v_zero*(gamma_zero)**(3/2)/omega_p
L = 200.  # input: length of the domain

delta_t = L/(4.*N_grid)

charge_to_mass = 1.  # q_p/m_p
thermal_velocity = 0.1 # input
rho_zero = 1. # background charge density
nppc = 50. # input: Number of particles per cell

gamma = [2, 4, 6, 8, 10] # => L > [10.88, 34.41, 64.38, 99.7, 139.79]

#v_zero = sqrt(1- (1./gamma)**2)

for i in gamma:
    v_zero = np.sqrt(1. - (1./i)**2)
 
    #define the different particle species
    Electrons = S.Species(-1., thermal_velocity, rho_zero, nppc, color='red')
    Positrons = S.Species(1., thermal_velocity, rho_zero, nppc, color='k')    

    PIC = P1D.PIC_1D_rel(N_grid, L, delta_t, Species = [Electrons, Positrons], time_steps = 2000, folder='gamma_' + str(i))
    PIC.relativistic_two_Stream(i)
    PIC.Initiate_Relativistic_Plasma(snapshot = 10)
#"""