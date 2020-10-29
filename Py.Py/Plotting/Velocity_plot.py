import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

label_size = 20
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size

#"""
###############################################################################
#                                                                             #
#   2-Stream: Variations of v_0: N_grid=512, nppc=100 and delta_t=delta_x/4   #
#                                                                             #
###############################################################################

v = [0.01, 0.05, 0.1] 
Output_folder = os.getcwd() + '/Electro_Static/'

fig = plt.figure(figsize=(10,7))
plt.xlabel(r'Time [$\omega_{p}^{-1}$]', fontsize = 20) 
plt.ylabel(r'v$_{max}$', fontsize = 20) 
plt.xlim([0, 40])
#plt.ylim([10**-8.5, 10**-1.5])
 
for i in range(len(v)):
    
    Input_folder=os.getcwd() + '/Electro_Static/Velocity Variation/v_zero_' + str(v[i])
            
    Time = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[0]
    Velocity = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[3]  
    
    plt.semilogy(Time, Velocity, linewidth=2.0, label=r'$v_b$: ' +str(v[i]))
    
plt.title(r"Velocity Study of $v_b$",  fontsize=25)        
plt.legend(loc='lower right')
plt.grid(linestyle='--', linewidth='1.', color='black')
fig.tight_layout()
plt.savefig(Output_folder + 'V_zero_velocity_study.png')    
plt.close(fig)
#"""
#"""
###############################################################################
#                                                                             #
#  2-Stream: Variations of gamma_0: N_grid=64, nppc=50 and delta_t=delta_x/4  #
#                                                                             #
###############################################################################

gamma = [2, 4, 6] #, 8, 10]
Output_folder = os.getcwd() + '/Electro_Static/'

fig = plt.figure(figsize=(10,7))
plt.xlabel(r'Time [$\omega_{p}^{-1}$]', fontsize = 20) 
plt.ylabel(r'$\gamma_{max}$', fontsize = 20) 
plt.xlim([0, 440])
#plt.ylim([10**-5.5, 10**2.5]) 

for i in range(len(gamma)):
    
    Input_folder=os.getcwd() + '/Electro_Static/Gamma Variation/gamma_' + str(gamma[i])
            
    Time = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[0]
    Velocity = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[3] 
    
    plt.semilogy(Time, Velocity, linewidth=2.0, label=r'$\gamma_0$: ' +str(gamma[i]))

        
plt.title(r"Velocity Study of $\gamma_0$",  fontsize=25)        
plt.legend(loc='lower right')
plt.grid(linestyle='--', linewidth='1.', color='black')
fig.tight_layout()
plt.savefig(Output_folder + 'gamma_zero_velocity_study.png')      
        
plt.close(fig)
#"""
#"""
###############################################################################
#                                                                             #
#   2-Stream: Variations of v_0: N_grid=128, nppc=500 and delta_t=delta_x/4   #
#                                                                             #
###############################################################################

v = [0.01, 0.05, 0.1] 
Output_folder = os.getcwd() + '/Electro_Magnetic/'

Time_begin = [150, 30, 10]
Time_end = [370, 90, 40]

fig = plt.figure(figsize=(10,7))
plt.xlabel(r'Time [$\omega_{p}^{-1}$]', fontsize = 20) 
plt.ylabel(r'$v_{max}$', fontsize = 20) 
#plt.xlim([0, 400])
#plt.ylim([10**-8.5, 10**-1.5])
 
for i in range(len(v)):
    
    Input_folder=os.getcwd() + '/Electro_Magnetic/Velocity Variation/v_zero_' + str(v[i])
            
    Time = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[0]
    Velocity_x = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[3]
    Velocity_y = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[4]
    Velocity_z = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[5]

    
    plt.semilogy(Time, Velocity_x, linewidth=2.0, label=r'$v_b$: ' + str(v[i]))
    
plt.title(r"Velocity $v_x$ Study of $v_b$",  fontsize=25)        
plt.legend(loc='lower right')
plt.grid(linestyle='--', linewidth='1.', color='black')
fig.tight_layout()
plt.savefig(Output_folder + 'V_zero_velocityx_study.png')      
        
plt.close(fig)
#"""
#"""
###############################################################################
#                                                                             #
#   2-Stream: Variations of v_0: N_grid=128, nppc=500 and delta_t=delta_x/4   #
#                                                                             #
###############################################################################

v = [0.01, 0.05, 0.1] 
Output_folder = os.getcwd() + '/Electro_Magnetic/'

Time_begin = [150, 30, 10]
Time_end = [370, 90, 40]

fig = plt.figure(figsize=(10,7))
plt.xlabel(r'Time [$\omega_{p}^{-1}$]', fontsize = 20) 
plt.ylabel(r'$v_{max}$', fontsize = 20) 
#plt.xlim([0, 40])
#plt.ylim([10**-8.5, 10**-1.5])
 
for i in range(len(v)):
    
    Input_folder=os.getcwd() + '/Electro_Magnetic/Velocity Variation/v_zero_' + str(v[i])
            
    Time = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[0]
    Velocity_x = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[3]
    Velocity_y = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[4]
    Velocity_z = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[5]

    
    plt.semilogy(Time, Velocity_y, linewidth=2.0, label=r'$v_b$: ' + str(v[i]))
    
plt.title(r"Velocity $v_y$ Study of $v_b$",  fontsize=25)        
plt.legend(loc='lower right')
plt.grid(linestyle='--', linewidth='1.', color='black')
fig.tight_layout()
plt.savefig(Output_folder + 'V_zero_velocityy_study.png')      
        
plt.close(fig)
#"""
#"""
###############################################################################
#                                                                             #
#   2-Stream: Variations of v_0: N_grid=128, nppc=500 and delta_t=delta_x/4   #
#                                                                             #
###############################################################################

v = [0.01, 0.05, 0.1] 
Output_folder = os.getcwd() + '/Electro_Magnetic/'

Time_begin = [150, 30, 10]
Time_end = [370, 90, 40]

fig = plt.figure(figsize=(10,7))
plt.xlabel(r'Time [$\omega_{p}^{-1}$]', fontsize = 20) 
plt.ylabel(r'$v_{max}$', fontsize = 20) 
#plt.xlim([0, 40])
#plt.ylim([10**-8.5, 10**-1.5])
 
for i in range(len(v)):
    
    Input_folder=os.getcwd() + '/Electro_Magnetic/Velocity Variation/v_zero_' + str(v[i])
            
    Time = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[0]
    Velocity_x = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[3]
    Velocity_y = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[4]
    Velocity_z = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[5]

    
    plt.semilogy(Time, Velocity_z, linewidth=2.0, label=r'$v_b$: ' + str(v[i]))
    
plt.title(r"Velocity $v_z$ Study of $v_b$",  fontsize=25)        
plt.legend(loc='lower right')
plt.grid(linestyle='--', linewidth='1.', color='black')
fig.tight_layout()
plt.savefig(Output_folder + 'V_zero_velocityz_study.png')      
        
plt.close(fig)
#"""
