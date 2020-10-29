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
omega = 1.0/2.0 
b = [8.4, 7.0, 6.55]

Output_folder = os.getcwd() + '/Electro_Static/'

fig = plt.figure(figsize=(10,7))
plt.xlabel(r'Time [$\omega_{p}^{-1}$]', fontsize = 20) 
plt.ylabel(r'E$_{Elec}$ []', fontsize = 20) 
plt.xlim([0, 40])
plt.ylim([10**-8.5, 10**-1.5])

for i in range(len(v)):
    
    Input_folder=os.getcwd() + '/Electro_Static/Velocity Variation/v_zero_' + str(v[i])

    Time = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[0]
    Energy = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[1]
    
    x_theory = Time[750:3000]
    y_theory = np.exp(2*omega*x_theory/np.sqrt(2))/(10**b[i]) #Energy[i][Indices]

    plt.semilogy(Time, Energy, linewidth=2.0, label=r'$v_b$: ' + str(v[i]))   
    plt.semilogy(x_theory, y_theory, linewidth=3.0, color='black')

        
plt.title(r"Study of $v_b$",  fontsize=25)        
plt.legend(loc='lower right')
plt.grid(linestyle='--', linewidth='1.', color='black')
fig.tight_layout()
plt.savefig(Output_folder + 'V_zero_variation_theory.png')    
plt.close(fig)
#"""
#"""
###############################################################################
#                                                                             #
#  2-Stream: Variations of gamma_0: N_grid=64, nppc=50 and delta_t=delta_x/4  #
#                                                                             #
###############################################################################

gamma = [2, 4, 6] #, 8, 10]
omega = 1./(2.0*np.power(gamma, 3.0/2.0))
b = [3.2, 4.55, 5.42] #, 5.6, 6.2] 

begin = [20, 70, 200] #, 150, 200]
end = [50, 190, 400] #, 600, 800]

Output_folder = os.getcwd() + '/Electro_Static/'

fig = plt.figure(figsize=(10,7))
plt.xlabel(r'Time [$\omega_{p}^{-1}$]', fontsize = 20) 
plt.ylabel(r'E$_{Elec}$', fontsize = 20) 
plt.xlim([0, 600])
plt.ylim([10**-5.5, 10**2.5]) 

for i in range(len(gamma)):
    
    Input_folder=os.getcwd() + '/Electro_Static/Gamma Variation/gamma_' + str(gamma[i])
            
    Time = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[0]
    Energy = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[1] 
    
    x_theory = Time[begin[i]:end[i]]
    y_theory = np.exp(2*omega[i]*x_theory/np.sqrt(2))/(10**b[i])
    
    plt.semilogy(Time, Energy, linewidth=2.0, label=r'$\gamma_0$: ' +str(gamma[i]))
    plt.semilogy(x_theory, y_theory, linewidth=3.0, color='black')

        
plt.title(r"Study of $\gamma_0$",  fontsize=25)        
plt.legend(loc='lower right')
plt.grid(linestyle='--', linewidth='1.', color='black')
fig.tight_layout()
plt.savefig(Output_folder + 'gamma_zero_variation_theory.png')      
        
plt.close(fig)
#"""
sys.exit
#"""
###############################################################################
#                                                                             #
#   2-Stream: Variations of v_0: N_grid=128, nppc=500 and delta_t=delta_x/4   #
#                                                                             #
###############################################################################

v = [0.01, 0.05, 0.1] 
omega = [np.sqrt(2)*i for i in v] 
b = [7.9, 6.83, 6.18]

begin = [8000, 1000, 300]
end = [20000, 4500, 2300]

Output_folder = os.getcwd() + '/Electro_Magnetic/'

fig = plt.figure(figsize=(10,7))
plt.xlabel(r'Time [$\omega_{p}^{-1}$]', fontsize = 20) 
plt.ylabel(r'$E_{Mag}$', fontsize = 20) 
#plt.xlim([0, 40])
plt.ylim([10**-8.5, 10**-1.5])
 
for i in range(len(v)):
    
    Input_folder=os.getcwd() + '/Electro_Magnetic/Velocity Variation/v_zero_' + str(v[i])
            
    Time = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[0]
    Energy = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[2]  
    
    x_theory = Time[begin[i]:end[i]]
    y_theory = np.exp(2*omega[i]*x_theory/np.sqrt(2))/(10**b[i])
    
    plt.semilogy(Time, Energy, linewidth=2.0, label=r'$v_b$: ' + str(v[i]))
    plt.semilogy(x_theory, y_theory, linewidth=3.0, color='black')
    
plt.title(r"Study of $v_b$",  fontsize=25)        
plt.legend(loc='lower right')
plt.grid(linestyle='--', linewidth='1.', color='black')
fig.tight_layout()
plt.savefig(Output_folder + 'V_zero_variation_theory.png')      
        
plt.close(fig)
#"""
