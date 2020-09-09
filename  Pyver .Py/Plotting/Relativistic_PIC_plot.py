import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

label_size = 20
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size

#"""
###############################################################################
#                                                                             #
#                  relativistic two Stream: Vary Parameters                   #
#                                                                             #
###############################################################################

gamma_zero = [2.0, 4.0, 6.0]

L=200.
N_grid = [32., 48., 64., 96.] #, 128.
nppc =  [10., 50., 100., 200.] 
t_factor = [1./2., 1./4., 1./8.]

fig, axes = plt.subplots(3, 3, sharex=True, sharey=True, figsize=(20,20))
axes[0,0].set_xlim([0, 440])
axes[0,0].set_ylim([10**-5.5, 10**2.5])

folder=os.getcwd() + '/Electro_Static/Variation of Parameters for gamma_zero/'

for j in range(len(gamma_zero)):
    
    #variation of delta_t
    for k in t_factor:
        delta_t = L*k/64.
        
        Input_folder=folder + 'N_grid_64.0_nppc_50.0_delta_t_' + str(delta_t) + '_gamma_zero_' + str(gamma_zero[j])
        
        Time = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[0]
        Energy = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[1]

        axes[j,0].semilogy(Time, Energy, linewidth=2.0, label=r'$\Delta$t: $\Delta$x/' + str(int(1./k) ) )
        axes[j,0].legend(fontsize = 17, loc='lower right')
        axes[j,0].grid(linestyle='--', linewidth='1.', color='black')
                
    #variations of N_grid
    if j == 2:
        N_grid = [32., 48., 64.]
        
    for k in N_grid:
        delta_t = L/(4.*64)
        Input_folder=folder + 'N_grid_' + str(k) + '_nppc_50.0_delta_t_' + str(delta_t) + '_gamma_zero_' + str(gamma_zero[j])
        
        Time = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[0]
        Energy = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[1]
        
        axes[j,1].semilogy(Time, Energy, linewidth=2.0, label=r'$N_{grid}$: ' + str(k))
        axes[j,1].legend(fontsize = 17, loc='lower right')
        axes[j,1].grid(linestyle='--', linewidth='1.', color='black') 
        
    #variations of nppc
    for k in nppc:
        delta_t = L/(4.*64)
        Input_folder=folder + 'N_grid_64.0_nppc_' + str(k) + '_delta_t_' + str(delta_t)  + '_gamma_zero_' + str(gamma_zero[j])
        
        Time = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[0]
        Energy = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[1]

        axes[j,2].semilogy(Time, Energy, linewidth=2.0, label='nppc: ' + str(k))
        axes[j,2].legend(fontsize = 17, loc='lower right')
        axes[j,2].grid(linestyle='--', linewidth='1.0', color='black')

axes[0,0].set_ylabel(r'E$_{Elec}$, $\gamma_0$ = 2', fontsize = 25)
axes[1,0].set_ylabel(r'E$_{Elec}$, $\gamma_0$ = 4', fontsize = 25)
axes[2,0].set_ylabel(r'E$_{Elec}$, $\gamma_0$ = 6', fontsize = 25)

axes[2,0].set_xlabel(r'Time [$\omega_{p}^{-1}$]', fontsize = 25)
axes[2,1].set_xlabel(r'Time [$\omega_{p}^{-1}$]', fontsize = 25)
axes[2,2].set_xlabel(r'Time [$\omega_{p}^{-1}$]', fontsize = 25)

axes[0,0].set_title(r"Variation of $\Delta$t with" + '\n' + r"L=200, $N_{grid}$=64," + '\n' + r"nppc=50 and $v_{th}$=0.1", fontsize=30)
axes[0,1].set_title(r"Variation of $N_{grid}$ with" + '\n' +  "L=200, nppc=50," + '\n' + r"$\Delta$t=$\Delta$x/4 and $v_{th}$=0.1", fontsize=30)
axes[0,2].set_title(r"Variation of nppc with" + '\n' + r"L=200, $N_{grid}$=64," + '\n' + r"$\Delta$t=$\Delta$x/4 and $v_{th}$=0.1", fontsize=30)

fig.tight_layout()
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)
plt.savefig(os.getcwd() + '/Electro_Static/Parameter_variations_gamma.png')
plt.close(fig)

#"""
#"""
###############################################################################
#                                                                             #
#  2-Stream: Variations of gamma_0: N_grid=64, nppc=50 and delta_t=delta_x/4  #
#                                                                             #
###############################################################################

gamma = [2, 4, 6] #, 8, 10]
Time = [[None] for k in range(len(gamma))]
Energy = [[None] for k in range(len(gamma))]
Output_folder = os.getcwd() + '/Electro_Static/'


fig = plt.figure(figsize=(10,7))
plt.xlabel(r'Time [$\omega_{p}^{-1}$]', fontsize = 20) 
plt.ylabel(r'E$_{Elec}$', fontsize = 20) 
plt.xlim([0, 600])
plt.ylim([10**-5.5, 10**2.5]) 

for i in range(len(gamma)):
    
    Input_folder=os.getcwd() + '/Electro_Static/Gamma Variation/gamma_' + str(gamma[i])
            
    Time[i] = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[0]
    Energy[i] = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[1]  

    plt.semilogy(Time[i], Energy[i], linewidth=2.0, label=r'$\gamma_0$: ' +str(gamma[i]))

        
plt.title(r"Study of $\gamma_0$",  fontsize=25)        
plt.legend(loc='lower right')
plt.grid(linestyle='--', linewidth='1.', color='black')
fig.tight_layout()
plt.savefig(Output_folder + 'gamma_zero_variation.png')      
        
plt.close(fig)
#"""
#"""
###############################################################################
#                                                                             #
#  2-Stream: Variations of gamma_0: N_grid=64, nppc=50 and delta_t=delta_x/4  #
#                                                                             #
###############################################################################

gamma = [2, 4, 6, 8, 10]
Time = [[None] for k in range(len(gamma))]
Energy = [[None] for k in range(len(gamma))]
Output_folder = os.getcwd() + '/Electro_Static/'


fig = plt.figure(figsize=(10,7))
plt.xlabel(r'Time [$\omega_{p}^{-1}$]', fontsize = 20) 
plt.ylabel(r'E$_{Elec}$', fontsize = 20) 
plt.xlim([0, 800])
plt.ylim([10**-5.5, 10**2.5]) 

for i in range(len(gamma)):
    
    Input_folder=os.getcwd() + '/Electro_Static/Gamma Variation/gamma_' + str(gamma[i])
            
    Time[i] = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[0]
    Energy[i] = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[1]  

    plt.semilogy(Time[i], Energy[i], linewidth=2.0, label=r'$\gamma_0$: ' +str(gamma[i]))

        
plt.title(r"Study of $\gamma_0$",  fontsize=25)        
plt.legend(loc='lower right')
plt.grid(linestyle='--', linewidth='1.', color='black')
fig.tight_layout()
plt.savefig(Output_folder + 'gamma_zero_variation2.png')      
        
plt.close(fig)
#"""