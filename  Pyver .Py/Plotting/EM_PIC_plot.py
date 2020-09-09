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
#                         two Stream: Vary Parameters                         #
#                                                                             #
###############################################################################

v_zero = [0.01, 0.05, 0.1]

L=10.
N_grid = [64., 128., 256.]
nppc =  [10., 50., 100., 200.] #[10., 50., 100.] 
t_factor = [1./2., 1./4., 1./8.]

fig, axes = plt.subplots(3, 3, sharex='row', sharey=True, figsize=(20,20))
axes[0,0].set_xlim([0, 750])
axes[1,0].set_xlim([0, 200])
axes[2,0].set_xlim([0, 200])
axes[0,0].set_ylim([10**-8.5, 10**-1.5])

folder=os.getcwd() + '/Electro_Magnetic/Variation of Parameters for V_zero/'

for j in range(len(v_zero)):
    
    #variation of delta_t
    for k in t_factor:
        delta_t = 10.*k/128.
        
        Input_folder=folder + 'N_grid_128.0_nppc_50.0_delta_t_' + str(delta_t) + '_v_zero_' + str(v_zero[j])
        
        Time = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[0]
        Energy = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[2]

        axes[j,0].semilogy(Time, Energy, linewidth=2.0, label=r'$\Delta$t: $\Delta$x/' + str(int(1./k) ))
        axes[j,0].legend(fontsize=17, loc='lower right')
        axes[j,0].grid(linestyle='--', linewidth='1.', color='black')
                
    #variations of N_grid
    for k in N_grid:
        delta_t = 10./(4.*128)
        Input_folder=folder + 'N_grid_' + str(k) + '_nppc_50.0_delta_t_' + str(delta_t) + '_v_zero_' + str(v_zero[j])
        
        Time = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[0]
        Energy = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[2]

        axes[j,1].semilogy(Time, Energy, linewidth=2.0, label=r'$N_{grid}$: ' + str(k))
        axes[j,1].legend(fontsize=17, loc='lower right')
        axes[j,1].grid(linestyle='--', linewidth='1.', color='black')        
      
   
    if j != 0:
        nppc = [10., 50., 100.]            
    #variations of nppc
    for k in nppc:
        
        delta_t = 10./(4.*128)
        Input_folder=folder + 'N_grid_128.0_nppc_' + str(k) + '_delta_t_' + str(delta_t)  + '_v_zero_' + str(v_zero[j])
        
        Time = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[0]
        Energy = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[2]

        axes[j,2].semilogy(Time, Energy, linewidth=2.0, label='nppc: ' + str(k))
        axes[j,2].legend(fontsize=17, loc='lower right')
        axes[j,2].grid(linestyle='--', linewidth='1.0', color='black')

axes[0,0].set_ylabel(r'$E_{Mag}$, v$_0$ = 0.01', fontsize = 25)
axes[1,0].set_ylabel(r'$E_{Mag}$, v$_0$ = 0.05', fontsize = 25)
axes[2,0].set_ylabel(r'$E_{Mag}$, v$_0$ = 0.1', fontsize = 25)

axes[2,0].set_xlabel(r'Time [$\omega_{p}^{-1}$]', fontsize = 25)
axes[2,1].set_xlabel(r'Time [$\omega_{p}^{-1}$]', fontsize = 25)
axes[2,2].set_xlabel(r'Time [$\omega_{p}^{-1}$]', fontsize = 25)

axes[0,0].set_title(r"Variation of $\Delta$t with" + '\n' + r"L=10, $N_{grid}$=128," + '\n' +  r"nppc=50 and $v_{th}$=$v_b$/10", fontsize=30)
axes[0,1].set_title(r"Variation of $N_{grid}$ with" + '\n' +  "L= 10, nppc=50," + '\n' +  r"$\Delta$t=$\Delta$x/4 and $v_{th}$=$v_b$/10", fontsize=30)
axes[0,2].set_title(r"Variation of nppc with" + '\n' + r"L=10, $N_{grid}$=128," + '\n' +  r"$\Delta$t=$\Delta$x/4 and $v_{th}$=$v_b$/10", fontsize=30)

fig.tight_layout()
#fig.subplots_adjust(hspace=0)
#fig.subplots_adjust(wspace=0)
plt.savefig(os.getcwd() + '/Electro_Magnetic/Parameter_variations_velocity.png')
plt.close(fig)

#"""
#"""
###############################################################################
#                                                                             #
#   2-Stream: Variations of v_0: N_grid=128, nppc=500 and delta_t=delta_x/4   #
#                                                                             #
###############################################################################

v = [0.01, 0.05, 0.1] 
Time = [[None] for k in range(len(v))]
Energy = [[None] for k in range(len(v))]
Output_folder = os.getcwd() + '/Electro_Magnetic/'


fig = plt.figure(figsize=(10,7))
plt.xlabel(r'Time [$\omega_{p}^{-1}$]', fontsize = 20) 
plt.ylabel(r'$E_{Mag}$', fontsize = 20) 
#plt.xlim([0, 40])
plt.ylim([10**-8.5, 10**-1.5])
 
for i in range(len(v)):
    
    Input_folder=os.getcwd() + '/Electro_Magnetic/Velocity Variation/v_zero_' + str(v[i])
            
    Time[i] = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[0]
    Energy[i] = np.loadtxt(os.path.normpath(Input_folder + '/Parameters.txt'), skiprows=1, unpack=True)[2]  

    plt.semilogy(Time[i], Energy[i], linewidth=2.0, label=r'$v_b$: ' +str(v[i]))

        
plt.title(r"Study of $v_b$",  fontsize=25)        
plt.legend(loc='lower right')
plt.grid(linestyle='--', linewidth='1.', color='black')
fig.tight_layout()
plt.savefig(Output_folder + 'V_zero_variation.png')      
        
plt.close(fig)
#"""