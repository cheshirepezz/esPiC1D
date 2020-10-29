import os
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

label_size = 20
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size


class PIC_1D_rel:
    
    def __init__(self, N_grid, Domain_Length, delta_t, Species, time_steps = 100, folder=None, rel=True):
        
        self.N_grid = float(N_grid)  # the amount of grid cells 
        self.L = float(Domain_Length)  # the length of the domain
        
        self.delta_x = Domain_Length/N_grid  # size of each of the cells
        self.x_grid = np.linspace(0., Domain_Length, num=N_grid+1)  # node positions of the grid        

        self.Species = Species  # the different particle species 
        
        self.nppc = 0.
        for i in self.Species:
            self.nppc += i.nppc  # the total amount of particles in each cell          
        self.delta_x_p = self.delta_x/self.nppc  # particle spatial step
        
        self.delta_t = float(delta_t)  # delta_t = delta_x/2 The Time step (has to satisfy thermal_velocity*Cell_size/time_step < 1)
        self.time_steps = time_steps
        self.time_array = np.arange(0, (self.time_steps + 1.)*self.delta_t, self.delta_t)
        
        self.output_folder = os.getcwd() + '/Output/' + str(folder) + '/'
#        if not os.path.exists(self.output_folder + 'Elec_E/'): os.makedirs(self.output_folder + 'Elec_E/')
#        if not os.path.exists(self.output_folder + 'Kin_E/'): os.makedirs(self.output_folder + 'Kin_E/')
        if not os.path.exists(self.output_folder + 'Phase_Space/'): os.makedirs(self.output_folder + 'Phase_Space/')
        
        self.E_grid = np.zeros(len(self.x_grid))
        self.J_grid = np.zeros(len(self.x_grid))
        
        # Values needed for two stream
        self.Elec_E = []
        self.Kin_E = []
        self.v_max = []
        
        #Values needed for Relativistic two stream
        self.gamma_max = []
        
        # initial position and velocities
        self.initialize_positions()   
        
        for i in self.Species:
            i.initialize_velocities(self.N_grid)
            i.set_charge(self.delta_x)
        
    def initialize_positions(self):  

        # determine where the particles should be placed
        x_pos = np.arange(self.delta_x_p/2, self.L, self.delta_x_p) 
            
        #determine the average distance between particles of the same species
        pp = np.array([self.nppc/(j.nppc) for j in self.Species])
        
        # create approximated index array
        for j in range(len(self.Species)):
            self.Species[j].index_array = list(np.arange(j, len(x_pos+1.), pp[j]))
            self.Species[j].index_array.append(np.infty)
            
        # find the closest approximated index and put the particle in its spot 
        for i in range(len(x_pos)):
            
            index = np.array([j.index_array[0] for j in self.Species])
            closest = (index-i).argmin()
            
            self.Species[closest].x_pos.append(x_pos[i])
            del self.Species[closest].index_array[0]
        
    def two_Stream(self, v_zero):
        
        self.v_zero = v_zero
     
        for i in self.Species:
            
            v_beam = [None]*len(i.v_pos)
            
            v_beam[::2] = [v + v_zero for v in i.v_pos[::2]]
            v_beam[1::2] = [v - v_zero for v in i.v_pos[1::2]]
    
            i.v_pos = v_beam
            
    def relativistic_two_Stream(self, gamma_zero):
        
        #gamma_zero = 1./np.sqrt(1 - v_zero**2)
        v_zero =  np.sqrt(1. - (1./gamma_zero)**2)
        
        for i in self.Species:          
            u_beam = [None]*len(i.v_pos)
            
            u_beam[::2] = [v/np.sqrt(1. - v**2) + v_zero*gamma_zero for v in i.v_pos[::2]]
            u_beam[1::2] = [v/np.sqrt(1. - v**2) - v_zero*gamma_zero for v in i.v_pos[1::2]]
    
            i.u_pos = u_beam
            
    def take_snapshot(self, current_time_step):

        if current_time_step == self.time_steps: # or !=0
        
            fig = plt.figure(figsize=(15, 10))
            plt.semilogy(self.time_array[0:len(self.Elec_E)], self.Elec_E)
            plt.xlabel('Time', fontsize=25)
            plt.ylabel(r'E$_{elec}$', fontsize=25)
            fig.tight_layout()
            plt.savefig(self.output_folder + 'Elec_E.png')
            plt.close()
            
            fig = plt.figure(figsize=(15, 10))
            plt.semilogy(self.time_array[0:len(self.Kin_E)], self.Kin_E)
            plt.ylim([np.min(self.Kin_E), np.max(self.Kin_E)])
            plt.xlabel('Time', fontsize=25)
            plt.ylabel(r'E$_{kin}$', fontsize=25)
            fig.tight_layout()
            plt.savefig(self.output_folder + 'Kin_E.png')
            plt.close()
            
            fig = plt.figure(figsize=(15, 10))
            plt.plot(self.time_array[0:len(self.v_max)], self.v_max)
            plt.xlabel('Time', fontsize=25)
            plt.ylabel(r'v$_{max}$', fontsize=25)
            fig.tight_layout()
            plt.savefig(self.output_folder + 'v_max.png')
            plt.close()
        
        fig = plt.figure(figsize=(15, 10))
        
        for l in range(len(self.Species)):
                                  
            plt.scatter(self.Species[l].x_pos, self.Species[l].v_pos, color=self.Species[l].color)
            plt.axvline(0., color='red')
            plt.axvline(self.L, color='red')
            plt.xlim([-0.1*self.L, 1.1*self.L])
            plt.ylim([-2.0*self.v_zero, 2.0*self.v_zero])
        
        plt.xlabel('x', fontsize=25)
        plt.ylabel('v', fontsize=25) 
        plt.tight_layout()                        
        plt.savefig(self.output_folder + 'Phase_Space/%0*d'%(5, current_time_step) + '.png') 
        plt.close()
        
    def take_rel_snapshot(self, current_time_step):

        if current_time_step == self.time_steps: # or  != 0
            
            
            fig = plt.figure(figsize=(15, 10))
            plt.semilogy(self.time_array[0:len(self.Elec_E)], self.Elec_E)
            plt.xlabel('Time', fontsize=25)
            plt.ylabel(r'E$_{elec}$', fontsize=25)
            fig.tight_layout()
            plt.savefig(self.output_folder + 'Elec_E.png')
            plt.close()
            
            fig = plt.figure(figsize=(15, 10))
            plt.semilogy(self.time_array[0:len(self.Kin_E)], self.Kin_E)
            plt.ylim([np.min(self.Kin_E), np.max(self.Kin_E)])
            plt.xlabel('Time', fontsize=25)
            plt.ylabel(r'E$_{kin}$', fontsize=25)
            fig.tight_layout()
            plt.savefig(self.output_folder + 'Kin_E.png')
            plt.close()
            
            fig = plt.figure(figsize=(15, 10))
            plt.plot(self.time_array[0:len(self.gamma_max)], self.gamma_max)
            plt.xlabel('Time', fontsize=25)
            plt.ylabel(r'$\gamma_{max}$', fontsize=25)
            fig.tight_layout()
            plt.savefig(self.output_folder + 'gamma_max.png')
            plt.close()
        
        fig = plt.figure(figsize=(15, 10))
        
        for l in range(len(self.Species)):
                       
            plt.scatter(self.Species[l].x_pos, self.Species[l].u_pos/np.sqrt(1+np.square(self.Species[l].u_pos)), color=self.Species[l].color)
            plt.axvline(0., color='red')
            plt.axvline(self.L, color='red')
            plt.xlim([-0.1*self.L, 1.1*self.L])
        
        plt.xlabel('x', fontsize=25)
        plt.ylabel('v', fontsize=25)
        fig.tight_layout()
        plt.savefig(self.output_folder + 'Phase_Space/%0*d'%(5, current_time_step) + '.png') 
        plt.close()
            
    def Initiate_Plasma(self, snapshot=0, Output_xv=0):
                         
        for i in range(0, self.time_steps + 1):
            
            t1 = time.time()
            
            for l in self.Species:
                    
                for j in range(len(l.x_pos)):
                        
                    # calculate the Electric field at the particle position with a b-spline
                    k = int(np.floor(l.x_pos[j]/self.delta_x))
                    E_pos = self.E_grid[k]*(1. - (l.x_pos[j]-self.x_grid[k])/self.delta_x) + self.E_grid[k+1]*(l.x_pos[j]-self.x_grid[k])/self.delta_x
                        
                    # update the positions and check boundaries
                    l.x_pos[j] += self.delta_t*l.v_pos[j]
                        
                    if l.x_pos[j] > self.L:
                        l.x_pos[j] -= self.L
                    if l.x_pos[j] < 0:
                        l.x_pos[j] += self.L
                        
                    # update the velocity  
                    l.v_pos[j] += l.charge_to_mass*E_pos*self.delta_t
                        
            self.J_grid = np.zeros(len(self.x_grid))
            
            for l in self.Species:
                    
                for j in range(len(l.x_pos)):
                        
                    #calculate the current
                    h = int(np.floor(l.x_pos[j]/self.delta_x))
                        
                    self.J_grid[h] += (1.-(l.x_pos[j]-self.x_grid[h])/self.delta_x)*l.charge_p*l.v_pos[j]/self.delta_x
                    self.J_grid[h+1] += (l.x_pos[j]-self.x_grid[h])*l.charge_p*l.v_pos[j]/self.delta_x**2
            
            # reflective boundary conditions
            self.J_grid[0] += self.J_grid[-1]
            self.J_grid[-1] = self.J_grid[0]
            
            self.E_grid -= self.delta_t*self.J_grid
            
            # Store Electric Energy
            self.Elec_E.append(np.sum(np.square(self.E_grid[:-1])*self.delta_x/2))
            
            # Store Kinetic Energy
            kin_energy = 0
            v_max = 0
            for l in self.Species:
                if v_max < np.amax(np.abs(l.v_pos)):
                    v_max = np.amax(np.abs(l.v_pos))

                kin_energy += np.sum(np.square(l.v_pos)*l.determine_mass()/2.)
                
            self.Kin_E.append(kin_energy)

            
            # Store Maximum Velocity
            self.v_max.append(v_max)
            
            if snapshot != 0 and snapshot <= self.time_steps and i%(self.time_steps/snapshot) == 0:
                self.take_snapshot(i)
            
            if Output_xv != 0 and i%(self.time_steps/Output_xv) == 0:
                for l in range(len(self.Species)):
                    with open(self.output_folder + 'txt_files/x_positions_' + str(l) + '.txt', "a") as output_file:                
                        np.savetxt(output_file, [self.Species[l].x_pos], delimiter=" \t", fmt="%.7E")
                        
                for l in range(len(self.Species)):    
                    with open(self.output_folder + 'txt_files/v_positions_' + str(l) + '.txt', "a") as output_file:    
                        np.savetxt(output_file, [self.Species[l].v_pos], delimiter=" \t", fmt="%.7E")
                    
            t2 = time.time()
            if i%10 == 0:
                print("(Iteration : {}) {}s (Estimated time left: {}s)".format(i, t2-t1, (t2-t1)*(self.time_steps - i)))          

        Header = [r'Time ($\omega_p^{-1}$)', r'E_{Elec}', r'E_{Kin}', r'v_{max}']
        DATA = np.column_stack((self.time_array, self.Elec_E, self.Kin_E, self.v_max))
        data = np.row_stack((Header, DATA))
        with file(self.output_folder + 'Parameters.txt', 'w') as output_file:
            np.savetxt(output_file, data, delimiter = " \t", fmt = "%s")
              

    def Initiate_Relativistic_Plasma(self, snapshot=0, Output_xv=0):
        
        """
        instead of using the velocity, we now use the momentum: u = gamma*v
        where gamma = 1./np.sqrt(1.-v**2) = np.sqrt(1. + u**2) 
        
        The formulas then become:
            
            du/dt = (Q/m) * E 
            
            dx/dt = u/gamma
            
            J = Sum_p  (u/gamma) * Q_p * W/delta_x    # the sum is over al the particles and W is the interpolation function 
            
        Also not that the initialization will change to:
            
            U_p = u_b + u_th = V_b*gamma_b + V_th*gamma_th
            
        And of course the Formula for the kinetic energy also changes:
            
            E_kin = Sum_p m_p * (gamma - 1.)
            
        and gamma_max is saved instead of v_max
        
        """
                      
        for i in range(0, self.time_steps + 1):
            
            t1 = time.time()
            
            for l in self.Species:
                    
                for j in range(len(l.x_pos)):
                        
                    # calculate the Electric field at the particle position with a b-spline
                    k = int(np.floor(l.x_pos[j]/self.delta_x))
                    E_pos = self.E_grid[k]*(1. - (l.x_pos[j]-self.x_grid[k])/self.delta_x) + self.E_grid[k+1]*(l.x_pos[j]-self.x_grid[k])/self.delta_x
                        
                    # update the positions and check boundaries
                    l.x_pos[j] += self.delta_t*l.u_pos[j]/np.sqrt(1. + l.u_pos[j]**2)
                    
                    if l.x_pos[j]> self.L:
                        l.x_pos[j] -= self.L
                    if l.x_pos[j] < 0:
                        l.x_pos[j] += self.L
                        
                    # update the velocity  
                    l.u_pos[j] += l.charge_to_mass*E_pos*self.delta_t
                        
            self.J_grid = np.zeros(len(self.x_grid))
            
            for l in self.Species:
                    
                for j in range(len(l.x_pos)):
                        
                    #calculate the current
                    h = int(np.floor(l.x_pos[j]/self.delta_x))
                    
                    gamma = np.sqrt(1. + l.u_pos[j]**2) 
                        
                    self.J_grid[h] += (1.-(l.x_pos[j]-self.x_grid[h])/self.delta_x)*l.charge_p*l.u_pos[j]/(gamma*self.delta_x)
                    self.J_grid[h+1] += (l.x_pos[j]-self.x_grid[h])*l.charge_p*l.u_pos[j]/(gamma*self.delta_x**2)
            
            # reflective boundary conditions
            self.J_grid[0] += self.J_grid[-1]
            self.J_grid[-1] = self.J_grid[0]
            
            self.E_grid -= self.delta_t*self.J_grid
            
            # Store Electric Energy
            self.Elec_E.append(np.sum(np.square(self.E_grid[:-1])*self.delta_x/2))
            
            # Store Kinetic Energy
            kin_energy = 0 
            gamma_max = 0
            for l in self.Species:
                
                if gamma_max < np.amax(np.sqrt(1. + np.square(l.u_pos))):
                    gamma_max = np.amax(np.sqrt(1. + np.square(l.u_pos)))
                
                for j in range(len(l.u_pos)):
                    kin_energy += np.sum(np.sqrt(np.square(l.u_pos[j]) + 1.)*l.determine_mass())
                    
            #print kin_energy +  np.sum(np.square(self.E_grid[:-1])*self.delta_x/2)   
            self.Kin_E.append(kin_energy)        
            # Store Maximum gamma
            self.gamma_max.append(gamma_max)
            
            if snapshot != 0 and snapshot <= self.time_steps and i%(self.time_steps/snapshot) == 0:
                self.take_rel_snapshot(i)
            
            if Output_xv != 0 and i%(self.time_steps/Output_xv) == 0:
                for l in range(len(self.Species)):
                    with open(self.output_folder + 'txt_files/x_positions_' + str(l) + '.txt', "a") as output_file:                
                        np.savetxt(output_file, [self.Species[l].x_pos], delimiter=" \t", fmt="%.7E")
                        
                for l in range(len(self.Species)):    
                    with open(self.output_folder + 'txt_files/v_positions_' + str(l) + '.txt', "a") as output_file:    
                        np.savetxt(output_file, [self.Species[l].v_pos], delimiter=" \t", fmt="%.7E")
                  
            t2 = time.time()
            if i%10 == 0:
                print("(Iteration : {}) {}s (Estimated time left: {}s)".format(i, t2-t1, (t2-t1)*(self.time_steps - i)))          
        
        Header = [r'Time ($\omega_p^{-1}$)', r'E_{Elec}', r'E_{Kin}', r'gamma_{max}']
        DATA = np.column_stack((self.time_array, self.Elec_E, self.Kin_E, self.gamma_max))
        data = np.row_stack((Header, DATA))
        with file(self.output_folder + 'Parameters.txt', 'w') as output_file:
            np.savetxt(output_file, data, delimiter = " \t", fmt = "%s")