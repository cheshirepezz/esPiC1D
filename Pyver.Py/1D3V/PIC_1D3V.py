import os
import numpy as np
import matplotlib.pyplot as plt
import time
import matplotlib as mpl

label_size = 20
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size

class PIC_1D3V:
    
    def __init__(self, N_grid, Domain_Length, delta_t, Species, time_steps = 100, folder=None, rel=True):
        
        self.N_grid = float(N_grid)  # the amount of grid cells 
        self.L = float(Domain_Length)  # the length of the domain
        
        self.delta_x = Domain_Length/N_grid  # size of each of the cells
        self.x_grid_node = np.linspace(0., Domain_Length, num=N_grid+1)  # node positions of the grid 
        self.x_grid_center = np.linspace(self.delta_x/2, Domain_Length-self.delta_x/2, num=N_grid)  # node positions of the grid

        self.Species = Species  # the different particle species 
        
        self.nppc = 0.
        for i in self.Species:
            self.nppc += i.nppc  # the total amount of particles in each cell          
        self.delta_x_p = self.delta_x/self.nppc  # particle spatial step
        
        self.delta_t = float(delta_t)  # delta_t = delta_x/2 The Time step (has to satisfy thermal_velocity*Cell_size/time_step < 1)
        self.time_steps = time_steps
        self.time_array = np.arange(0, (self.time_steps + 0.99)*self.delta_t, self.delta_t)
        
        self.output_folder = os.getcwd() + '/Output/' + str(folder) + '/'
        if not os.path.exists(self.output_folder + 'Phase_Space/vx/'): os.makedirs(self.output_folder + 'Phase_Space/vx/')
        if not os.path.exists(self.output_folder + 'Phase_Space/vy/'): os.makedirs(self.output_folder + 'Phase_Space/vy/')
        if not os.path.exists(self.output_folder + 'Phase_Space/vz/'): os.makedirs(self.output_folder + 'Phase_Space/vz/')
        if not os.path.exists(self.output_folder + 'txt_files/'): os.makedirs(self.output_folder + 'txt_files/')
        
        self.E_grid_node = np.zeros((3, len(self.x_grid_node)))        
        self.B_grid_center = np.zeros((3, len(self.x_grid_center)))        
        self.B_grid_node = np.zeros((3, len(self.x_grid_node)))        
        self.J_grid_node = np.zeros((3, len(self.x_grid_node)))

        # Values needed for two stream
        self.Elec_E = []
        self.Mag_E = []
        self.v_max = [[],[],[]]
        
        # initial position and velocities
        self.initialize_positions()   
        
        for i in self.Species:
            i.initialize_velocity_vector(self.N_grid)
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
            
    def EM_two_Stream(self, v_zero):
        
        self.v_zero =v_zero
     
        for i in self.Species:
            
            v_beam = [None]*len(i.v_pos[1])
            
            v_beam[::2] = [v + v_zero for v in i.v_pos[1][::2]]
            v_beam[1::2] = [v - v_zero for v in i.v_pos[1][1::2]]
    
            i.v_pos[1] = v_beam
            print i.v_pos[1]
            
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
            plt.semilogy(self.time_array[0:len(self.Mag_E)], self.Mag_E)
            plt.ylim([np.min(self.Mag_E), np.max(self.Mag_E)])
            plt.xlabel('Time', fontsize=25)
            plt.ylabel(r'E$_{Mag}$', fontsize=25)
            fig.tight_layout()
            plt.savefig(self.output_folder + 'Mag_E.png')
            plt.close()
            
            fig = plt.figure(figsize=(15, 10))
            plt.plot(self.time_array[0:len(self.v_max[0])], self.v_max[0])
            plt.xlabel('Time', fontsize=25)
            plt.ylabel(r'v$_{max,x}$', fontsize=25)
            fig.tight_layout()
            plt.savefig(self.output_folder + 'vx_max.png')
            plt.close()
            
            fig = plt.figure(figsize=(15, 10))
            plt.plot(self.time_array[0:len(self.v_max[1])], self.v_max[1])
            plt.xlabel('Time', fontsize=25)
            plt.ylabel(r'v$_{max,y}$', fontsize=25)
            fig.tight_layout()
            plt.savefig(self.output_folder + 'vy_max.png')
            plt.close()
            
            fig = plt.figure(figsize=(15, 10))
            plt.plot(self.time_array[0:len(self.v_max[2])], self.v_max[2])
            plt.xlabel('Time', fontsize=25)
            plt.ylabel(r'v$_{max,z}$', fontsize=25)
            fig.tight_layout()
            plt.savefig(self.output_folder + 'vz_max.png')
            plt.close()
            
        fig = plt.figure(figsize=(15, 10))
        
        for l in range(len(self.Species)):
                                  
            plt.scatter(self.Species[l].x_pos, self.Species[l].v_pos[0], color=self.Species[l].color)
            plt.axvline(0., color='red')
            plt.axvline(self.L, color='red')
            plt.xlim([-0.1*self.L, 1.1*self.L])
            plt.ylim([-2.0*self.v_zero, 2.0*self.v_zero])
                              
        plt.xlabel('x', fontsize=25)
        plt.ylabel(r'v$_{x}$', fontsize=25) 
        plt.savefig(self.output_folder + 'Phase_Space/vx/' + str(current_time_step) + '.png') 
        plt.close()
        
        fig = plt.figure(figsize=(15, 10))
        
        for l in range(len(self.Species)):
                                  
            plt.scatter(self.Species[l].x_pos, self.Species[l].v_pos[1], color=self.Species[l].color)
            plt.axvline(0., color='red')
            plt.axvline(self.L, color='red')
            plt.xlim([-0.1*self.L, 1.1*self.L])
            plt.ylim([-2.0*self.v_zero, 2.0*self.v_zero])
              
        plt.xlabel('x', fontsize=25)    
        plt.ylabel(r'v$_{y}$', fontsize=25)                      
        plt.savefig(self.output_folder + 'Phase_Space/vy/' + str(current_time_step) + '.png') 
        plt.close()
        
        fig = plt.figure(figsize=(15, 10))
        
        for l in range(len(self.Species)):
                                  
            plt.scatter(self.Species[l].x_pos, self.Species[l].v_pos[2], color=self.Species[l].color)
            plt.axvline(0., color='red')
            plt.axvline(self.L, color='red')
            plt.xlim([-0.1*self.L, 1.1*self.L])
            plt.ylim([-2.0*self.v_zero, 2.0*self.v_zero])
         
        plt.xlabel('x', fontsize=25)
        plt.ylabel(r'v$_{z}$', fontsize=25)                           
        plt.savefig(self.output_folder + 'Phase_Space/vz/' + str(current_time_step) + '.png') 
        plt.close()
        
    def Initiate_Plasma(self, snapshot=0, Output_xv=0):
                         
        for i in range(0, self.time_steps + 1):
            
            t1 = time.time()
            
            for l in range(len(self.B_grid_node)):
                
                self.B_grid_node[l] = [(self.B_grid_center[l][g-1] + self.B_grid_center[l][g]) / 2 for g in range(len(self.B_grid_center[l]))] + [(self.B_grid_center[l][-1] + self.B_grid_center[l][0]) / 2]
            
            for l in self.Species:
                
                E_pos = np.zeros((3, len(l.x_pos)))
                B_pos = np.zeros((3, len(l.x_pos)))
                    
                for j in range(len(l.x_pos)):
                        
                    # calculate the Electric field at the particle position with a b-spline
                    k = int(np.floor(l.x_pos[j]/self.delta_x))                    
                    
                    for s in range(len(E_pos)):
                        E_pos[s][j] = self.E_grid_node[s][k]*(1. - (l.x_pos[j]-self.x_grid_node[k])/self.delta_x) + self.E_grid_node[s][k+1]*(l.x_pos[j]-self.x_grid_node[k])/self.delta_x
                        B_pos[s][j] = self.B_grid_node[s][k]*(1. - (l.x_pos[j]-self.x_grid_node[k])/self.delta_x) + self.B_grid_node[s][k+1]*(l.x_pos[j]-self.x_grid_node[k])/self.delta_x

                # update the positions and check boundaries
                l.x_pos += self.delta_t*l.v_pos[0]
                   
                index_right = np.argwhere(l.x_pos > self.L)
                l.x_pos[index_right] -= self.L
              
                index_left = np.argwhere(l.x_pos < 0)
                l.x_pos[index_left] += self.L
                        
                # update the velocity                 
                l.v_pos += l.charge_to_mass*self.delta_t*(E_pos + [l.v_pos[1]*B_pos[2] - l.v_pos[2]*B_pos[1], l.v_pos[2]*B_pos[0] - l.v_pos[0]*B_pos[2], l.v_pos[0]*B_pos[1] - l.v_pos[1]*B_pos[0]])
                                 
#                    
            #initialise the coponents of the current as zero 
            self.J_grid_node = np.zeros((3, len(self.x_grid_node)))
            
            for l in self.Species:
                    
                for j in range(len(l.x_pos)):
                        
                    #calculate the current
                    h = int(np.floor(l.x_pos[j]/self.delta_x)) 
                    
                    for s in range(len(self.J_grid_node)):                                       
                        self.J_grid_node[s][h] += (1. - (l.x_pos[j] - self.x_grid_node[h]) / self.delta_x) * l.charge_p * l.v_pos[s][j] / self.delta_x
                        self.J_grid_node[s][h+1] += (l.x_pos[j] - self.x_grid_node[h]) * l.charge_p * l.v_pos[s][j] / self.delta_x**2
            
            # reflective boundary conditions (for each component)
            for s in range(len(self.J_grid_node)): 
                self.J_grid_node[s][0] += self.J_grid_node[s][-1]
                self.J_grid_node[s][-1] = self.J_grid_node[s][0]
            
            
            B = np.zeros((3, len(self.E_grid_node[0])))
            # do this for By and Bz not nessacary fo Bx since that one is constant (and zero)
            for l in range(1, len(self.B_grid_node)):
                B[l] = [(self.B_grid_center[l][g] - self.B_grid_center[l][g-1]) / self.delta_x for g in range(len(self.B_grid_center[l]))] + [(self.B_grid_center[l][0] - self.B_grid_center[l][-1]) / self.delta_x]

            # Evolve the fields
            self.E_grid_node[0] -= self.delta_t*self.J_grid_node[0]
            self.E_grid_node[1] = self.E_grid_node[1] - self.delta_t*np.array(B[2]) - self.delta_t*self.J_grid_node[1]
            self.E_grid_node[2] = self.E_grid_node[2] + self.delta_t*np.array(B[1]) - self.delta_t*self.J_grid_node[2]
            
            
            # calculate the new magentic field        
            self.B_grid_center[1] = [self.B_grid_center[1][g] + self.delta_t*(self.E_grid_node[2][g+1]-self.E_grid_node[2][g])/self.delta_x for g in range(len(self.B_grid_center[0]))]
            self.B_grid_center[2] = [self.B_grid_center[2][g] - self.delta_t*(self.E_grid_node[1][g+1]-self.E_grid_node[1][g])/self.delta_x for g in range(len(self.B_grid_center[0]))]
                    
            # Store Electric Energy
            self.Elec_E.append(np.sum((np.square(self.E_grid_node[0][:-1]) + np.square(self.E_grid_node[1][:-1]) + np.square(self.E_grid_node[2][:-1]) )*self.delta_x/2))
            
            # Store Magnetic Energy
            self.Mag_E.append(np.sum((np.square(self.B_grid_node[0][:-1]) + np.square(self.B_grid_node[1][:-1]) + np.square(self.B_grid_node[2][:-1]) )*self.delta_x/2))
                  
            #store the maximum velocity in each direction                       
            for s in range(3):
                v_max = 0
                for l in self.Species:                                   
                    if v_max < np.amax(np.abs(l.v_pos[s])):
                        v_max = np.amax(np.abs(l.v_pos[s]))

                self.v_max[s].append(v_max)

            if snapshot != 0 and snapshot <= self.time_steps and i%(self.time_steps/snapshot) == 0:
                self.take_snapshot(i)
            
            if Output_xv != 0 and i%(self.time_steps/Output_xv) == 0:
                for l in range(len(self.Species)):
                    with open(self.output_folder + 'txt_files/x_positions_' + str(l) + '.txt', "a") as output_file:                
                        np.savetxt(output_file, [self.Species[l].x_pos], delimiter=" \t", fmt="%.7E")
                        
                for l in range(len(self.Species)):    
                    with open(self.output_folder + 'txt_files/vx_positions_' + str(l) + '.txt', "a") as output_file:    
                        np.savetxt(output_file, [self.Species[l].v_pos[0]], delimiter=" \t", fmt="%.7E")
                        
                for l in range(len(self.Species)):    
                    with open(self.output_folder + 'txt_files/vy_positions_' + str(l) + '.txt', "a") as output_file:    
                        np.savetxt(output_file, [self.Species[l].v_pos[1]], delimiter=" \t", fmt="%.7E")
                        
                for l in range(len(self.Species)):    
                    with open(self.output_folder + 'txt_files/vz_positions_' + str(l) + '.txt', "a") as output_file:    
                        np.savetxt(output_file, [self.Species[l].v_pos[2]], delimiter=" \t", fmt="%.7E")
                    
            t2 = time.time()
            if i%10 == 0:
                print("(Iteration : {}) {}s (Estimated time left: {}s)".format(i, t2-t1, (t2-t1)*(self.time_steps - i)))          
                
        Header = [r'Time ($\omega_p^{-1}$)', r'E_{Elec}', r'E_{Mag}', r'vx_{max}', r'vy_{max}', r'vz_{max}']
        DATA = np.column_stack((self.time_array, self.Elec_E, self.Mag_E, self.v_max[0], self.v_max[1], self.v_max[2]))
        data = np.row_stack((Header, DATA))
        with file(self.output_folder + 'Parameters.txt', 'w') as output_file:
            np.savetxt(output_file, data, delimiter = " \t", fmt = "%s")
