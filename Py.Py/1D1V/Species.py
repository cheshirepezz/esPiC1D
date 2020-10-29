import numpy as np

class Species:

    def __init__(self, charge_to_mass, thermal_velocity, rho_zero, nppc, name='particle', color='blue'):
        
        self.charge_to_mass = charge_to_mass
        self.thermal_velocity = thermal_velocity
        self.rho_zero = rho_zero  # background charge
        self.nppc = nppc  # number of particles per cell
        self.name = name
        self.color = color
        
        self.x_pos = []
        self.v_pos = []
        self.u_pos = []

    def set_charge(self, delta_x):
        self.charge_p = self.rho_zero*delta_x*np.sign(self.charge_to_mass)/self.nppc #rho_zero*delta_x/nppc should prolly be in species class
        
    def initialize_velocities(self, N_grid):
        self.v_pos = (self.thermal_velocity*np.random.standard_normal(int(N_grid*self.nppc))).tolist() # initial particle velocities
        
    def determine_mass(self):
        return np.abs(self.charge_p/(self.charge_to_mass)) #np.sign(self.charge_to_mass)*delta_x*self.rho_zero/(self.nppc*self.charge_to_mass)
    