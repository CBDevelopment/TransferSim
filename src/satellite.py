from astropy import units as u
from astropy.constants import G, M_sun
from planet import Planet
import numpy as np

class Satellite():
    def __init__(self, name, home_planet: Planet, transfer_planet: Planet, tof, orbit_height, graph_ax):
        self.name = name
        self.home_planet = home_planet
        self.transfer_planet = transfer_planet
        self.tof = tof * u.d
        self.orbit_height = orbit_height * u.km
        self.mu = G * M_sun
        self.graph_ax = graph_ax

    def __str__(self):
        return f'{self.name} is a satellite of {self.home_planet.name}'
    
    def calc_escape_velocity(self):
        """
        ### Returns
        - The escape velocity of the home planet in m/s
        """
        return np.sqrt((2 * G * self.home_planet.mass) / self.orbit_height.to(u.m))