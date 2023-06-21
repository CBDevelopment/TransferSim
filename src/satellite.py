from astropy import units as u
from astropy.constants import G, M_sun
from planet import Planet
import numpy as np

class Satellite():
    def __init__(self, name, home_planet: Planet, transfer_planet: Planet, orbit_height, graph_ax):
        self.name = name
        self.home_planet = home_planet
        self.transfer_planet = transfer_planet
        self.orbit_height = orbit_height * u.km
        self.mu = G * M_sun
        self.delta_v_total = self.calc_delta_v_total()
        self.graph_ax = graph_ax

    def __str__(self):
        return f'{self.name} is a satellite of {self.home_planet.name}'
    
    def calc_escape_velocity(self):
        """
        ### Returns
        - The escape velocity of the home planet in m/s
        """
        return np.sqrt((2 * G * self.home_planet.mass) / self.orbit_height.to(u.m))
    
    # https://en.wikipedia.org/wiki/Hohmann_transfer_orbit#Calculation
    def calc_delta_v_1(self):
        """
        ### Returns
        - The delta v required to leave the orbit of the home planet in m/s
        """
        p1 = self.home_planet
        r1 = p1.calc_distance_from_sun(p1.right_ascension.to(u.rad)).to(u.m)
        p2 = self.transfer_planet
        r2 = p2.calc_distance_from_sun(p2.right_ascension.to(u.rad)).to(u.m)
        return np.sqrt(self.mu / r1) * (np.sqrt((2 * r2) / (r1 + r2)) - 1)
    
    def calc_delta_v_2(self):
        """
        ### Returns
        - The delta v required to enter the orbit of the transfer planet in m/s
        """
        p1 = self.home_planet
        r1 = p1.calc_distance_from_sun(p1.right_ascension.to(u.rad)).to(u.m)
        p2 = self.transfer_planet
        r2 = p2.calc_distance_from_sun(p2.right_ascension.to(u.rad)).to(u.m)
        return np.sqrt(self.mu / r2) * (1 - np.sqrt((2 * r1) / (r1 + r2)))
    
    def calc_delta_v_total(self):
        """
        ### Returns
        - The total delta v required to transfer from the home planet to the transfer planet in m/s
        """
        self.calc_delta_v_1() + self.calc_delta_v_2()
    
    def calc_semimajor_transfer_axis(self):
        """
        ### Returns
        - The semimajor axis of the transfer orbit in meters
        """
        p1 = self.home_planet
        p2 = self.transfer_planet
        return (p1.semimajor + p2.semimajor).to(u.m) / 2
    
    def calc_transfer_period(self):
        """
        ### Returns
        - The transfer period in seconds
        """
        k = (4 * np.pi**2) / self.mu
        return np.sqrt(k * self.calc_semimajor_transfer_axis()**3)