from astropy.constants import G, M_sun
import astropy.units as u
import numpy as np

class Planet():
    def __init__(self, semimajor, semiminor, orbital_inclination):
        self.semimajor = semimajor * u.AU
        self.semiminor = semiminor * u.AU
        self.orbital_inclination = orbital_inclination * u.deg
        self.eccentricity = self.calc_eccentricity()
        self.perihelion = self.calc_perihelion()
        self.aphelion = self.calc_aphelion()
        self.orbital_period = self.calc_orbital_period()

    def calc_eccentricity(self):
        return (1 - self.semiminor**2 / self.semimajor**2)**0.5 * u.dimensionless_unscaled
    
    def calc_perihelion(self):
        """
        Calculates the planet's position closest to the sun.
        ### Returns
        - perihelion: float, in AU
        """
        return self.semimajor * (1 - self.eccentricity)

    def calc_aphelion(self):
        """
        Calculates the planet's position furthest from the sun.
        ### Returns
        - aphelion: float, in AU
        """
        return self.semimajor * (1 + self.eccentricity)
    
    def calc_orbital_period(self):
        """
        ### Returns
        - orbital_period: float, in seconds
        """
        return np.sqrt(((4 * np.pi**2) / (G * M_sun)) * self.semimajor.to(u.m)**3)