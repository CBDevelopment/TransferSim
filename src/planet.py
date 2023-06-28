from astropy.constants import G, M_sun
import astropy.units as u
import numpy as np
import copy

class Planet():
    def __init__(self, name: str, semimajor: float, semiminor: float, orbital_inclination: float, radius: float, right_ascension: float, mass):
        self.name = name
        self.semimajor = semimajor * u.AU
        self.semiminor = semiminor * u.AU
        self.orbital_inclination = (90 - (90 + orbital_inclination)) * u.deg
        self.right_ascension = right_ascension # degrees
        self.radius = radius * u.km
        self.eccentricity = self.calc_eccentricity()
        self.perihelion = self.calc_perihelion()
        self.aphelion = self.calc_aphelion()
        self.orbital_period = self.calc_orbital_period()
        self.velocity = self.calc_velocity(self.right_ascension)
        self.mass = mass * u.kg

    def __str__(self):
        return f"Planet: {self.name}\nRA: {self.right_ascension:.3f},\nVelocity: {self.velocity:.3f}\n"

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
    
    def calc_distance_from_sun(self, theta):
        """
        ### Parameters
        - theta: float, in radians

        ### Returns
        - distance: float, in AU
        """
        return self.semimajor * (1 - self.eccentricity**2) / (1 + self.eccentricity * np.cos(theta))

    def calc_velocity(self, theta):
        """
        ### Parameters
        - theta: float, in radians

        ### Returns
        - speed: float, in m / s
        """
        return np.sqrt((G * M_sun * self.semimajor.to(u.m) * (1 - self.eccentricity**2)) / (self.calc_distance_from_sun(theta).to(u.m)**2))
    
    def calc_position_after_days(self, days_past):
        velocity = copy.copy(self.velocity)
        angular_position = copy.copy(self.right_ascension)
        for j in range(int(days_past)):
            # planet velocity is in m/s
            # planet right ascension is in degrees
            velocity = self.calc_velocity(angular_position.to(u.rad).value)
            angular_velocity = velocity / self.calc_distance_from_sun(angular_position.to(u.rad)).to(u.m) * u.rad
            # angular velocity is in rad/s
            seconds_in_day = 86400 * u.s
            if days_past > 0:
                angular_position += (angular_velocity * seconds_in_day).to(u.deg)
            else:
                angular_position -= (angular_velocity * seconds_in_day).to(u.deg)
            if angular_position >= 360 * u.deg:
                angular_position -= 360 * u.deg
        return angular_position