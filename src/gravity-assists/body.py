from astropy import units as u
from astropy.constants import G, M_sun
import numpy as np
from skyfield import elementslib

class Body():
    def __init__(self, name, time, body_eph, radius, mass):
        self.name = name
        self.body_eph = body_eph
        self.radius = radius.to(u.m)
        self.mass = mass * u.kg

        self.time = time
        # Depends on time
        self.orbital_elements = elementslib.osculating_elements_of(self.get_body_position(self.time))
        self.orbital_period = self.orbital_elements.period_in_days
        self.semimajor_axis = self.orbital_elements.semi_major_axis
        self.position = self.calc_position() # [x, y]
        self.velocity = self.calc_velocity() # [x, y]
        self.r_soi = self.calc_r_soi() # au

    def __str__(self):
        return f"Body: {self.name}"
    
    def get_body_position(self, time):
        """
        ### Parameters
        - time: Skyfield time object

        ### Returns
        - position: Skyfield position object
        """
        return self.body_eph.at(time)

    def calc_position(self):
        pos = self.get_body_position(self.time).position.au
        return [pos[0] * u.au, pos[1] * u.au]

    def calc_velocity(self):
        """
        Calculates a planetary-style body's velocity given the orbital period
        ### Returns
        - [x_vel, y_vel]
        """
        angle = self.orbital_elements.true_anomaly.to(u.rad)
        velocity_mag = (((2 * np.pi * self.semimajor_axis.to(u.km).value) * u.km) / ((self.orbital_period * 24 * 60 * 60) * u.s)) # km/s
        x = velocity_mag * np.sin(angle)
        y = velocity_mag * np.cos(angle)
        return [x, y]
    
    def calc_r_soi(self):
        """
        Calculates the radius of the sphere of influence of a body
        Uses the Hill Radius formula
        ### Returns
        - r_soi: float au
        """
        r_soi = self.semimajor_axis.to(u.m) * (self.mass / (3 * (M_sun + self.mass)))**(1/3)
        return r_soi.to(u.au)
    
    def set_time(self, time):
        self.time = time
        self.orbital_elements = elementslib.osculating_elements_of(self.get_body_position(self.time))
        self.orbital_period = self.orbital_elements.period_in_days
        self.semimajor_axis = self.orbital_elements.semi_major_axis
        self.position = self.calc_position()
        self.velocity = self.calc_velocity()
        self.r_soi = self.calc_r_soi()