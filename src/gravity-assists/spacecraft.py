from astropy import units as u
from astropy.constants import G
import numpy as np
from body import Body

class Spacecraft():
    def __init__(self, name, mass, position_vec, velocity_vec, bodies: list[Body]):
        self.name = name
        self.bodies = bodies
        self.mass = mass * u.kg
        # Depends on time
        self.velocity = velocity_vec # [x, y] km / s
        self.position = position_vec # [x, y] au
        self.acceleration = self.calc_acceleration() # [x, y] m/s^2
        self.past_positions = [[self.position[0].value, self.position[1].value]]
        
    
    def calc_acceleration(self):
        """
        Calculates the acceleration of the spacecraft due to gravity from all bodies
        ### Returns
        - [x_accel, y_accel]
        """
        softening = 0.0001 * u.m
        acceleration = [0, 0]
        for body in self.bodies:
            distance = np.sqrt((self.position[0].to(u.m) - body.position[0].to(u.m))**2 + (self.position[1].to(u.m) - body.position[1].to(u.m))**2 + softening**2)
            acceleration[0] += (G * body.mass) / (distance**2)
            acceleration[1] += (G * body.mass) / (distance**2)
        return acceleration
    
    def move(self, time_step):
        """
        ### Parameters
        - time_step: int days
        """
        self.acceleration = self.calc_acceleration()
        time_step = time_step * u.day
        self.velocity[0] += self.acceleration[0].to(u.km / u.s**2) * time_step.to(u.s)
        self.velocity[1] += self.acceleration[1].to(u.km / u.s**2) * time_step.to(u.s)

        self.position[0] += self.velocity[0].to(u.au / u.s) * time_step.to(u.s)
        self.position[1] += self.velocity[1].to(u.au / u.s) * time_step.to(u.s)

        if time_step > 0:
            self.past_positions.append([self.position[0].value, self.position[1].value])
        else:
            self.past_positions.pop()