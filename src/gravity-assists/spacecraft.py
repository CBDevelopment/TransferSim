from astropy import units as u
from astropy.constants import G, M_sun
import numpy as np
from body import Body

class Spacecraft():
    def __init__(self, name, mass, radius, position_vec, velocity_vec, bodies: list[Body], sun):
        self.name = name
        self.bodies = bodies
        self.mass = mass
        self.radius = radius.to(u.au)
        self.sun = sun

        self.velocity = velocity_vec # [x, y] km / s
        self.position = position_vec # [x, y] au
        self.acceleration = self.calc_acceleration() # [x, y] m/s^2
        self.past_positions = [[self.position[0].value, self.position[1].value]]
        
    
    def calc_acceleration(self):
        """
        ### Returns
        - [x_accel, y_accel]
        """
        acceleration = [0, 0]
        in_soi = False
        soi = None
        softening = 1 * u.m
        for body in self.bodies:
            # Distance from spacecraft to body
            dist = np.sqrt((self.position[0].to(u.m) - body.position[0].to(u.m))**2 + (self.position[1].to(u.m) - body.position[1].to(u.m))**2)
            # Check if spacecraft is within body's SOI
            if dist <= body.r_soi.to(u.m):
                in_soi = True
                soi = body
                break
        if in_soi:
            # Acceleration due to body
            print(soi.name)
            centered_pos = [self.position[0].to(u.m) - soi.position[0].to(u.m), self.position[1].to(u.m) - soi.position[1].to(u.m)]
            angle = np.arctan2(centered_pos[1], centered_pos[0])
            distance = np.sqrt((self.position[0].to(u.m) - soi.position[0].to(u.m))**2 + (self.position[1].to(u.m) - soi.position[1].to(u.m))**2 + softening**2)
            accel = (G * soi.mass) / (distance**2)

            acceleration = [0, 0]
            acceleration[0] -= accel * np.cos(angle)
            acceleration[1] -= accel * np.sin(angle)
        else:
            # Acceleration due to the Sun
            print("Sun")
            angle = np.arctan2(self.position[1].to(u.m), self.position[0].to(u.m))
            distance = np.sqrt((self.position[0].to(u.m) - self.sun.position[0].to(u.m))**2 + (self.position[1].to(u.m) - self.sun.position[1].to(u.m))**2 + softening**2)
            accel = (G * M_sun) / (distance**2)
            
            acceleration[0] -= accel * np.cos(angle)
            acceleration[1] -= accel * np.sin(angle)
        return acceleration
    
    def move(self, time_step):
        """
        ### Parameters
        - time_step: int days
        """
        if time_step > 0:
            for i in range(time_step):
                self.acceleration = self.calc_acceleration()
                step = 1 * u.day
                self.velocity[0] += self.acceleration[0].to(u.km / u.s**2) * step.to(u.s)
                self.velocity[1] += self.acceleration[1].to(u.km / u.s**2) * step.to(u.s)

                self.position[0] += ((self.velocity[0].to(u.au / u.s).value * step.to(u.s).value) + (0.5 * self.acceleration[0].to(u.au / u.s**2).value * step.to(u.s).value)**2) * u.au
                self.position[1] += ((self.velocity[1].to(u.au / u.s).value * step.to(u.s).value) + (0.5 * self.acceleration[1].to(u.au / u.s**2).value * step.to(u.s).value)**2) * u.au
                self.past_positions.append([self.position[0].value, self.position[1].value])
        else:
            self.position[0] = self.past_positions[time_step][0] * u.au
            self.position[1] = self.past_positions[time_step][1] * u.au
            for i in range(abs(time_step)):
                self.past_positions.pop()