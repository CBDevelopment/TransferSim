import math
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from astropy.constants import G, M_sun
import astropy.units as u
from planet import Planet
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Button
from planet import Planet
import datetime
import sys
from tqdm import tqdm
from satellite import Satellite

FIG = plt.figure()
FIG.suptitle("Solar System", fontsize=16)
FIG.set_size_inches(8, 8)
FIG.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
AX = FIG.add_subplot(111, projection="3d")

class Simulator():
    """
    By default, creates a simulation controllable by buttons.
    """
    def __init__(self, planet_dict: dict[str: Planet], planets_to_sim: list[str], planet_scale: int, sun: bool, sun_scale: int, step_days: int, year: str, month: str, day: str, satellite: Satellite=None, dim3: bool=False, stride: int=50, buttons: bool=True):
        self.planet_dict = planet_dict
        self.planets_to_sim = planets_to_sim
        self.planet_scale = planet_scale
        self.sun = sun
        self.sun_scale = sun_scale
        self.step_days = step_days
        self.satellite = satellite
        self.dim3 = dim3
        self.stride = stride
        self.buttons = buttons
        self.total_days = 0
        self.start_date = datetime.datetime(int(year), int(month), int(day))

    def right_ascension_to_deg(right_ascension: str):
        """
        ### Parameters
        - right_ascension: str, in the form of "hh mm ss.ss"

        ### Returns
        - deg: float, in degrees
        """
        hours, minutes, seconds = right_ascension.split(" ")
        return ((float(hours) + float(minutes) / 60 + float(seconds) / 3600) * 15) * u.deg

    def plot_sun(self):
        # https://www.space.com/17001-how-big-is-the-sun-size-of-the-sun.html
        r = (696000 * u.km).to(u.AU).value * self.sun_scale
        # print(f"Sun Radius: {r:.4f}")
        theta = np.linspace(0, 2 * np.pi, 1000)

        if self.dim3:
            v = np.linspace(0, np.pi, 1000)

            x = r * np.outer(np.cos(theta), np.sin(v))
            y = r * np.outer(np.sin(theta), np.sin(v))
            z = r * np.outer(np.ones(np.size(theta)), np.cos(v))
            AX.plot_surface(x, y, z, rstride=self.stride, cstride=self.stride, color="orange")
        else:
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            AX.plot(x, y, 0, color="orange")

        
    def plot_orbit(self, planet: Planet):
        theta = np.linspace(0, 2 * np.pi, 1000) * u.rad

        r = planet.calc_distance_from_sun(theta)
        rotX = planet.orbital_inclination.to(u.rad).value

        x = r * np.cos(theta) * np.cos(rotX)
        y = r * np.sin(theta)
        z = r * np.cos(theta) * np.sin(rotX)

        AX.plot3D(x, y, z, linewidth=0.75, alpha=0.5, label=planet)

        periX = planet.perihelion.value * math.cos(0) * math.cos(rotX)
        periY = planet.perihelion.value * math.sin(0)
        periZ = planet.perihelion.value * math.cos(0) * math.sin(rotX)
        AX.text(periX, periY, periZ, planet.name[0], color="black", fontsize=8)
    
    def plot_multiple_orbits(self, planet_list: list[str]):
        for planet in planet_list:
            self.plot_orbit(self.planet_dict[planet])

    def plot_planet(self, planet: Planet, angular_pos=0):
        color = "black"
        if planet.name == "Earth":
            color = "blue"
        elif planet.name == "Mars":
            color = "red"

        angular_pos = angular_pos.to(u.rad)
            
        theta = np.linspace(0, 2 * np.pi, 1000)
        radius = planet.radius.to(u.AU).value * self.planet_scale

        rotX = planet.orbital_inclination.to(u.rad).value

        r = planet.calc_distance_from_sun(angular_pos)
        incline_X = r * math.cos(angular_pos.value) * math.cos(rotX)
        incline_Y = r * math.sin(angular_pos.value)
        incline_Z = r * math.cos(angular_pos.value) * math.sin(rotX)

        if self.dim3:
            v = np.linspace(0, np.pi, 1000)
            x = (radius * np.outer(np.cos(theta), np.sin(v))) + incline_X.value
            y = (radius * np.outer(np.sin(theta), np.sin(v))) + incline_Y.value
            z = (radius * np.outer(np.ones(np.size(theta)), np.cos(v))) + incline_Z.value
            AX.plot_surface(x, y, z, rstride=self.stride, cstride=self.stride, color=color)
        else:
            x = (radius * np.cos(theta)) + incline_X.value
            y = (radius * np.sin(theta)) + incline_Y.value
            AX.plot(x, y, 0, color=color)

    def draw_model(self, time_interval_days: int):
        # Set up the graph axis for display
        AX.set_xlabel("X (AU)")
        AX.set_ylabel("Y (AU)")
        AX.set_zlabel("Z (AU)")
        if self.planets_to_sim == ["Earth", "Mars"]:
            # For Earth and Mars
            AX.view_init(15, -90)
            AX.set_zlim(-1, 1)
        else:
            # General
            max_semimajor = max([self.planet_dict[planet].semimajor for planet in self.planets_to_sim])
            AX.set_zlim(-max_semimajor.value, max_semimajor.value)

        # Plot celestial objects
        self.plot_multiple_orbits(self.planets_to_sim)

        for planet in self.planets_to_sim:
            planet_obj = self.planet_dict[planet]
            
            for j in range(abs(time_interval_days)):
                # planet velocity is in m/s
                # planet right ascension is in degrees
                planet_obj.velocity = planet_obj.calc_velocity(planet_obj.right_ascension.to(u.rad).value)
                angular_velocity = planet_obj.velocity / planet_obj.calc_distance_from_sun(planet_obj.right_ascension.to(u.rad)).to(u.m) * u.rad
                # angular velocity is in rad/s
                seconds_in_day = 86400 * u.s
                if time_interval_days > 0:
                    planet_obj.right_ascension += (angular_velocity * seconds_in_day).to(u.deg)
                else:
                    planet_obj.right_ascension -= (angular_velocity * seconds_in_day).to(u.deg)
                if planet_obj.right_ascension >= 360 * u.deg:
                    planet_obj.right_ascension -= 360 * u.deg

            # print(f"{planet_obj.name} right ascension: {planet_obj.right_ascension}")
            self.plot_planet(planet_obj, angular_pos=planet_obj.right_ascension)
        if self.sun:
            self.plot_sun()

        if self.satellite is not None:
            pass
        
        # Display the date of the planetary positions above the sun
        self.total_days += time_interval_days
        if self.total_days >= 0:
            date = self.start_date + datetime.timedelta(self.total_days)
            AX.text(0, 0, 0.08, f"{date.strftime('%b %d, %Y')}", color="black", fontsize=8)
        else:
            date = self.start_date - datetime.timedelta(abs(self.total_days))
            AX.text(0, 0, 0.08, f"{date.strftime('%b %d, %Y')}", color="black", fontsize=8)

    def run(self, anim_length: int=36, anim_save_path: str="../gifs/test.gif", days_from_start: int=0):
        # if self.buttons is true, runs in manual mode, otherwise it will create a GIF output
        if self.buttons:
            self.draw_model(time_interval_days=days_from_start)

            def increase(i):
                AX.clear()
                self.draw_model(time_interval_days=self.step_days)
                plt.draw()

            def decrease(i):
                AX.clear()
                self.draw_model(time_interval_days=-self.step_days)
                plt.draw()

            axes_increase = plt.axes([0.71, 0.000001, 0.2, 0.075])
            bincrease = Button(axes_increase, label=f"+{self.step_days} Days", color="yellow")
            bincrease.on_clicked(increase)
            axes_decrease = plt.axes([0.51, 0.000001, 0.2, 0.075])
            bdecrease = Button(axes_decrease, label=f"-{self.step_days} Days", color="yellow")
            bdecrease.on_clicked(decrease)

            plt.show()

        else:
            def increase(i):
                AX.clear()
                self.draw_model(time_interval_days=self.step_days)
            
            ani = animation.FuncAnimation(FIG, increase, frames=tqdm(range(anim_length), file=sys.stdout))
            writergif = animation.PillowWriter(fps=30)
            ani.save(anim_save_path, writer=writergif)