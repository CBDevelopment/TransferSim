import math
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from astropy.constants import G, M_sun
import astropy.units as u
from planet import Planet
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Button

FIG = plt.figure()
FIG.suptitle("Solar System", fontsize=16)
FIG.set_size_inches(5, 5)
FIG.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
AX = FIG.add_subplot(111, projection="3d")
AX.set_xlabel("X (AU)")
AX.set_ylabel("Y (AU)")
AX.set_zlabel("Z (AU)")
INFERIOR_PLANET_LIST = ["Mercury", "Venus", "Earth"]
EARTH_MARS = ["Earth", "Mars"]
SUPERIOR_PLANET_LIST = ["Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]
ALL_PLANET_LIST = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]

def right_ascension_to_deg(right_ascension: str):
    """
    ### Parameters
    - right_ascension: str, in the form of "hh mm ss.ss"

    ### Returns
    - deg: float, in degrees
    """
    hours, minutes, seconds = right_ascension.split(" ")
    return ((float(hours) + float(minutes) / 60 + float(seconds) / 3600) * 15) * u.deg


# https://en.wikipedia.org/wiki/Semi-major_and_semi-minor_axes
# https://www.physicsforums.com/threads/calculating-plutos-orbital-time-above-below-ecliptic.612487/
# https://en.wikipedia.org/wiki/Orbital_inclination
# https://nssdc.gsfc.nasa.gov/planetary/factsheet/
# Right Ascension data from https://ssd.jpl.nasa.gov/horizons/app.html#/ pulled for June 15, 2023 from Worcester, MA
PLANETS = {
    "Mercury": Planet("Mercury", 0.38700, 0.37870, 7.01, 4879/2, right_ascension_to_deg("04 16 16.43")),
    "Venus": Planet("Venus", 0.72300, 0.72298, 3.39, 12104/2, right_ascension_to_deg("08 44 55.49")),
    "Earth": Planet("Earth", 1.00000, 0.99986, 0, 12756/2, (((5.75*np.pi)/4) * u.rad).to(u.deg)),
    "Mars": Planet("Mars", 1.52400, 1.51740, 1.85, 6792/2, right_ascension_to_deg("09 09 06.16")),
    "Jupiter": Planet("Jupiter", 5.20440, 5.19820, 1.31, 142984/2, right_ascension_to_deg("02 16 09.35")),
    "Saturn": Planet("Saturn", 9.58260, 9.56730, 2.49, 120536/2, right_ascension_to_deg("22 36 47.66")),
    "Uranus": Planet("Uranus", 19.21840, 19.19770, 0.77, 51118/2, right_ascension_to_deg("03 13 06.04")),
    "Neptune": Planet("Neptune", 30.11000, 30.10870, 1.77, 49528/2, right_ascension_to_deg("23 52 00.38")),
    "Pluto": Planet("Pluto", 39.48200, 34.031, 17.14, 2376/2, right_ascension_to_deg("20 09 27.76"))
}

def plot_sun(ax=AX, dim3: bool = True, scale=1):
    # https://www.space.com/17001-how-big-is-the-sun-size-of-the-sun.html
    r = (696000 * u.km).to(u.AU).value * scale
    # print(f"Sun Radius: {r:.4f}")
    theta = np.linspace(0, 2 * np.pi, 1000)

    if dim3:
        stride = 50

        v = np.linspace(0, np.pi, 1000)

        x = r * np.outer(np.cos(theta), np.sin(v))
        y = r * np.outer(np.sin(theta), np.sin(v))
        z = r * np.outer(np.ones(np.size(theta)), np.cos(v))
        ax.plot_surface(x, y, z, rstride=stride, cstride=stride, color="orange")
    else:
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        ax.plot(x, y, 0, color="orange")


def plot_orbit(planet: Planet, ax=AX):
    # print(planet)
    # print(f"Perihelion: {planet.perihelion:.4f}")
    # print(f"Aphelion: {planet.aphelion:.4f}")

    theta = np.linspace(0, 2 * np.pi, 1000) * u.rad

    r = planet.calc_distance_from_sun(theta).value
    rotX = planet.orbital_inclination.to(u.rad).value

    x = r * np.cos(theta) * np.cos(rotX)
    y = r * np.sin(theta)
    z = r * np.cos(theta) * np.sin(rotX)

    ax.plot3D(x, y, z, linewidth=0.75, alpha=0.5, label=planet)

    periX = planet.perihelion.value * math.cos(0) * math.cos(rotX)
    periY = planet.perihelion.value * math.sin(0)
    periZ = planet.perihelion.value * math.cos(0) * math.sin(rotX)
    ax.text(periX, periY, periZ, planet.name[0], color="black", fontsize=8)
    return ax

def plot_multiple_orbits(ax=AX, planets: dict[str: Planet]=PLANETS, planet_list: list[str]=ALL_PLANET_LIST):
    for planet in planet_list:
        plot_orbit(planets[planet], ax)

def plot_planet(planet: Planet, ax=AX, dim3: bool = True, angular_pos=0, scale: int=1):
    color = "black"
    if planet.name == "Earth":
        color = "blue"
    elif planet.name == "Mars":
        color = "red"

    angular_pos = angular_pos.to(u.rad)
        
    theta = np.linspace(0, 2 * np.pi, 1000)
    radius = planet.radius.to(u.AU).value * scale

    rotX = planet.orbital_inclination.to(u.rad).value

    r = planet.calc_distance_from_sun(angular_pos)
    incline_X = r * math.cos(angular_pos.value) * math.cos(rotX)
    incline_Y = r * math.sin(angular_pos.value)
    incline_Z = r * math.cos(angular_pos.value) * math.sin(rotX)

    if dim3:
        stride = 50

        v = np.linspace(0, np.pi, 1000)
        x = (radius * np.outer(np.cos(theta), np.sin(v))) + incline_X.value
        y = (radius * np.outer(np.sin(theta), np.sin(v))) + incline_Y.value
        z = (radius * np.outer(np.ones(np.size(theta)), np.cos(v))) + incline_Z.value
        ax.plot_surface(x, y, z, rstride=stride, cstride=stride, color=color)
    else:
        x = (radius * np.cos(theta)) + incline_X.value
        y = (radius * np.sin(theta)) + incline_Y.value
        ax.plot(x, y, 0, color=color)

# ----- UTILIZATION -----

planets_to_sim = ["Earth", "Mars"]

def run_sim(planet_scale, planets_to_plot: list, planets: dict[str: Planet]=PLANETS, sun=False, sun_scale=1, time_interval_days=0):
    # General
    # max_semimajor = max([PLANETS[planet].semimajor for planet in planets_to_plot])
    # AX.set_zlim(-max_semimajor.value, max_semimajor.value)
    # For Earth and Mars
    AX.view_init(15, -90)
    AX.set_zlim(-1, 1)

    plot_multiple_orbits(AX, planets, planets_to_plot)

    for planet in planets_to_plot:
        planet_obj = planets[planet]
        
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

        print(planet_obj)
        plot_planet(planet_obj, angular_pos=planet_obj.right_ascension, scale=planet_scale)
    if sun:
        plot_sun(scale=sun_scale)


no_pluto = ALL_PLANET_LIST[:-1]
STEP_DAYS = 10
SUN_SCALE = 10
PLANET_SCALE = 1000

def increase(i):

    AX.clear()
    run_sim(planet_scale=PLANET_SCALE, planets_to_plot=planets_to_sim, planets=PLANETS, sun=True, sun_scale=SUN_SCALE, time_interval_days=STEP_DAYS)
    plt.draw()

# ani = animation.FuncAnimation(FIG, increase, frames=36)
# writergif = animation.PillowWriter(fps=60)
# ani.save("../gifs/test.gif", writer=writergif)


def decrease(i):
    AX.clear()
    run_sim(planet_scale=PLANET_SCALE, planets_to_plot=planets_to_sim, planets=PLANETS, sun=True, sun_scale=SUN_SCALE, time_interval_days= -1 * STEP_DAYS)
    plt.draw()

run_sim(planet_scale=PLANET_SCALE, planets_to_plot=planets_to_sim, planets=PLANETS, sun=True, sun_scale=SUN_SCALE, time_interval_days=0)

axes_increase = plt.axes([0.71, 0.000001, 0.2, 0.075])
bincrease = Button(axes_increase, label=f"+{STEP_DAYS} Days", color="yellow")
bincrease.on_clicked(increase)
axes_decrease = plt.axes([0.51, 0.000001, 0.2, 0.075])
bdecrease = Button(axes_decrease, label=f"-{STEP_DAYS} Days", color="yellow")
bdecrease.on_clicked(decrease)

plt.show()