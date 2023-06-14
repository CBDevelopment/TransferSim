import math
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from astropy.constants import G, M_sun
import astropy.units as u
from planet import Planet
import mpl_toolkits.mplot3d

PLANET = "Mars"
THETA_OF_PLANET = np.pi * u.rad
FIG = plt.figure()
AX = FIG.add_subplot(111, projection="3d")
AX.set_xlabel('X')
AX.set_ylabel('Y')
AX.set_zlabel('Z')
# AX.set_zlim(-1, 1)
AX.view_init(5, -75)
INFERIOR_PLANET_LIST = ["Mercury", "Venus", "Earth"]
SUPERIOR_PLANET_LIST = ["Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]
ALL_PLANET_LIST = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]

# https://en.wikipedia.org/wiki/Semi-major_and_semi-minor_axes
# https://www.physicsforums.com/threads/calculating-plutos-orbital-time-above-below-ecliptic.612487/
# https://en.wikipedia.org/wiki/Orbital_inclination
PLANETS = {
    "Mercury": Planet(0.38700, 0.37870, 7.01),
    "Venus": Planet(0.72300, 0.72298, 3.39),
    "Earth": Planet(1.00000, 0.99986, 0),
    "Mars": Planet(1.52400, 1.51740, 1.85),
    "Jupiter": Planet(5.20440, 5.19820, 1.31),
    "Saturn": Planet(9.58260, 9.56730, 2.49),
    "Uranus": Planet(19.21840, 19.19770, 0.77),
    "Neptune": Planet(30.11000, 30.10870, 1.77),
    "Pluto": Planet(39.48200, 34.031, 17.14)
}


def calc_distance_from_sun(theta, planets=PLANETS, planet=PLANET):
    """
    ### Parameters
    - theta: float, in radians

    ### Returns
    - distance: float, in AU
    """
    SEMIMAJOR = planets[planet].semimajor # AU
    ECCENTRICITY = planets[planet].eccentricity # dimensionless
    return SEMIMAJOR * (1 - ECCENTRICITY**2) / (1 + ECCENTRICITY * np.cos(theta))

def calc_speed(theta, planets=PLANETS, planet=PLANET):
    """
    ### Parameters
    - theta: float, in radians

    ### Returns
    - speed: float, in m / s
    """
    SEMIMAJOR = planets[planet].semimajor # AU
    ECCENTRICITY = planets[planet].eccentricity # dimensionless
    return np.sqrt((G * M_sun * SEMIMAJOR.to(u.m) * (1 - ECCENTRICITY**2)) / (calc_distance_from_sun(theta, planets, planet).to(u.m)**2))

def plot_sun(ax=AX, dim3: bool = True):
    # https://www.space.com/17001-how-big-is-the-sun-size-of-the-sun.html
    r = (696000 * u.km).to(u.AU).value
    print(f"Sun Radius: {r:.4f}")
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


def plot_orbit(ax=AX, planets=PLANETS, planet=PLANET):
    print(planet)
    print(f"Perihelion: {planets[planet].perihelion:.4f}")
    print(f"Aphelion: {planets[planet].aphelion:.4f}")

    theta = np.linspace(0, 2 * np.pi, 1000) * u.rad

    r = calc_distance_from_sun(theta, planets, planet)
    rotX = planets[planet].orbital_inclination.to(u.rad).value

    x = r * np.cos(theta) * np.cos(rotX)
    y = r * np.sin(theta)
    z = r * np.cos(theta) * np.sin(rotX)

    ax.plot3D(x, y, z, linewidth=0.75, alpha=0.5, label=planet)

    periX = planets[planet].perihelion.value * math.cos(0) * math.cos(rotX)
    periY = planets[planet].perihelion.value * math.sin(0)
    periZ = planets[planet].perihelion.value * math.cos(0) * math.sin(rotX)
    print(periX, periY, periZ)
    ax.text(periX, periY, periZ, planet[0], color="black", fontsize=8)
    return ax

def plot_multiple_orbits(ax=AX, planets=PLANETS, planet_list=ALL_PLANET_LIST):
    for planet in planet_list:
        plot_orbit(ax, planets, planet)

plot_sun(dim3=False)
# plot_orbit(planet="Earth")
no_pluto = ALL_PLANET_LIST[:-1]
plot_multiple_orbits(AX, PLANETS, ALL_PLANET_LIST)
plt.show()