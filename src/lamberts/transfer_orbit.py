from sunpy.coordinates import get_body_heliographic_stonyhurst
from astropy.time import Time
import numpy as np
from astropy import units as u
from planet import Planet
from simulator import Simulator
from satellite import Satellite

START_YEAR = "2020"
START_MONTH = "07"
START_DAY = "30"
obstime = Time(f'{START_YEAR}-{START_MONTH}-{START_DAY}T00:00:00.000')

# https://en.wikipedia.org/wiki/Semi-major_and_semi-minor_axes
# https://nssdc.gsfc.nasa.gov/planetary/factsheet/
PLANETS = {
    "Mercury": Planet("Mercury", 0.38700, 0.37870, 7.01, 4879/2, get_body_heliographic_stonyhurst("Mercury", time=obstime).lon, 0.330*10**24),
    "Venus": Planet("Venus", 0.72300, 0.72298, 3.39, 12104/2, get_body_heliographic_stonyhurst("Venus", time=obstime).lon, 4.87*10**24),
    "Earth": Planet("Earth", 1.00000, 0.99986, 0, 12756/2, get_body_heliographic_stonyhurst("Earth", time=obstime).lon, 5.97*10**24),
    "Mars": Planet("Mars", 1.52400, 1.51740, 1.85, 6792/2, get_body_heliographic_stonyhurst("Mars", time=obstime).lon, 0.642*10**24),
    "Jupiter": Planet("Jupiter", 5.20440, 5.19820, 1.31, 142984/2, get_body_heliographic_stonyhurst("Jupiter", time=obstime).lon, 1898*10**24),
    "Saturn": Planet("Saturn", 9.58260, 9.56730, 2.49, 120536/2, get_body_heliographic_stonyhurst("Saturn", time=obstime).lon, 568*10**24),
    "Uranus": Planet("Uranus", 19.21840, 19.19770, 0.77, 51118/2, get_body_heliographic_stonyhurst("Uranus", time=obstime).lon, 86.8*10**24),
    "Neptune": Planet("Neptune", 30.11000, 30.10870, 1.77, 49528/2, get_body_heliographic_stonyhurst("Neptune", time=obstime).lon, 102*10**24)
}

ALL_PLANETS = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]
EARTH_MARS = ["Earth", "Mars"]
planets_to_sim = EARTH_MARS
STEP_DAYS = 10
SUN_SCALE = 25
PLANET_SCALE = 2000
ANIMATE = False

satellite = Satellite("Loyalty", PLANETS["Earth"], PLANETS["Mars"], 1, 37000, graph_ax=None)
# satellite = None
sim = Simulator(planet_dict=PLANETS, planets_to_sim=planets_to_sim, planet_scale=PLANET_SCALE, sun=True, sun_scale=SUN_SCALE, step_days=STEP_DAYS, year=START_YEAR, month=START_MONTH, day=START_DAY, satellite=satellite, dim3=True, stride=50, buttons=(not ANIMATE))

if ANIMATE:
    sim.run(anim_length=36, anim_save_path="./gifs/EXAMPLE.gif")
else:
    sim.run(days_from_start=0)

# https://solarsystem.nasa.gov/basics/chapter5-1/#:~:text=the%20Clarke%20orbit.-,Geosynchronous%20Transfer%20Orbit,Geosynchronous%20Transfer%20Orbit%20(GTO).
# sat = Satellite("Loyalty", PLANETS["Earth"], PLANETS["Mars"], 37000)

# print(sat.calc_transfer_period())