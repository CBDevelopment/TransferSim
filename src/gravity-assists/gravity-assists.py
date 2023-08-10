
import numpy as np
from astropy import units as u
import matplotlib.pyplot as plt
from skyfield.api import load
from skyfield import elementslib
from body import Body
from spacecraft import Spacecraft
from matplotlib.widgets import Button
from datetime import datetime, timedelta
from astropy.constants import G, M_sun
from matplotlib import animation
from tqdm import tqdm
import sys

YEAR = 2023
MONTH = 8
DAY = 8
date = datetime(YEAR, MONTH, DAY)
BUTTONS = True

# Load the planetary ephemeris
eph = load('de440.bsp')

# Load objects for all of the planets
sun = eph['sun']
mercury = eph['mercury barycenter']
venus = eph['venus barycenter']
earth = eph['earth barycenter']
mars = eph['mars barycenter']
jupiter = eph['jupiter barycenter']
saturn = eph['saturn barycenter']
uranus = eph['uranus barycenter']
neptune = eph['neptune barycenter']
pluto = eph['pluto barycenter']

# Generate a range of times
ts = load.timescale()
t = ts.utc(date.year, date.month, date.day)

# Dictionary of all the planets
BODIES = {
    "Sun": Body("Sun", t, sun, 695700 * u.km, 1.989*10**30),
    "Mercury": Body("Mercury", t, mercury, 4879/2 * u.km, 0.330*10**24),
    "Venus": Body("Venus", t, venus, 12104/2 * u.km, 4.87*10**24),
    "Earth": Body("Earth", t, earth, 12756/2 * u.km, 5.97*10**24),
    "Mars": Body("Mars", t, mars, 6792/2 * u.km, 0.642*10**24),
    "Jupiter": Body("Jupiter", t, jupiter, 142984/2 * u.km, 1898*10**24),
    "Saturn": Body("Saturn", t, saturn, 120536/2 * u.km, 568*10**24),
    "Uranus": Body("Uranus", t, uranus, 51118/2 * u.km, 86.8*10**24),
    "Neptune": Body("Neptune", t, neptune, 49528/2 * u.km, 102*10**24)
}

# Running simulator
FIG = plt.figure()
FIG.suptitle("Solar System", fontsize=16)
FIG.set_size_inches(8, 8)
AX = FIG.add_subplot(111)

def plot_body_position(body_name, color=None, scale=1):
    body: Body = BODIES[body_name]
    body_circ = plt.Circle((body.position[0].value, body.position[1].value), scale * body.radius.to(u.au).value, color=color, label=body_name)
    AX.add_patch(body_circ)
    if body_name != "Sun":
        AX.text(body.position[0].value, body.position[1].value, body_name, color="black", fontsize=8)

def plot_soi(body_name):
    body: Body = BODIES[body_name]
    soi = plt.Circle((body.position[0].value, body.position[1].value), body.r_soi.value, color="black", fill=False)
    AX.add_patch(soi)

def plot_orbit(body_name, length_days=None):
    global date
    body: Body = BODIES[body_name]

    orbit_ts = load.timescale()
    if length_days is None:
        orbit_t = orbit_ts.utc(date.year, date.month, range(date.day, date.day + int(body.orbital_elements.period_in_days) + 5))
    else:
        orbit_t = orbit_ts.utc(date.year, date.month, range(date.day, date.day + length_days))

    planet_pos = body.get_body_position(orbit_t).position.au
    AX.plot(planet_pos[0], planet_pos[1], label=body_name + " Orbit", color="black", alpha=0.3)

def plot_spacecraft_position(spacecraft: Spacecraft, scale=1):
    craft_circ = plt.Circle((spacecraft.position[0].value, spacecraft.position[1].value), scale * spacecraft.radius.to(u.au).value, color="red", label=spacecraft.name)
    AX.add_patch(craft_circ)
    AX.plot([pos[0] for pos in spacecraft.past_positions], [pos[1] for pos in spacecraft.past_positions], alpha=0.5, label=spacecraft.name + " Path")
    AX.text(spacecraft.position[0].value, spacecraft.position[1].value, spacecraft.name, color="black", fontsize=8)

DAYS = 200
STEP = 10
SUN_SCALE = 10
PLANET_SCALE = 1000
CRAFT_SCALE = 10000

ALL_PLANETS = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]
INNER_PLANETS = ["Mercury", "Venus", "Earth", "Mars"]
JOVIAN_PLANETS = ["Jupiter", "Saturn", "Uranus", "Neptune"]
EARTH_MARS = ["Earth", "Mars"]
TO_DRAW = INNER_PLANETS

earth = BODIES["Earth"]
e = earth.orbital_elements.eccentricity
ta = earth.orbital_elements.true_anomaly.to(u.rad)
# flight_angle = np.arctan((e * np.sin(ta)) / (1 + (e * np.cos(ta))))
flight_angle = 2 * u.rad
print(flight_angle)
delta_v = 8 * (u.km / u.s)

craft = Spacecraft(
            "Ingenuity", 
            1000 * u.kg,
            10 * u.m, 
            [earth.position[0] + (0.03 * u.au), earth.position[1] + (0.03 * u.au)],
            [delta_v * np.cos(flight_angle), delta_v * np.sin(flight_angle)],
            [BODIES[planet] for planet in TO_DRAW],
            BODIES["Sun"]
            )
print(craft.velocity)

def draw(bodies, time=t):
    global date, velocities
    AX.set_xlabel("X (AU)")
    AX.set_ylabel("Y (AU)")
    AX.grid()
    plot_body_position("Sun", color="orange", scale=SUN_SCALE)
    for body in bodies:
        body_obj: Body = BODIES[body]
        body_obj.set_time(time)
        plot_orbit(body)
        plot_body_position(body, scale=PLANET_SCALE)
        plot_soi(body)
    plot_spacecraft_position(craft, scale=CRAFT_SCALE)
    AX.text(0, 0, f"{date.strftime('%b %d, %Y')}", color="black", fontsize=8)

draw(TO_DRAW)

if BUTTONS:
    def increase(i):
        global date
        AX.clear()
        date = date + timedelta(days=STEP)
        craft.move(STEP)
        draw(TO_DRAW, time=ts.utc(date.year, date.month, date.day))
        plt.draw()

    def decrease(i):
        global date
        AX.clear()
        date = date - timedelta(days=STEP)
        craft.move(-STEP)
        draw(TO_DRAW, time=ts.utc(date.year, date.month, date.day))
        plt.draw()

    axes_increase = plt.axes([0.5, 0, 0.2, 0.05])
    bincrease = Button(axes_increase, label=f"+{STEP} Days", color="yellow")
    bincrease.on_clicked(increase)
    axes_decrease = plt.axes([0.7, 0, 0.2, 0.05])
    bdecrease = Button(axes_decrease, label=f"-{STEP} Days", color="yellow")
    bdecrease.on_clicked(decrease)

    plt.show()
else:
    def increase(i):
        global date
        AX.clear()
        date = date + timedelta(days=STEP)
        craft.move(STEP)
        draw(TO_DRAW, time=ts.utc(date.year, date.month, date.day))
        plt.draw()
    
    ani = animation.FuncAnimation(FIG, increase, frames=tqdm(range(10), file=sys.stdout))
    writergif = animation.PillowWriter(fps=60)
    ani.save("D:\WPI\Junior Year\Summer Classes\Solar Systems\TransferSim\src\gravity-assists\gifs\parabola.gif", writer=writergif)