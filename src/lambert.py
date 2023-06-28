import matplotlib.pyplot as plt
import numpy as np
from astropy.constants import G, M_sun, M_earth
import astropy.units as u
from time import sleep

class Point():
    def __init__(self, deg, semimajor, semiminor, rotX):
        """
        ### Parameters
        - deg: float, in degrees
        - semimajor: float, in AU
        - semiminor: float, in AU
        - rotX: float, in degrees
        """
        self.deg = deg * u.deg
        self.rad = self.deg.to(u.rad)

        self.rotX = rotX
        self.rotX_rad = ((90 - (90 + self.rotX)) * u.deg).to(u.rad)

        self.x = (semimajor * np.cos(self.rad) * np.cos(self.rotX_rad))
        self.y = (semiminor * np.sin(self.rad))
        self.z = u.au * (np.cos(self.rad) * np.sin(self.rotX_rad))

    def __str__(self):
        return f"Point: ({self.x:.3f}, {self.y:.3f}, {self.z:.3f})"

FIG = plt.figure()
FIG.set_size_inches(6, 6)
AX = FIG.add_subplot(111, projection="3d")
AX.set_xlabel("X")
AX.set_ylabel("Y")
AX.set_zlabel("Z")
AX.set_xlim(-5, 5)
AX.set_ylim(-5, 5)
AX.set_zlim(-0.5, 0.5)
AX.view_init(90, -90)

# https://www.youtube.com/watch?v=XltNtuw6P44
# http://control.asu.edu/Classes/MAE462/462Lecture10.pdf
def modern_lamberts(a, alpha, beta):
    """
    ### Parameters
    - a: float, in m
    - alpha: float
    - beta: float
    ### Returns
    - dt: Time of flight in seconds
    """
    mu = M_sun * G
    return np.sqrt(a**3 / mu) * (alpha - beta - (np.sin(alpha) - np.sin(beta)))

def alpha(s, a):
    """
    ### Returns
    - Angle in radians, dimensionless
    """
    return (2 * np.arcsin(np.sqrt(s / (2 * a)))).value

def beta(s, c, a):
    """
    ### Returns
    - Angle in radians, dimensionless
    """
    return (2 * np.arcsin(np.sqrt((s - c) / (2 * a)))).value

def min_tof(s, c):
    """
    ### Returns
    - Time of flight in seconds
    """
    mu = M_sun * G
    return (np.sqrt(2) / 3) * np.sqrt(s**3 / mu) * (1 - ((s - c) / s)**1.5)

def max_tof(a_min, s, c):
    """
    ### Returns
    - Time of flight in seconds
    """
    mu = M_sun * G
    return np.sqrt(a_min**3 / mu) * (alpha(s, a_min) - beta(s, c, a_min) - (np.sin(alpha(s, a_min)) - np.sin(beta(s, c, a_min))))

def draw_circle(r, rotX, cX, cY):
    """
    ### Parameters
    - r: float, in AU
    - rotX: float, in degrees
    - cX: float, in AU
    - cY: float, in AU
    - cZ: float, in AU

    ### Returns
    - r: float, in AU
    - cX: float, in AU
    - cY: float, in AU
    """
    theta = np.linspace(0, 2 * np.pi, 1000) * u.rad

    x = (r * np.cos(theta) * np.cos(rotX)) + cX
    y = (r * np.sin(theta)) + cY
    # z = (np.cos(theta) * np.sin(rotX)) * u.au + cZ

    AX.plot3D(x, y, color="black", linestyle="--", linewidth=0.5)
    return r, cX, cY

# Solving Lambert's problem - Bisection
# https://www.youtube.com/watch?v=a3BWBYXi-IM
def bisection(s, c, tof):
    """
    ### Parameters
    - s: float, semiperimeter in m
    - c: float, chord distance in m
    - tof: float, time of flight in seconds
    """
    a_min = s / 2
    a_max = 2 * s
    dt = (1000 * u.yr).to(u.s)

    while dt != tof:
        # print(round((dt * u.s).to(u.d).value, 3))
        prev = dt
        a = (a_min + a_max) / 2
        A = alpha(s, a)
        B = beta(s, c, a)

        dt = modern_lamberts(a, A, B)
        if dt > tof:
            a_min = a
        else:
            a_max = a
        if prev - dt < 0.0001 * u.s:
            break
    return a, dt

# http://ambrnet.com/TrigoCalc/Circles2/circle2intersection/CircleCircleIntersection.htm#jscriptAreaBetween2circles
def calc_dist_2_points(c1x, c1y, c2x, c2y):
    """
    ### Parameters
    - c1x: float, in AU
    - c1y: float, in AU
    - c2x: float, in AU
    - c2y: float, in AU

    ### Returns
    - distance: float, in AU
    """
    return (np.sqrt((c1x - c2x)**2 + (c1y - c2y)**2))

def calc_semiperimeter(chord_distance, r1, r2):
    """
    ### Parameters
    - chord_distance: float, in AU
    - r1: float, in AU
    - r2: float, in AU

    ### Returns
    - s: float, in AU
    """

    perimeter = (chord_distance + r1 + r2) # AU
    s = perimeter / 2 # semiperimeter, AU
    return s

def calc_area_herons(r1, c1x, c1y, r2, c2x, c2y):
    """
    ### Parameters
    - r1: float, in AU
    - c1x: float, in AU
    - c1y: float, in AU
    - r2: float, in AU
    - c2x: float, in AU
    - c2y: float, in AU

    ### Returns
    - area: float, in AU
    """
    chord_distance = calc_dist_2_points(c1x, c1y, c2x, c2y) # AU
    s = calc_semiperimeter(chord_distance, r1, r2) # AU

    return np.sqrt(s * (s - r1) * (s - r2) * (s - chord_distance))

def calc_x_intersections(r1, c1x, c1y, r2, c2x, c2y):
    """
    ### Parameters
    - r1: float, in AU
    - c1x: float, in AU
    - c1y: float, in AU
    - r2: float, in AU
    - c2x: float, in AU
    - c2y: float, in AU

    ### Returns
    - x1: float, in AU
    - x2: float, in AU
    """
    a = c1x
    b = c1y
    c = c2x
    d = c2y

    chord_length = calc_dist_2_points(c1x, c1y, c2x, c2y)
    area = calc_area_herons(r1, c1x, c1y, r2, c2x, c2y)

    term1 = ((a + c) / 2) + (((c - a) * (r1**2 - r2**2)) / (2 * chord_length**2))
    term2 = 2 * ((b - d) / chord_length**2) * area
    return term1 + term2, term1 - term2

def calc_y_intersections(r1, c1x, c1y, r2, c2x, c2y):
    """
    ### Parameters
    - r1: float, in AU
    - c1x: float, in AU
    - c1y: float, in AU
    - r2: float, in AU
    - c2x: float, in AU
    - c2y: float, in AU
    
    ### Returns
    - y1: float, in AU
    - y2: float, in AU
    """
    a = c1x
    b = c1y
    c = c2x
    d = c2y

    chord_length = calc_dist_2_points(c1x, c1y, c2x, c2y)
    area = calc_area_herons(r1, c1x, c1y, r2, c2x, c2y)

    term1 = ((b + d) / 2) + (((d - b) * (r1**2 - r2**2)) / (2 * chord_length**2))
    term2 = 2 * ((a - c) / chord_length**2) * area
    return term1 - term2, term1 + term2

def calc_semimajor_axis(planetX, planetY, f1x, f1y, f2x, f2y):
    """
    ### Parameters
    - planetX: float, in AU
    - planetY: float, in AU
    - f1x: float, in AU
    - f1y: float, in AU
    - f2x: float, in AU
    - f2y: float, in AU

    ### Returns
    - a: float, in AU
    """
    return (calc_dist_2_points(planetX, planetY, f1x, f1y) + calc_dist_2_points(planetX, planetY, f2x, f2y)) / 2

def calc_midpoint(c1x, c1y, c2x, c2y):
    """
    ### Parameters
    - c1x: float, in AU
    - c1y: float, in AU
    - c2x: float, in AU
    - c2y: float, in AU

    ### Returns
    - midpointX: float, in AU
    - midpointY: float, in AU
    """
    return (c1x + c2x) / 2, (c1y + c2y) / 2

def calc_semiminor_axis(semimajor, f1x, f1y, f2x, f2y):
    """
    ### Parameters
    - semimajor: float, in AU
    - f1x: float, in AU
    - f1y: float, in AU
    - f2x: float, in AU
    - f2y: float, in AU

    ### Returns
    - b: float, in AU
    """
    midX, midY = calc_midpoint(f1x, f1y, f2x, f2y)
    c = calc_dist_2_points(f1x, f1y, midX, midY)
    b = np.sqrt(semimajor**2 - c**2)
    return b

def calc_slope_angle_foci(f1x, f1y, f2x, f2y):
    """
    ### Parameters
    - f1x: float, in AU
    - f1y: float, in AU
    - f2x: float, in AU
    - f2y: float, in AU

    ### Returns
    - angle: float, in radians
    """
    slope = (f2y - f1y) / (f2x - f1x)
    return np.arctan(slope)

def plot_transfer(planet1, planet2, semimajor_axis):
    sunX = 0
    sunY = 0

    planet1_position = planet1.rad
    planet2_position = planet2.rad

    planet1_to_sun = calc_dist_2_points(planet1.x, planet1.y, sunX, sunY)
    planet2_to_sun = calc_dist_2_points(planet2.x, planet2.y, sunX, sunY)
    r1 = (2*semimajor_axis).to(u.au) - planet1_to_sun
    r2 = (2*semimajor_axis).to(u.au) - planet2_to_sun

    x1, x2 = calc_x_intersections(r1, planet1.x, planet1.y, r2, planet2.x, planet2.y)
    y1, y2 = calc_y_intersections(r1, planet1.x, planet1.y, r2, planet2.x, planet2.y)

    semimajor1 = calc_semimajor_axis(planet1.x, planet1.y, sunX, sunY, x1, y1)
    semimajor2 = calc_semimajor_axis(planet1.x, planet1.y, sunX, sunY, x2, y2)

    semiminor1 = calc_semiminor_axis(semimajor1, sunX, sunY, x1, y1)
    semiminor2 = calc_semiminor_axis(semimajor2, sunX, sunY, x2, y2)

    angle1 = calc_slope_angle_foci(sunX, sunY, x1, y1)
    angle2 = calc_slope_angle_foci(sunX, sunY, x2, y2)

    center1x, center1y = calc_midpoint(sunX, sunY, x1, y1)
    center2x, center2y = calc_midpoint(sunX, sunY, x2, y2)

    theta = np.linspace(0, 2 * np.pi, 1000) * u.rad

    def x(centerx, semimajor, semiminor, theta, angle):
        return centerx + ((semimajor * np.cos(theta) * np.cos(angle)) - (semiminor * np.sin(theta) * np.sin(angle)))

    def y(centery, semimajor, semiminor, theta, angle):
        return centery + ((semimajor * np.cos(theta) * np.sin(angle)) + (semiminor * np.sin(theta) * np.cos(angle)))

    x1 = x(center1x, semimajor1, semiminor1, theta, angle1)
    y1 = y(center1y, semimajor1, semiminor1, theta, angle1)

    AX.plot(x1, y1, color="green")

    x2 = x(center2x, semimajor2, semiminor2, theta, angle2)
    y2 = y(center2y, semimajor2, semiminor2, theta, angle2)
    
    AX.plot(x2, y2, color="green")

    return x1, x2, y1, y2

theta = np.linspace(0, 2 * np.pi, 1000) * u.rad

sunX, sunY, sunZ = (0 * u.au), (0 * u.au), (0 * u.au)
AX.scatter(sunX.value, sunY.value, sunZ.value, color="orange")
AX.text(sunX.value, sunY.value, sunZ.value + 0.005, "Sun", size=8)

point1 = Point(0, (1 * u.au), (0.98 * u.au), 0)
print(point1)
AX.scatter(point1.x.value, point1.y.value, point1.z.value)
r1 = (np.sqrt((point1.x - sunX)**2 + (point1.y - sunY)**2)) # AU

point2 = Point(150, (1.52 * u.au), (1.51 * u.au), 7)
print(point2)
AX.scatter(point2.x.value, point2.y.value, point2.z.value)
r2 = (np.sqrt((point2.x - sunX)**2 + (point2.y - sunY)**2)) # AU

chord_distance = calc_dist_2_points(point1.x, point1.y, point2.x, point2.y).to(u.m)
s = calc_semiperimeter(chord_distance, r1, r2).to(u.m)

print(f"Min TOF: {(min_tof(s, chord_distance)).to(u.d):.3f}")
print(f"Max TOF: {(max_tof(s / 2, s, chord_distance)).to(u.d):.3f}")

delta_t = (200 * u.d).to(u.s)
if delta_t > min_tof(s, chord_distance):
    a, t = bisection(s, chord_distance, delta_t)

    r1 = (2*a).to(u.au) - r1
    r2 = (2*a).to(u.au) - r2
    
    draw_circle(r1, 0, point1.x, point1.y)
    draw_circle(r2, 0, point2.x, point2.y)

    x1, x2 = calc_x_intersections(r1, point1.x, point1.y, r2, point2.x, point2.y)
    y1, y2 = calc_y_intersections(r1, point1.x, point1.y, r2, point2.x, point2.y)

    AX.scatter(x1.value, y1.value, 0, color="red")
    AX.scatter(x2.value, y2.value, 0, color="red")

    plot_transfer(point1, point2, a)
    

    plt.show()
    
