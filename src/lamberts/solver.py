import numpy as np
from astropy.constants import G, M_sun
import astropy.units as u
import matplotlib.pyplot as plt

class Solver():
    def __init__(self, planet1, planet2):
        self.planet1 = planet1
        self.planet2 = planet2

    
    # Calculations
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
        return np.sqrt(a_min**3 / mu) * (Solver.alpha(s, a_min) - Solver.beta(s, c, a_min) - (np.sin(Solver.alpha(s, a_min)) - np.sin(Solver.beta(s, c, a_min))))

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
            A = Solver.alpha(s, a)
            B = Solver.beta(s, c, a)

            dt = Solver.modern_lamberts(a, A, B)
            if dt > tof:
                a_min = a
            else:
                a_max = a
            if prev - dt < 0.0001 * u.s:
                break
        return a, dt
    
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
        chord_distance = Solver.calc_dist_2_points(c1x, c1y, c2x, c2y) # AU
        s = Solver.calc_semiperimeter(chord_distance, r1, r2) # AU

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

        chord_length = Solver.calc_dist_2_points(c1x, c1y, c2x, c2y)
        area = Solver.calc_area_herons(r1, c1x, c1y, r2, c2x, c2y)

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

        chord_length = Solver.calc_dist_2_points(c1x, c1y, c2x, c2y)
        area = Solver.calc_area_herons(r1, c1x, c1y, r2, c2x, c2y)

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
        return (Solver.calc_dist_2_points(planetX, planetY, f1x, f1y) + Solver.calc_dist_2_points(planetX, planetY, f2x, f2y)) / 2

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
        midX, midY = Solver.calc_midpoint(f1x, f1y, f2x, f2y)
        c = Solver.calc_dist_2_points(f1x, f1y, midX, midY)
        b = np.sqrt(semimajor**2 - c**2)
        return b

    def calc_slope_angle(f1x, f1y, f2x, f2y):
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


    # Visualizations
    def draw_circle(r, rotX, cX, cY, AX):
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
    
    def plot_transfer(planet_x1, planet_y1, planet_x2, planet_y2, semimajor_axis, AX):
        sunX = 0
        sunY = 0

        planet1_to_sun = Solver.calc_dist_2_points(planet_x1, planet_y1, sunX, sunY)
        planet2_to_sun = Solver.calc_dist_2_points(planet_x2, planet_y2, sunX, sunY)
        r1 = (2*semimajor_axis).to(u.au) - planet1_to_sun
        r2 = (2*semimajor_axis).to(u.au) - planet2_to_sun

        x1, x2 = Solver.calc_x_intersections(r1, planet_x1, planet_y1, r2, planet_x2, planet_y2)
        y1, y2 = Solver.calc_y_intersections(r1, planet_x1, planet_y1, r2, planet_x2, planet_y2)

        semimajor1 = Solver.calc_semimajor_axis(planet_x1, planet_y1, sunX, sunY, x1, y1)
        semimajor2 = Solver.calc_semimajor_axis(planet_x1, planet_y1, sunX, sunY, x2, y2)

        semiminor1 = Solver.calc_semiminor_axis(semimajor1, sunX, sunY, x1, y1)
        semiminor2 = Solver.calc_semiminor_axis(semimajor2, sunX, sunY, x2, y2)

        angle1 = Solver.calc_slope_angle(sunX, sunY, x1, y1)
        angle2 = Solver.calc_slope_angle(sunX, sunY, x2, y2)

        center1x, center1y = Solver.calc_midpoint(sunX, sunY, x1, y1)
        center2x, center2y = Solver.calc_midpoint(sunX, sunY, x2, y2)

        theta = np.linspace(0, 2 * np.pi, 1000) * u.rad

        def x(centerx, semimajor, semiminor, theta, angle):
            return centerx + ((semimajor * np.cos(theta) * np.cos(angle)) - (semiminor * np.sin(theta) * np.sin(angle)))

        def y(centery, semimajor, semiminor, theta, angle):
            return centery + ((semimajor * np.cos(theta) * np.sin(angle)) + (semiminor * np.sin(theta) * np.cos(angle)))

        x1 = x(center1x, semimajor1, semiminor1, theta, angle1)
        y1 = y(center1y, semimajor1, semiminor1, theta, angle1)

        AX.plot(x1, y1, color="green", linestyle="--", linewidth=0.5)

        x2 = x(center2x, semimajor2, semiminor2, theta, angle2)
        y2 = y(center2y, semimajor2, semiminor2, theta, angle2)
        
        AX.plot(x2, y2, color="green", linestyle="--", linewidth=0.5)

        # # Calculate angle to planets for ellipse 1
        # x1_prime = planet_x1 - center1x
        # y1_prime = planet_y1 - center1y
        # x1_pp = (x1_prime.value * np.cos(angle1) + y1_prime.value * np.sin(angle1)) * u.au
        # y1_pp = (-x1_prime.value * np.sin(angle1) + y1_prime.value * np.cos(angle1)) * u.au
        # # print(x1_pp, y1_pp)
        
        # x2_prime = planet_x2 - center1x
        # y2_prime = planet_y2 - center1y
        # x2_pp = (x2_prime.value * np.cos(angle1) + y2_prime.value * np.sin(angle1)) * u.au
        # y2_pp = (-x2_prime.value * np.sin(angle1) + y2_prime.value * np.cos(angle1)) * u.au
        # # print(x2_pp, y2_pp)

        # start = np.arctan2(y1_pp, x1_pp)
        # end = np.arctan2(y2_pp, x2_pp)
        # print(start, end)

        # theta1 = np.linspace(start, end, 50)

        # x_opt1 = x(center1x, semimajor1, semiminor1, theta1, angle1)
        # y_opt1 = y(center1y, semimajor1, semiminor1, theta1, angle1)

        # AX.plot(x_opt1, y_opt1, color="purple")

        return x1, x2, y1, y2