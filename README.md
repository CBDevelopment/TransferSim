# TransferSim

- A Python program simulating transfer orbits between planets in the solar system
- Developed by Cutter Beck in collaboration with Chris Baxter and Sean Nuzio at Worcester Polytechnic Institute (WPI)

![gifs\transfer_no_circles.gif](https://github.com/CBDevelopment/TransferSim/blob/main/gifs/transfer_no_circles.gif)

## Getting Started
- Install Python 3.11 or higher
- Install all dependencies found in `requirements.txt` using the command `pip install -r requirements.txt`


## Running Gravity Assists

- Open `src/gravity-assists` in a terminal
- Run `python gravity-assists.py`
    - In the current configuration, clicking +10 Days will increment the simulation by 10 days, and -10 Days will decrement the simulation by 10 days
- The following variables can be changed to set how the simulation runs
    - `YEAR` will set the starting year of the simulation
    - `MONTH` will set the starting month of the simulation
    - `DAY` will set the starting day of the simulation
    - `BUTTONS` determines whether the simulation will run with interactive buttons or if it will run automatically in the background and save a GIF
    - `STEP` determines how many days at a time the simulation will increase, this changes the value on the buttons as well
    - `SUN_SCALE` determines the relative size of the Sun in the simulation
    - `PLANET_SCALE` determines the relative size of the planets in the simulation
    - `CRAFT_SCALE` determines the relative size of the spacecraft in the simulation
    - `TO_DRAW` determines which planets will be drawn in the simulation
    - Changing `flight_angle` and `delta_v` will change the flight path of the spacecraft
    - Changing `craft` will change starting values of the spacecraft

## Solving Lambert's Problem

- Open `src/lamberts` in a terminal
- Run `python lambert.py`
    - This calculates the 4 possible paths between two bodies, like Earth and Mars, based on a solved Lambert's problem
    - It prints the following in the current configuration:
        ```
        Point: (1.000 AU, 0.000 AU, 0.000 AU)
        Point: (-1.307 AU, 0.755 AU, 0.106 AU)
        Min TOF: 106.022 d
        Max TOF: 250.115 d
        ```
        - The points are the X, Y, and Z coordinates of each of the bodies (i.e. Earth and Mars)
        - The minimum and maximum time of flight (TOF) are the minimum and maximum time it would take to travel between the two bodies
- To edit the configuration, edit lines 348 and 353, which are Point objects, which represent the locations of the bodies
    ```
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
    ```
    ```
    348 point1 = Point(0, (1 * u.au), (0.98 * u.au), 0)
    353 point2 = Point(150, (1.52 * u.au), (1.51 * u.au), 7)
    ```

## Viewing a Transfer Orbit Solved by Lambert's Problem

- Open `src/lamberts` in a terminal
- Run `python transfer_orbit.py`
    - This will show a 3D plot of the transfer orbit between two bodies, like Earth and Mars, based on a solved Lambert's problem
    - It also prints out the maximum and minimum time of flight (TOF) between the two bodies
- To edit the configuration,
    - Change the string `START_YEAR` to the starting year of the simulation
    - Change the string `START_MONTH` to the starting month of the simulation
    - Change the string `START_DAY` to the starting day of the simulation
    - Change `planets_to_sim` to be a list of names of the planets you want to simulate
        - The list of available planets to simulate can be seen in the dictionary `PLANETS`
    - Change `STEP_DAYS` to be the number of days you want to increment the simulation by
    - Change `SUN_SCALE` to be the relative size of the Sun in the simulation
    - Change `PLANET_SCALE` to be the relative size of the planets in the simulation
- To run the simulation manually, set `ANIMATE` to `False`, to run it in the background and save a GIF, set `ANIMATE` to `True`
    - `anim_length` determines the number of frames in the GIF; the number of days simulated can be inferred by `anim_length` * `STEP_DAYS`
