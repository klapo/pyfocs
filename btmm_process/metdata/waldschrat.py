import numpy as np


# Program for calculation of the wind speed & direction from x and y input data
# data might be single values or vectors

# written by Christoph Thomas, Dept. of Micrometeorology, Univ. of Bayreuth
# July 2003, lastly modified 02.24.2006

# Mode 1: Gill Solent R2
# Convention of the Gill Solent R2: x is positive for flow N --> S
# Convention of the Gill Solent R2: y is positive for flow E --> W

# Mode 2: Metek USA
# Convention of the Metek USAT: x is positive for flow S --> N, opposite to Solent R2
# Convention of the Metek USAT: y is positive for flow W --> E, opposite to Solent R2

# Mode 3: Gill Solent R3-50
# Convention of the Gill Solent R3-50: x is positive for flow S --> N
# Convention of the Gill Solent R3-50: y is positive for flow E --> W

# Mode 4: Young 81000
# Convention of the Young 81000: x is positive for flow E --> W
# Convention of the Young 81000: y is positive for flow N --> S

# Mode 5: Campbell CSAT3
# Convention of the CSAT3: x is positive for flow N --> S
# Convention of the CSAT3: y is positive for flow W --> E


# Azimuth is the direction in degrees to which the north arrow of the sonic shows

def winddirection(x, y, mode, azimuth=0, correction=1):

    if ((azimuth >= 0) & (azimuth <= 360)) & ((correction == 0) | (correction == 1)):

        # Applying sonic dependent rotations and offsets
        if mode == 1:
            x = x
            y = y
            offset = 30
        elif mode == 2:
            x = -x
            y = -y
            offset = 0
        elif mode == 3:
            x = -x
            y = y
            offset = 30
        elif mode == 4:
            x = -x
            y = y
            offset = 90
        elif mode == 5:
            x = x
            y = -y
            offset = 0

        # Calculating the horizontal wind speed
        speed = (np.round(np.sqrt(x**2 + y**2) * 100)) / 100

        # Preallocating wind direction
        direction = np.zeros_like(x)

        # Calculating the wind direction
        for i, xi in enumerate(x):
            if (y[i] < 0):  # Quadrant 1 and 2
                # +180 for conversion from vector direction to meteorological
                # wind direction convention
                direction[i] = 90 - (np.arctan(x[i] / y[i])) * (360 / (2*np.pi)) + 180 - offset
            elif (y[i] > 0):  # Quadrant 3 and 4
                direction[i] = 270 - (np.arctan(x[i] / y[i])) * (360 / (2*np.pi)) + 180 - offset
            elif (y[i] == 0):
                if (x[i] > 0):  # case of wind from north
                    direction[i] = 360 - offset
                elif (x[i] < 0):  # case of wind from south
                    direction[i] = 180 - offset
                elif (x[i] == 0) & (not i == 0):  # assign value of preceeding data
                    direction[i] = direction[i-1]
            else:
                direction[i] = np.nan

            # Correcting for azimut angle
            direction[i] = direction[i] + azimuth

            if (correction == 1):

                # Not applying this correction leads to time series with
                # hysteresis but outputs non meteorological values!
                if (direction[i] < 0):
                    direction[i] = direction[i] + 360
                elif (direction[i] > 720):
                    direction[i] = direction[i] - 720
                elif (direction[i] > 360):
                    direction[i] = direction[i] - 360
        direction = (np.round(direction * 100) / 100)

        return (speed, direction)
