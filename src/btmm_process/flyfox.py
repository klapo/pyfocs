import numpy as np
import sys
import xarray as xr


def hypsometric_height(T_v, p_1, p_2):
    # Equation 1, the hypsometric height
    # T_v is in Kelvin
    # pressures are in hPa (must be the same units, can vary)
    # p_1 is the lower point's pressure (e.g., the surface)
    # p_2 is the higher point's pressure (e.g., the tethersonde)
    R = 287.058  # [J/(kg K)]
    g = 9.81  # [m/s^2]

    # Make sure we are in Kelvin
    if np.max(T_v) < 200:
        T_v = T_v + 273.15

    h = (R * T_v / g) * np.log(p_1 / p_2)
    return h


def virtual_temp(T, r_v, r_l=0):
    # Equation 2, virtual temperature

    # Make sure we are in celsius
    if np.max(T) < 200:
        T = T + 273.15
    T_v = T * (1 + 0.61 * r_v - r_l)
    return T_v


def hypsometric_press(p_ref, z, z_ref, T_v):
    # Equation 4, hypsometric equation for pressure.
    # T_v is in Kelvin and is the average virtual temperature
    # between layers.
    # p_ref is the pressure at the reference height (e.g., the surface)

    # Make sure we are in Kelvin
    if np.max(T_v) < 200:
        T_v = T_v + 273.15

    R = 287.058  # [J/(kg K)]
    g = 9.81  # [m/s^2]

    p = p_ref * np.exp(-g * (z - z_ref) / (R * T_v))
    return p


def mixing_ratio(e, p):
    # Equation 3a, mixing ratio of water vapor
    # e is water vapor pressure in hPa
    # p is air pressure in hPa
    r_v = 0.622 * e / (p - 0.378 * e)
    return r_v


def e_s_wmo(T):
    # Equation 3b Saturated water vapor pressure

    # Make sure we are in celsius
    if np.max(T) > 200:
        T = T - 273.15

    # Constants
    c1 = 6.112
    c2 = 17.62
    c3 = 243.12

    # WMO approximation
    e_s = c1 * np.exp(c2 * T / (c3 + T))
    return e_s


def potential_temperature(T, p_0, p):
    # R = 287.058  # [J/(kg K)]
    # cp = 1004  # [J / (K kg)]
    R_cp = 0.286  # Use the approximate value of R/cp

    # Make sure we are in Kelvin
    if np.max(T) < 200:
        T = T + 273.15

    theta = T * (p_0 / p) ** (R_cp)
    return theta


def height_interp(section, z_bottom, num_layers=None, flip=False):
    '''
    Determine the linear height to interpolate along for each time step.
    '''
    # Number of points within each profile, this does not change with time.
    if num_layers:
        n_laf = num_layers
    else:
        n_laf = section.LAF.size

    # Create the z coordinates for this time step using similar
    # assumptions as above, but for just this time step.
    if not flip:
        z_temp = np.linspace(z_bottom, section['top'], num=n_laf)
    else:
        z_temp = np.linspace(section['top'], z_bottom, num=n_laf)
    # Interpolate the profile to the average height spacing found above.
    section.coords['z'] = ('z', z_temp)

    return section


def height_interp_time(section, z_bottom):
    '''
    Create a uniform height coordinate through time for a tethered balloon
    with a time varying height.
    '''
    # Find just the section to be labeled.
    LAF = section.LAF.values

    # Determine which direction the fiber is oriented
    if LAF[0] > LAF[1]:
        flip = True
    else:
        flip = False

    # --------
    # Determine the linear height to interpolate along for each time step.

    # Number of points within each profile, this does not change with time.
    n_laf = section.LAF.size

    # The average height of the profile. This is used for creating the
    # height step for the linear interpolation.
    z_top_avg = section['top'].mean(dim='time').values

    # The maximum height of the layer. This sets the top value that the
    # profile is interpolated to.
    z_top_max = section['top'].max(dim='time').values

    # The average height step for the profile.
    delta_z_avg = (z_top_avg - z_bottom) / n_laf

    # The interpolated height, using the average height step and the
    # maximum height of the flight.
    z_int = np.arange(z_bottom, z_top_max, delta_z_avg)

    # For each time step interpolate the the observed profile to the
    # heights derived above. The assumptions used mean that the
    # interpolated points are closer together than the delta_LAF
    # points. This is effectively downsampling.
    ds_to_label = []
    for n, t in enumerate(section.time):
        sys.stdout.write('\r' + 'Time step ' + str(n + 1) + ' of '
                         + str(np.size(section.time.values)))

        # Select this time step.
        temp = section.sel(time=t)

        # Create the z coordinates for this time step using similar
        # assumptions as above, but for just this time step.
        if not flip:
            z_temp = np.linspace(z_bottom, temp['top'], num=n_laf)
        else:
            z_temp = np.linspace(temp['top'], z_bottom, num=n_laf)
        # Interpolate the profile to the average height spacing found above.
        temp.coords['z'] = ('LAF', z_temp)
        temp = temp.swap_dims({'LAF': 'z'}).interp(z=z_int)

        # Append this time step to the container for the section
        ds_to_label.append(temp)

    # Concatenate along the time dimension
    out = xr.concat(ds_to_label, 'time')
    return out
