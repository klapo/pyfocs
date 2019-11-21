import numpy as np


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
