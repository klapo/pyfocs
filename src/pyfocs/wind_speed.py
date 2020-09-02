import pandas as pd
import numpy as np


def prandtl(temp):
    '''
    Calculate Prandlt number of air (pr) from temperature (degree K)
    Source:  http://www.engineeringtoolbox.com/air-properties-d_156.html
    '''
    # Units check
    if (np.nanmax(temp) and np.nanmin(temp)) < 150:
        temp = temp + 273.15
        # print('Converted air temperature from Celsius to Kelvin.')

    temp = temp / 1000
    return(-0.1 * temp + 0.7423)


def kinematicViscosity(temp):
    '''
    Calculate the Kinematic viscosity of air (m2/s) from temperature (degree K)
    Source:  http://www.engineeringtoolbox.com/air-properties-d_156.html
    '''

    # Units check
    if (np.nanmax(temp) and np.nanmin(temp)) < 150:
        temp = temp + 273.15
        # print('Converted air temperature from Celsius to Kelvin.')

    temp = temp / 1000
    return((75 * temp**2 + 47.827 * temp - 5.3589) * 0.000001)


def thermalConductivity(temp):
    '''
    This function gives the thermal conductivity of air for a given temperature
        k = 1.5207 x 10^(-11) x T^3 - 4.8574 x 10^(-8) x T^2 + 1.0184 x 10^(-4) x T - 0.00039333
    Where k is in W/m/K and T is in K
    '''

    # Units check
    if (np.nanmax(temp) and np.nanmin(temp)) < 150:
        temp = temp + 273.15
        # print('Converted air temperature from Celsius to Kelvin.')

    k = 1.5207 * 10**(-11) * temp**3 - 4.8574 * 10**(-8) * temp**2 + 1.0184 * 10**(-4) * temp - 0.00039333
    return(k)


def calculate(heated, unheated, power):
    '''
    Inputs:
        heated - xarray object of dts observations. Assumed to have be in physical units, not LAF.
        unheated - xarray object of dts observations. Assumed to have be in physical units, not LAF.
        power - Power applied to the fiber in W/m. Must either be a single
            number or a DataArray with the same coordinates as heated and
            unheated.
    Outputs:
    '''

    # Power should either be a float or an array that is a function of LAF
    # @ check that the above statement is True.

    # Constants coefficients
    offset = 0  # unused
    # Cable radius in meters
    rad = 0.000625
    # Spatial (instrument) resolution of dts observations
    L = 0.127
    # Heat capacity (in ? units)
    crv = 1
    # Stefan-Boltzmann constant
    sigma = 5.67 * 10**(-8)
    # Emissivity of the fiber surface
    emissivity = 0.95  # From table for PE
    # Coefficients for modeling the convective heat flux assuming that Re > 40
    c = 0.51
    m = 0.5
    npr = 0.37

    ########
    # Diffeerence between the heated and unheated fibers.
    # mean(Heated(t) + Heated(t+1) - Unheated(t) - Unheated(t+1)) / 2
    # Future development could include using a smoother/more sophisticated derivative estimate.
    delta = (heated - unheated)
    # size of the dataset
    # len_time = dts.time.size
    # Resampling/averaging/smoothing should be done by the user.
    #     delta = (heated.isel(time=slice(0, len_time-1)).values
    #              + heated.isel(time=slice(1, len_time)).values
    #              - unheated.isel(time=slice(0, len_time-1)).values
    #              - unheated.isel(time=slice(1, len_time))).values

    ########
    # Time-derivative
    # Determine the size of the time step from the data (assumes a regular interval).
    dt = heated.time.diff(dim='time').median(dim='time').values
    # Convert to seconds for SI.
    dt = pd.Timedelta(dt).total_seconds()
    # Difference between each time step for the provided dataset.
    deltaT_dt_heat = (heated.diff(dim='time')) / dt
    deltaT_dt_unheat = (unheated.diff(dim='time')) / dt

    ########
    # Align the heated/unheated datasets with the differenced datasets
    heated = heated.reindex_like(heated.diff(dim='time'))
    unheated = unheated.reindex_like(unheated.diff(dim='time'))
    delta = delta.reindex_like(delta.diff(dim='time'))

    ########
    # Radiation components
    # By using identical fibers we can avoid most of the radiation modeling.
    # Only component that is different between the heated/unheated fibers is
    # the longwave cooling.
    # @ Check for units
    lw_out_heat = sigma * emissivity * heated.values**4
    lw_out_unheat = sigma * emissivity * unheated.values**4

    ########
    # Temperature dependent properties of the convective heat transfer
    # Prandtl numbers for each fiber.
    prandtl_heat = prandtl(heated.values)
    prandtl_unheat = prandtl(unheated.values)

    # Thermal conductivity (k) and kinematic viscosity (nu) are calculated using the "air" temperature, which we get from the unheated fiber.
    k = thermalConductivity(unheated.values)
    nu = kinematicViscosity(unheated.values)

    ########
    # Solving for wind speed
    # We are using a different approach than in Sayde et al., 2015. We are instead doing
    # the difference in the energy balance between the two fibers. This allows the radiation
    # terms to drop out, regardless of night/day time. The below solution, however, relies on
    # assuming the heat flux off of the unheated fibers is approximately zero. I'm not sure
    # if that is reasonable.
    # Use the coefficients for Re > 40 (defined above) to model the convective heat flux.
    # The original expression for reference
    # fd=c*((2*rad/kv).^m)*(pr.^npr)*((pr/prs).^0.25)*k*3.14*DTch2ch1;
    demon = c * (2 * rad)**(m - 1) * prandtl_unheat**npr * (prandtl_unheat / prandtl_heat)**(0.25) * k * nu**(-m) * delta

    # The original expression for reference
    # fn=crv*(((-dT1)/dt))+dT4*e*2*rad*3.14;
    numer = -1000 * crv * (rad / 2) * (deltaT_dt_heat - deltaT_dt_unheat) + (lw_out_heat - lw_out_unheat)

    fiber_wind = ((0.5*power / (2 * np.pi * rad) + numer) / demon)**(1 / m)
    return fiber_wind
