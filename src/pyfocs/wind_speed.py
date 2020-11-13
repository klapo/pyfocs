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


def calculate(
    heated,
    unheated,
    power,
    lw_in,
    method='S15',
    params=None,
    deltaT_dt_smoothing=None,
    return_EB=False):
    '''
    Inputs:
        heated - xarray object of dts observations. Assumed to have be in physical units, not LAF.
        unheated - xarray object of dts observations. Assumed to have be in physical units, not LAF.
        power - Power applied to the fiber in W/m. Must either be a single
            number or a DataArray with the same coordinates as heated and
            unheated.
        method - {['S15'], 'vR20'}, string indicating which method to use when
            calculating the convective heat coefficient. Choices are the
            original Sayde et al., 2015 and the van Ramshorst et al., 2020.
    Outputs:
        fiber_wind - an xarray dataset with the properties from subtracting
         `heated` and `unheated`.
    '''

    if 'S15' not in method and 'vR20' not in method:
        raise KeyError('method must be either S15 or vR20.')

    # Power should either be a float or an array that is a function of LAF
    # @ check that the above statement is True.

    # @ The params dictionary should not be optional -- there should be no
    # default values
    # Constants and coefficients
    if params is None:
        params = {}
    # Cable radius (in m)
    if 'rad' in params:
        rad = params['rad']
    else:
        rad = 0.000625
    # Heat capacity (in J kg^-1 K^-1)
    if 'crv' in params:
        crv = params['crv']
    else:
        crv = 1
    # Emissivity of the fiber surface (-)
    if 'emissivity' in params:
        emissivity = params['emissivity']
    else:
        emissivity = 0.95  # From table for PE
    # Fiber density (kg m^-3)
    if 'density' in params:
        density = params['density']
    else:
        density = 1000

   # Stefan-Boltzmann constant
    sigma = 5.67 * 10**(-8)

    # Coefficients for modeling the convective heat flux assuming that Re > 40
    if method == 'S15':
        c = 0.51
        m = 0.5
        npr = 0.37
    elif method == 'vR20':
        c = 0.683
        m = 0.466
        npr = 1 / 3

    # Basic plausibility check for units
    if heated.max() < 100:
        print('Fiber temperatures appear to be in Celcius, not Kelvin')
    if unheated.max() < 100:
        print('Fiber temperatures appear to be in Celcius, not Kelvin')

    ########
    # Diffeerence between the heated and unheated fibers.
    # Resampling/averaging/smoothing should be done by the user.
    # Future development could include using a smoother/more
    # sophisticated differencing estimate.
    delta = (heated - unheated)

    ########
    # Time-derivative
    # Determine the size of the time step from the data (assumes a regular interval).
    dt = heated.time.diff(dim='time').median(dim='time').values
    # Convert to seconds for SI.
    dt = pd.Timedelta(dt).total_seconds()
    # Difference between each time step for the provided dataset.
    deltaT_dt_heat = (heated.diff(dim='time')) / dt

    # Align the time derivative with the original time series.
    deltaT_dt_heat = deltaT_dt_heat.reindex_like(
        heated, method='nearest', tolerance=2 * dt)

    # This option "smooths" the noisy storage term. However, this term is orders of
    # magnitude smaller than all the others so it has a negligible effect. As such,
    # it is a candidate for removal in future versions.
    if deltaT_dt_smoothing is not None:
        deltaT_dt_heat = deltaT_dt_heat.resample(time=deltaT_dt_smoothing).mean()
        deltaT_dt_heat = deltaT_dt_heat.interp_like(heated)

    ########
    # Radiation components
    # By using identical fibers we can avoid some of the radiation modeling from S15.
    # We need to derive the incoming and outgoing longwave.
    lw_out_heat = sigma * emissivity * heated**4
    # @ Make the user do whatever treatment of the incoming longwave outside
    # pyfocs.
    lw_in = lw_in / 2 * emissivity

    ########
    # Temperature dependent properties of the convective heat transfer
    # Prandtl numbers for each fiber.
    prandtl_heat = prandtl(heated.values)
    prandtl_unheat = prandtl(unheated.values)

    # Thermal conductivity (k) and kinematic viscosity (nu) are calculated using
    # the "air" temperature, which we get from the unheated fiber.
    k = thermalConductivity(unheated.values)
    nu = kinematicViscosity(unheated.values)

    ########
    # Solving for wind speed
    # We are using a different approach than in Sayde et al., 2015. We are instead doing
    # the difference in the energy balance between the two fibers. This allows the radiation
    # terms to drop out, regardless of night/day time. The below solution, however, relies on
    # assuming the heat flux off of the unheated fibers is approximately zero. I'm not sure
    # if that is reasonable.
    if method == 'S15':
        # Use the coefficients for Re > 40 (defined above) to model the
        # convective heat flux. The original expression from the S15 matlab code
        # for reference:
        # fd=c*((2*rad/kv).^m)*(pr.^npr)*((pr/prs).^0.25)*k*3.14*DTch2ch1;
        demon = (c * (2 * rad)**(m - 1) * prandtl_unheat**npr
                 * (prandtl_unheat / prandtl_heat)**(0.25)
                 * k * nu**(-m) * delta)
    elif method == 'vR20':
        # This version assumes that prandtl_heat and prandtl_unheat are
        # identical as delta-T is never more than 6K. We find detla-T is on
        # average 8.8K and ranges to 31K for LOVE19. This still creates a
        # small difference (0.708 vs 0.711 assuming 60C and 30C) but the
        # assumption may need to be re-examined for certain applications.
        demon = (c * (2 * rad)**(m - 1) * prandtl_unheat**(npr)
                 * k * nu**(-m) * delta)

    # The original expression from the S15 code for reference
    # fn=crv*(((-dT1)/dt))+dT4*e*2*rad*3.14;
    dQdt = -density * crv * (rad / 2) * (deltaT_dt_heat)
    netLW = (lw_in - lw_out_heat)
    numer = (dQdt + netLW)

    # Original expression (pre 09.11.20)
    # fiber_wind = ((0.5 * power / (2 * np.pi * rad) + numer) / demon)**(1 / m)
    # I believe it should be
    fiber_wind = ((0.5 * power / (np.pi * rad) + numer) / demon)**(1 / m)
    if return_EB:
        return fiber_wind, dQdt, netLW, demon
    else:
        return fiber_wind
