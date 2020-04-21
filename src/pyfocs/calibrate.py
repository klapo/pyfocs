import os
import xarray as xr
import numpy as np
import sys


def matrixInversion(dsCal, cfg):
    refField1 = cfg['refField1']
    refField2 = cfg['refField2']
    refField3 = cfg['refField3']

    refLoc1 = cfg['refLoc1']
    refLoc2 = cfg['refLoc2']
    refLoc3 = cfg['refLoc3']

    # Assume that the PT100 data is in Celsius
    refT1 = dsCal[refField1] + 273.15
    refT2 = dsCal[refField2] + 273.15
    refT3 = dsCal[refField3] + 273.15

    section1 = dsCal.swap_dims({'LAF': 'calibration'}).sel(calibration=refLoc1)
    section2 = dsCal.swap_dims({'LAF': 'calibration'}).sel(calibration=refLoc2)
    section3 = dsCal.swap_dims({'LAF': 'calibration'}).sel(calibration=refLoc3)

    ref_z1 = section1.LAF.mean(dim='calibration')
    ref_z2 = section2.LAF.mean(dim='calibration')
    ref_z3 = section3.LAF.mean(dim='calibration')

    # Amplitudes of stokes/anti-stokes
    if 'logPsPas' in dsCal:
        stokesRatio1 = section1.logPsPas.mean(dim='calibration')
        stokesRatio2 = section2.logPsPas.mean(dim='calibration')
        stokesRatio3 = section3.logPsPas.mean(dim='calibration')
    else:
        dsCal['logPsPas'] = np.log(dsCal.Ps / dsCal.Pas)
        stokesRatio1 = np.log(section1.Ps / section1.Pas).mean(dim='calibration')
        stokesRatio2 = np.log(section2.Ps / section2.Pas).mean(dim='calibration')
        stokesRatio3 = np.log(section3.Ps / section3.Pas).mean(dim='calibration')

    # Allocate the calibration variables and manual temperature
    gamma = np.ones(np.shape(dsCal.time.values)) * -9999.
    C = np.ones(np.shape(dsCal.time.values)) * -9999.
    delta_alpha = np.ones(np.shape(dsCal.time.values)) * -9999.
    manualTemp = np.ones(np.shape(dsCal.temp.values)) * -9999.

    # Within each time step solve for the calibration parameters
    for n, t in enumerate(dsCal.time):
        sys.stdout.write('\r' + 'Time step ' + str(n + 1) + ' of '
                         + str(np.size(dsCal.time.values)))
        T1 = refT1.sel(time=t).values
        T2 = refT2.sel(time=t).values
        T3 = refT3.sel(time=t).values
        sR1 = stokesRatio1.sel(time=t).values
        sR2 = stokesRatio2.sel(time=t).values
        sR3 = stokesRatio3.sel(time=t).values

        # A matrix
        A = [[1, -T1, T1 * ref_z1],
             [1, -T2, T2 * ref_z2],
             [1, -T3, T3 * ref_z3],
             ]

        # b matrix
        b = [[T1 * sR1],
             [T2 * sR2],
             [T3 * sR3],
             ]

        x = np.linalg.solve(A, b)

        # Store the calibration coefficients
        gamma[n] = x[0]
        C[n] = x[1]
        delta_alpha[n] = x[2]

        # Recalculate temperature for this time step
        manualTemp[n] = gamma[n] / (dsCal.logPsPas.sel(time=t)
                                    + C[n] - delta_alpha[n] * dsCal.LAF)

    # Assign calibrated temperature to dataset
    dsCal['cal_temp'] = (('time', 'LAF'), manualTemp - 273.15)
    print('')
    print('Calibration done...')

    return(dsCal, gamma, C, delta_alpha)
