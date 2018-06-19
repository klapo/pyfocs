######################
## Matrix Inversion ##
######################
def matrixInversion(dsCal, refField1, refField2, refField3, refLoc1, refLoc2, refLoc3, direction='forward'):
    refT1 = dsCal[refField1] + 273.15
    refT2 = dsCal[refField2] + 273.15
    refT3 = dsCal[refField3] + 273.15

    # For double ended calibration, the reverse pulse needs to have LAF flipped.
    if direction == 'reverse':
        dsCal['LAF'] = np.flip(dsCal.LAF.values, 0)

    section1 = dsCal.swap_dims({'LAF': 'location'}).sel(location=refLoc1)
    section2 = dsCal.swap_dims({'LAF': 'location'}).sel(location=refLoc2)
    section3 = dsCal.swap_dims({'LAF': 'location'}).sel(location=refLoc3)

    ref_z1 = section1.LAF.mean(dim='location')
    ref_z2 = section2.LAF.mean(dim='location')
    ref_z3 = section2.LAF.mean(dim='location')

    # Amplitudes of stokes/anti-stokes
    if 'logPsPas' in dsCal:
        stokesRatio1 = section1.logPsPas.mean(dim='location')
        stokesRatio2 = section2.logPsPas.mean(dim='location')
        stokesRatio3 = section3.logPsPas.mean(dim='location')
    else:
        stokesRatio1 = np.log(section1.Ps / section1.Pas).mean(dim='location')
        stokesRatio2 = np.log(section2.Ps / section2.Pas).mean(dim='location')
        stokesRatio3 = np.log(section3.Ps / section3.Pas).mean(dim='location')

    # Allocate the calibration variables and manual temperature
    gamma = np.ones(np.shape(dsCal.time.values)) * -9999.
    C = np.ones(np.shape(dsCal.time.values)) * -9999.
    delta_alpha = np.ones(np.shape(dsCal.time.values)) * -9999.
    manualTemp = np.ones(np.shape(dsCal.temp.values)) * -9999.

    # Within each time step solve for the calibration parameters
    for n, t in enumerate(dsCal.time):
        print('Time step ' + str(n + 1) + ' of ' + str(np.size(dsCal.time.values)), end="\r", flush=True)
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

        gamma[n] = x[0]
        C[n] = x[1]
        delta_alpha[n] = x[2]

        #############################
        ## Recalculate temperature ##
        #############################
        manualTemp[n] = gamma[n] / (np.log(dsCal.Ps.sel(time=t) / dsCal.Pas.sel(time=t))
                                    + C[n] - delta_alpha[n] * dsCal.LAF)


    # For double ended calibration, the reverse pulse needs to have LAF flipped.
    if direction == 'reverse':
        dsCal['LAF'] = np.flip(dsCal.LAF.values, 0)

    # Assign calibrated temperature to dataset
    dsCal['manualTemp'] = (('time', 'LAF'), manualTemp - 273.15)
    print('')
    print('Calibration done...')

    return(dsCal, gamma, C, delta_alpha)