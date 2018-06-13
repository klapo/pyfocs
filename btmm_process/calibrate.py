import os
import xarray as xr
import numpy as np

def prepCalibrate(cfg):
    '''
    Prep the data for being calibrated!
    '''
    ####
    # Unpack the config file
    dirProcessed = cfg['directories']['dirProcessed']
    filePrefix = cfg['fileName']['filePrefix']
    fileSuffix = cfg['fileName']['fileSuffix']

    ####
    # Prep the data for calibration
    os.chdir(dirProcessed)

    # Find all netcdf files that match the experiment suffix and prefix. Load them
    # into a xarray Dataset.
    ncFiles = [nc for nc in os.listdir() if filePrefix in nc and fileSuffix in nc and '.nc' in nc]
    ncFiles.sort()
    ds = xr.open_mfdataset(ncFiles, chunks=cfg['dataProperties']['chunkSize'])

    # Call the matrix solver function.
    ds = matrixInversion(ds, cfg['calibration'])

    os.chdir(dirProcessed)
    ds.to_netcdf(filePrefix + '_calibrated_' + fileSuffix  + '.nc', 'w')
    return ds

def matrixInversion(ds, refConfig):
    '''
    Solve for the calibration parameters
    '''
    # Unpack the configuration file
    refField1 = refConfig['refField1']
    refField3 = refConfig['refField2']
    refField2 = refConfig['refField3']
    refLoc1 = refConfig['refLoc1']
    refLoc2 = refConfig['refLoc2']
    refLoc3 = refConfig['refLoc3']

    # Reference temperatures
    refT1 = ds[refField1] + 273.15
    refT2 = ds[refField2] + 273.15
    refT3 = ds[refField3] + 273.15

    section1 = ds.swap_dims({'LAF': 'location'}).sel(location=refLoc1)
    section2 = ds.swap_dims({'LAF': 'location'}).sel(location=refLoc2)
    section3 = ds.swap_dims({'LAF': 'location'}).sel(location=refLoc3)

    ref_z1 = section1.LAF.mean(dim='location')
    ref_z2 = section2.LAF.mean(dim='location')
    ref_z3 = section3.LAF.mean(dim='location')

    # Amplitudes of stokes/anti-stokes
    stokesRatio1 = np.log(section1.Ps / section1.Pas).mean(dim='location')
    stokesRatio2 = np.log(section2.Ps / section2.Pas).mean(dim='location')
    stokesRatio3 = np.log(section3.Ps / section3.Pas).mean(dim='location')

    # Allocate the calibration variables and manual temperature
    gamma = np.ones(np.shape(ds.time.values)) * -9999.
    C = np.ones(np.shape(ds.time.values)) * -9999.
    delta_alpha = np.ones(np.shape(ds.time.values)) * -9999.
    manualTemp = np.ones(np.shape(ds.temp.values)) * -9999.

    # Within each time step solve for the calibration parameters
    for n, t in enumerate(ds.time):
        print('Time step ' + str(n) + ' of ' + str(np.size(ds.time.values)), end="\r", flush=True)


        # A matrix
        A = [[1, -refT1.sel(time=t), refT1.sel(time=t) * ref_z1],
             [1, -refT2.sel(time=t), refT2.sel(time=t) * ref_z2],
             [1, -refT3.sel(time=t), refT3.sel(time=t) * ref_z3],
            ]

        # b matrix
        b = [[refT1.sel(time=t) * stokesRatio1.sel(time=t)],
             [refT2.sel(time=t) * stokesRatio2.sel(time=t)],
             [refT3.sel(time=t) * stokesRatio3.sel(time=t)],
            ]

        x = np.linalg.solve(A, b)

        gamma[n] = x[0]
        C[n] = x[1]
        delta_alpha[n] = x[2]

        #############################
        ## Recalculate temperature ##
        #############################
        manualTemp[n] = gamma[n] / (np.log(ds.Ps.sel(time=t) / ds.Pas.sel(time=t)) + C[n] - delta_alpha[n] * ds.LAF)

    # Assign calibrated temperature to dataset
    ds['manualTemp'] = (('time', 'LAF'), manualTemp - 273.15)

    return(ds)

        # ##############################
        # ## Plot Calibration Results ##
        # ##############################
        # plt.plot(gamma)
        # plt.gca().set_ylabel('$\gamma$')
        #
        # plt.figure()
        # plt.plot(C)
        # plt.gca().set_ylabel('$C(t)$')
        #
        # plt.figure()
        # plt.plot(delta_alpha)
        # plt.gca().set_ylabel(r'$\Delta  \alpha (z)$')
