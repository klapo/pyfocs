import dtscalibration as dtscal


def to_datastore(ds, cfg, ref):
    '''
    Convert the pyfocs version of an xarray dataset into a dtscalibration
    datastore object.
    '''

    phys_locs = cfg['location_library']

    # Dictionaries for converting to dtscalibration names.
    # Reverse options untested at the moment.
    varnames_conversion = {'Ps': 'st',
                           'Pas': 'ast',
                           'rPs': 'rst',
                           'rPas': 'rast',
                           'temp': 'TMP',
                           }

    # Convert to a DataStore object to take advantage of the
    # dtscalibration methods
    ds = dtscal.DataStore(ds)
    # Rename 'x' to 'LAF' as this dimension does not allow flexible names.
    ds = ds.rename({'LAF': 'x'})

    for dvar in ds.data_vars:
        # DataStore requires dimensions of 'x, time'
        try:
            dtscal.datastore.check_dims(ds, [dvar], correct_dims=('x', 'time'))
        except AssertionError:
            ds[dvar] = ds[dvar].T

        # Convert variable names to those dtscalibration expects.
        if dvar in varnames_conversion:
            ds = ds.rename({dvar: varnames_conversion[dvar]})

    probe_sections = {}
    # Insert the "sections" variable into the DataStore object.
    for fields in cfg['dataProperties']:
        if 'probe' in fields and 'Temperature' in fields:
            # Assign the reference variable to the datastore object.
            fn = cfg['dataProperties'][fields]
            ds[fn] = ref[fn]

            probe_sections[fn] = []
            # Iterate through the calibration reference sections
            for refsects in cfg['calibration']:
                if 'refLoc' not in refsects:
                    continue
                pl = phys_locs[cfg['calibration'][refsects]]['LAF']
                probe_sections[fn].append(slice(pl[0], pl[-1]))

    # Now do the datastore sections assignment
    ds.sections = probe_sections

    return ds
