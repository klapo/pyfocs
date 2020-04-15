import dtscalibration as dtscal
from dtscal.datastore_utils import suggest_cable_shift_double_ended, shift_double_ended, merge_double_ended


def to_datastore(ds, config, double):
    '''
    Convert the pyfocs version of an xarray dataset into a dtscalibration
    datastore object.
    '''

    # Dictionaries for converting to dtscalibration names.
    # Reverse options untested at the moment.
    varnames_conversion = {'Ps': 'st',
                           'Pas': 'ast',
                           'rPs': 'rst',
                           'rPas': 'rast',
                           'temp': 'tmp',
                          }

    # Convert to a DataStore object to take advantage of the dtscalibration methods
    dstore = dtscal.DataStore(ds)
    # Rename 'x' to 'LAF'
    dstore = dstore.rename({'LAF': 'x'})
    # DataStore requires dimensions of 'x, time'
    for dvar in dstore.data_vars:
        try:
            dtscal.datastore_utils.check_dims(dstore, [dvar], correct_dims=('x', 'time'))
        except AssertionError as e:
            dstore[dvar] = dstore[dvar].T
        if dvar in varnames_conversion:
            dstore = dstore.rename({dvar: varnames_conversion[dvar]})

    ds.attrs['cal_method'] = config['calibration']['method']

    # Impose a fake acquisition time here
    timeval = 0
    timeunit = 's'
    dt_np = np.timedelta64(timeval, timeunit)
    if double:
        dstore['userAcquisitionTimeFW'] = (('time'), np.ones(len(dstore['time'])) * dt_np)
    else:
        dstore['userAcquisitionTime'] = (('time'), np.ones(len(dstore['time'])) * dt_np)

    # Build the reference sections

    # Location library of reference sections
    cal_lib = config['calibration']['library']

    # Probe names
    probe_names = [cal_lib[cl]['ref_sensor'] for cl in cal_lib]
    probe_names = np.unique(probe_names)

    # Sections container for the DataStore Object. This object
    # is a collection of slices with the key given by the reference
    # sensor name.
    sections = {}
    for pr in probe_names:
        sections[pr] = []

    for cl in cal_lib:
        pr = cal_lib[cl]['ref_sensor']

        LAF1 = np.min(cal_lib[cl]['LAF'])
        LAF2 = np.max(cal_lib[cl]['LAF'])
        if np.isnan(LAF1) or np.isnan(LAF2):
            continue
        ref_section = slice(LAF1, LAF2)

        sections[pr].append(ref_section)

    dstore.sections = sections

    return dstore


def from_datastore(dstore,
                   double,
                   datavars,
                   coords,
                   return_cal_params=False):
    '''
    Convert the pyfocs version of an xarray dataset into a dtscalibration
    datastore object.

    '''

    # Dictionaries for converting to dtscalibration names.
    # Reverse options untested at the moment.
    if double:
        varn_conv = {
            'cal_temp': 'tmpw',
        }
    else:
        varn_conv = {
            'cal_temp': 'tmpf',
        }

# For later once I understand how the params are labeled.
#     if double and return_cal_params:
#         varnames_conversion.append()
#     elif not double and return_cal_params:
#         varnames_conversion.append()

    # Convert 'x' to 'LAF' for naming consistency
    dstore = dstore.rename({'x': 'LAF'})

    # Create the xr Dataset that will be returned.
    ds_cal = xr.Dataset({'cal_temp': (('LAF', 'time'), dstore[varn_conv['cal_temp']])},
                        coords={'time': dstore.time,
                                'LAF': dstore.LAF,
                               })


    # Pass out all datavars (typically these are the reference probes)
    for dvar in datavars:
        ds_cal[dvar] = dstore[dvar]

    # Pass out the ammended attributes
    ds_cal.attrs['cal_method'] = dstore.attrs['calibration_method']
    ds_cal.attrs['double_ended'] = double

    # Keep the specified coordinates
    for c in coords:
        ds_cal.coords[c] = dstore.coords[c]

    return ds_cal


def double_end_dv_clean(ds):
    acc_datavars = ['st', 'ast', 'userAcquisitionTimeFW']
    drop_vars = [dv for dv in ds.data_vars if dv not in acc_datavars and 'x' in ds[dv].coords]

    return ds.drop(drop_vars)


def merge_single(dstore_fw, dstore_bw, shift_window=20):
    '''
    Merge two single-ended channels to a single double-ended configuration.
    '''
    # Strip out any data values not expected by datastore
    dstore_fw = double_end_dv_clean(dstore_fw)
    dstore_bw = double_end_dv_clean(dstore_bw)

    # Assume that cable length is the same between the two channels.
    cable_length = dstore_fw.x.max().values

    # dtscal's plotting is broken due to some dependency stuff. Weirdly these exact lines work when I plot outside the dtscalibration package
    double = merge_double_ended(
        ds_fw=dstore_fw,
        ds_bw=dstore_bw,
        cable_length=cable_length,
        plot_result=False
    )

    shift_lims = np.arange(-shift_window, shift_window, 1, dtype=int)

    # Estimate the number of indices to shift the two channels to align them.
    # I use some overly generous limits for searching
    # This parameter should be made an optional argument passed to the function.
    shift1, shift2 = suggest_cable_shift_double_ended(
        double.isel(time=[0, -1]).compute(),
        shift_lims,
        plot_result=False,
    )

    # If the recommended shifts are identical
    if shift1 == shift2:
        shift = shift1
    # Otherwise use the smaller shift the two lengths should be close to each other.
    else:
        mess = ('When merging the single-ended observations into '
                'a double ended data set the optimal shift was '
                'found to have two values {s1} and {s2}.\nThe '
                'smallest value was autoselected.')
        print(mess.format(s1=shift1, s2=shift2))
        shift = np.min([shift1, shift2])
    double = shift_double_ended(double, shift)

    return double


def double_calibrate(double, method, keep_fw_bw=False):
    '''
    Calibrate a double-ended datastore object.
    '''

    # Calibrate
    double.calibration_double_ended(
        store_tmpw='tmpw',
        method=method,
    )

    # Since the confidence intervals, as described in example notebook 16,
    # are not relevant for atmospheric deployments, we drop everything except temperature.
    fields_to_drop = ['userAcquisitionTimeFW',
                      'userAcquisitionTimeBW',
                      'tmpf_mc_var',
                      'tmpb_mc_var',
                      'tmpw_mc_var',
                      'p_val',
                      'p_cov',
                      'st',
                      'rst',
                      'ast',
                      'rast',
                     ]

    # For internal testing it can be desirable to keep the
    # forward and backwards channels.
    if not keep_fw_bw:
        fields_to_drop.extend(
            ['tmpb',
             'tmpf',
            ]
        )

    # Drop the variables that conflict with pyfocs.
    double = double.drop_vars(fields_to_drop, errors='ignore')

    # pyfocs is generally used to process shorter time chunks
    # Loading the dask object into memory can alleviate slow
    # operations later when the data is accessed.
    double.load()

    return double
