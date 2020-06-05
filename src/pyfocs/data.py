# coding=utf-8
import dtscalibration
from dtscalibration.datastore_utils import suggest_cable_shift_double_ended, shift_double_ended, merge_double_ended
import numpy as np
import xarray as xr


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
    dstore = dtscalibration.DataStore(ds)
    # Rename 'x' to 'LAF'
    dstore = dstore.rename({'LAF': 'x'})
    # DataStore requires dimensions of 'x, time'
    for dvar in dstore.data_vars:
        try:
            dtscalibration.datastore_utils.check_dims(dstore, [dvar], correct_dims=('x', 'time'))
        except AssertionError as e:
            dstore[dvar] = dstore[dvar].T
        if dvar in varnames_conversion:
            dstore = dstore.rename({dvar: varnames_conversion[dvar]})

    # Impose a fake acquisition time here
    timeval = 1
    timeunit = 's'
    dt_np = np.timedelta64(timeval, timeunit)
    dstore['userAcquisitionTimeFW'] = (('time'), np.ones(len(dstore['time'])) * dt_np)

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
        if cal_lib[cl]['type'] == 'validation':
            continue
        pr = cal_lib[cl]['ref_sensor']

        LAF1 = np.min(cal_lib[cl]['LAF'])
        LAF2 = np.max(cal_lib[cl]['LAF'])

        # This section is poorly formatted
        if np.isnan(LAF1) or np.isnan(LAF2):
            continue
        # This section does not exist in the data.
        if LAF1 < ds.LAF.min().values or LAF2 > ds.LAF.max().values:
            continue
        ref_section = slice(LAF1, LAF2)

        sections[pr].append(ref_section)

    if not [ref for pr in sections for ref in sections[pr]]:
        mess = ('No reference sections of type calibration were found. '
                'Verify the calibration library.')
        print(mess)
        raise ValueError
    dstore.attrs.update(ds.attrs)
    dstore.sections = sections

    # Last missing attribute
    dstore.attrs['isDoubleEnded'] = '0'

    return dstore


def from_datastore(dstore,
                   datavars=None,
                   coords=None,
                   return_cal_params=False):
    '''
    Convert the pyfocs version of an xarray dataset into a dtscalibration
    datastore object.

    '''

    # Dictionaries for converting to dtscalibration names.
    # Reverse options untested at the moment.
    if dstore.attrs['isDoubleEnded'] == '1':
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
    ds_cal = xr.Dataset({'cal_temp': (('time', 'LAF'), dstore[varn_conv['cal_temp']].T)},
                        coords={'time': dstore.time,
                                'LAF': dstore.LAF,
                               })

    # Pass out the ammended attributes
    del dstore.attrs['_sections']
    ds_cal.attrs = dstore.attrs

    # Pass out all datavars (typically these are the reference probes)
    if datavars:
        for dvar in datavars:
            ds_cal[dvar] = dstore[dvar]

    # Keep the specified coordinates
    if coords:
        for c in coords:
            ds_cal.coords[c] = dstore.coords[c]

    return ds_cal


def double_end_dv_clean(ds):
    acc_datavars = ['st', 'ast', 'userAcquisitionTimeFW']
    drop_vars = [dv for dv in ds.data_vars if dv not in acc_datavars and 'x' in ds[dv].coords]

    return ds.drop(drop_vars)


def merge_single(dstore_fw, dstore_bw, shift_window=20, fixed_shift=None):
    '''
    Merge two single-ended channels to a single double-ended configuration.
    '''
    # Strip out any data values not expected by datastore
    dstore_fw = double_end_dv_clean(dstore_fw)
    dstore_bw = double_end_dv_clean(dstore_bw)

    # Assume that cable length is the same between the two channels.
    cable_length = dstore_fw.x.max().values

    double = merge_double_ended(
        ds_fw=dstore_fw,
        ds_bw=dstore_bw,
        cable_length=cable_length,
        plot_result=False
    )

    if not fixed_shift:
        shift_lims = np.arange(-shift_window, shift_window, 1, dtype=int)

        # Estimate the number of indices to shift the two channels to align them.
        # I use some overly generous limits for searching
        # This parameter should be made an optional argument passed to the function.
        shift1, shift2 = suggest_cable_shift_double_ended(
            double.mean(dim='time').compute(),
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
    else:
        shift = fixed_shift
    double = shift_double_ended(double, shift)

    return double


def single_calibrate(single, method):
    '''
    Calibrate a single-ended datastore object.
    '''

    if method == 'ols':
        # Calibrate
        single.calibration_single_ended(
            method=method,
        )
    elif method == 'wls':
        # Estimate the variances
        # @ allow for choosing between linear, exponential, and constant
        st_var, resid = single.variance_stokes_constant(st_label='st')
        ast_var, _ = single.variance_stokes_constant(st_label='ast')

        # And calibrate
        single.calibration_single_ended(
            st_var=st_var,
            ast_var=ast_var,
            method=method)

    # pyfocs is generally used to process shorter time chunks
    # Loading the dask object into memory can alleviate slow
    # operations later when the data is accessed.
    single.load()

    return single


def double_calibrate(double, method, keep_fw_bw=False):
    '''
    Calibrate a double-ended datastore object.
    '''

    if method == 'ols':
        # Calibrate
        double.calibration_double_ended(
            store_tmpw='tmpw',
            method=method,
        )
    elif method == 'wls':
        # Estimate the variances
        st_var, resid = double.variance_stokes(st_label='st')
        ast_var, _ = double.variance_stokes(st_label='ast')
        rst_var, _ = double.variance_stokes(st_label='rst')
        rast_var, _ = double.variance_stokes(st_label='rast')

        # And calibrate
        double.calibration_double_ended(
            st_var=st_var,
            ast_var=ast_var,
            rst_var=rst_var,
            rast_var=rast_var,
            store_tmpw='tmpw',
            method=method)

    # We only keep calibrated temperature
    # @ Add option for estimating noise variance and returning that value
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


def assign_ref_data(dstemp, cal, ref_data=None):
    # List of probe names
    probe_names = []

    # Add in external reference data. Interpolate to the DTS time.
    if cal['external_flag']:
        temp_ref_data = ref_data.reindex_like(dstemp.time,
                                              method='nearest')

        for ext_ref in cal['external_fields']:
            mess = ('{ef} was not found in the external reference data.')
            assert ext_ref in temp_ref_data, mess.format(ef=ext_ref)
            dstemp[ext_ref] = temp_ref_data[ext_ref]
            probe_names.append(ext_ref)

        # If the bath pt100s and dts do not line up in time,
        # notify the user.
        if not (np.size(np.flatnonzero(~np.isnan(dstemp.temp.values))) > 0):
            print('PT100 and DTS data do not line up in time for ' + raw_nc)

    # Rename built in probes if they are used.
    if cal['builtin_flag']:
        probe1 = cal['builtin_probe_names']['probe1Temperature']
        probe2 = cal['builtin_probe_names']['probe2Temperature']
        if probe1:
            dstemp = dstemp.rename({'probe1Temperature': probe1})
            probe_names.append(probe1)
        if probe2:
            dstemp = dstemp.rename({'probe2Temperature': probe2})
            probe_names.append(probe2)

    return dstemp, probe_names
