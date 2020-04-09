#!/usr/bin/env python
import numpy as np
from datetime import timedelta
import pandas as pd
import xarray as xr

# OS interaction
import os
import yaml
import tkinter as tk  # For the dialog that lets you choose your config file
from tkinter import filedialog
import copy
import csv
import sys

# UBT's package for handling dts data
import pyfocs

# Ignore the future compatibility warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)

# -----------------------------------------------------------------------------
# Configuration file
# -----------------------------------------------------------------------------
# Test if a file was provided to the script from the terminal.
try:
    filename_configfile = sys.argv[1]
    if not os.path.exists(filename_configfile):
        print('Config file' + sys.argv[1] + 'not found.')

# Prompt the user for the config file.
except IndexError:
    # Open the config file
    root = tk.Tk()
    # to close this stupid root window that doesn't let itself be closed.
    root.withdraw()
    ftypes = [
        ('yaml configuration file', '*.yml'),
        ]
    filename_configfile = filedialog.askopenfilename(filetypes=ftypes)

# Call the config file reader and formatter.
internal_config, lib = pyfocs.check.config(filename_configfile)

if internal_config is None:
    raise NameError('Configuration file generated an error.')

# Unpack items from the formatted configuration dictionaries.
# Processing flags
write_mode = internal_config['write_mode']
ref_temp_option = internal_config['flags']['ref_temp_option']
archiving_flag = internal_config['flags']['archiving_flag']
archive_read_flag = internal_config['flags']['archive_read_flag']
calibrate_flag = internal_config['flags']['calibrate_flag']
final_flag = internal_config['flags']['final_flag']

# Step loss corrections
step_loss_flag = internal_config['step_loss']['flag']
if step_loss_flag:
    splice_LAF = internal_config['step_loss']['LAF']
    step_loss_corrections = internal_config['step_loss']['correction']

# Channels/experiments/output file names
channelNames = internal_config['channelNames']
experiment_names = internal_config['experiment_names']
outname_suffix = internal_config['outname_suffix']

# -----------------------------------------------------------------------------
# Archive and create raw netcdfs
# -----------------------------------------------------------------------------
#%% Archive/read the raw xml files
for exp_name in experiment_names:
    print('-------------')
    print(exp_name)
    print('-------------')

    # archiving
    if archiving_flag:
        print('-------------')
        print('Archiving raw xml files.')
        print(' ')
        print('archiving ', exp_name)
        pyfocs.archiver(internal_config[exp_name])

    # write raw netCDF
    if archive_read_flag:
        print('-------------')
        print('Writing netcdfs from raw xml files.')
        print(' ')
        print('creating raw netcdf for experiment: ', exp_name)
        pyfocs.archive_read(internal_config[exp_name],
                            write_mode=write_mode)

for exp_name in experiment_names:
    # --------------------------------------------------------------------------
    # Process DTS files
    # --------------------------------------------------------------------------
    if calibrate_flag:

        # Grab the options detailing the calibration.
        cal = internal_config['calibration']

        # If we are doing a double ended calibration method we will need to
        # keep track of which files have been finished.
        double_ended_methods = ['ols double', 'wls double', 'temp_matching']
        if cal['method'] in double_ended_methods:
            finished_files = []

        # Get the external data if it was indicated.
        if cal['external_flag']:
            # Get the metdata
            ref_data = xr.open_dataset(internal_config['external_data'])
            ref_data = ref_data.resample(time=delta_t).interpolate('linear')

        # Find all 'raw' netcdfs within the processed directory,
        # sort them (by date), and process each individually.
        os.chdir(internal_config[exp_name]['directories']['dirRawNetcdf'])
        contents = os.listdir()
        ncfiles = [file for file in contents
                   if '.nc' in file
                   and 'raw' in file
                   and any(chN in file for chN in channelNames)]
        ncfiles.sort()
        ntot = np.size(ncfiles)
        for nraw, raw_nc in enumerate(ncfiles):
            outname_channel = raw_nc.split('.')[0].split('_')[-2]
            outname_date = raw_nc.split('.')[0].split('_')[-1]

            # Check of the files to create and determine if we are overwriting.
            if write_mode == 'preserve':
                dir_cal = internal_config[exp_name]['directories']['dirCalibrated']

                outname = '_'.join(filter(None, [exp_name, 'cal',
                                                 outname_channel,
                                                 outname_date,
                                                 outname_suffix,
                                                 ])) + '.nc'
                if os.path.exists(os.path.join(dir_cal, outname)):
                    skip_flag = True
                    # Determine if we are preserving or overwriting the
                    # existing files.
                    print('Preserving existing netcdf ' + outname
                          + ' (' + str(nraw + 1)
                          + ' of ' + str(ntot) + ')')
                else:
                    skip_flag = False

                # The calibrated file exists, do not process and calibrate.
                if skip_flag:
                    continue
            # Either the file does not exist or it is to be overwritten.
            print('Processing ' + outname_date + ' (' + str(nraw + 1)
                  + ' of ' + str(ntot) + ')')

            # Open up raw file for this interval.
            os.chdir(internal_config[exp_name]['directories']['dirRawNetcdf'])
            dstemp = xr.open_dataset(raw_nc)

            # Add in external reference data. Interpolate to the DTS time.
            if cal['external_flag']:
                temp_ref_data = ref_data.reindex_like(dstemp.time,
                                                      method='nearest')

                for ext_ref in cal['external_fields']:
                    dstemp[ext_ref] = temp_ref_data[ext_ref]

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
                if probe2:
                    dstemp = dstemp.rename({'probe2Temperature': probe2})

            # Step loss corrections if they are provided.
            if step_loss_flag:
                # Calculate the log power of the stokes/anti-stokes scattering
                dstemp['logPsPas'] = np.log(dstemp.Ps / dstemp.Pas)
                # Correct for any step-losses due to splicing.
                for spl_num, spl_LAF in enumerate(splice_LAF):
                    dstemp['logPsPas'] = dstemp.logPsPas.where((dstemp.LAF < spl_LAF),
                                                               dstemp.logPsPas + step_loss_corrections[spl_num])

            # Drop the unnecessary negative LAF indices
            LAFmin = internal_config['min_fiber_limit']
            LAFmax = internal_config['max_fiber_limit']
            dstemp = dstemp.sel(LAF=slice(LAFmin, LAFmax))
            dstemp.attrs['LAF_beg'] = dstemp.LAF.values[0]

            # Label the calibration locations
            dstemp = pyfocs.labelLoc_additional(dstemp,
                                                cal['library'],
                                                'calibration')

            # Calibrate the temperatures.
            if cal['method'] == 'single':
                dstemp, _, _, _ = pyfocs.matrixInversion(dstemp, cal)

            # Rename the instrument reported temperature
            dstemp = dstemp.rename({'temp': 'instr_temp'})

            # Output the calibrated dataset
            outname = '_'.join(filter(None, [exp_name, 'cal',
                                             outname_channel,
                                             outname_date,
                                             outname_suffix])) + '.nc'

            os.chdir(internal_config[exp_name]['directories']['dirCalibrated'])
            dstemp.to_netcdf(outname, engine='netcdf4')

            print('')

# -----------------------------------------------------------------------------
# Finalize with physical coordinates
# -----------------------------------------------------------------------------
#%% Final preparation of the DTS dataset

if final_flag:
    print('-------------')
    print('Final preparation: assigning physical coordinates, dropping unused fields and locations')
    print(' ')

    phys_locs = internal_config['phys_locs']

    # Time step for resampling to a uniform time step
    delta_t = internal_config['resampling_time']

    # When finalizing the dataset all extraneous coordinates and data
    # is dropped, leaving behind these variables.
    coords_to_keep = ['xyz', 'time', 'x', 'y', 'z', 'LAF']
    vars_to_keep = ['cal_temp']

    for exp_name in experiment_names:
        # Find all 'calibrated' netcdfs within the calibrated directory,
        # sort them by date, suffix, and process each individually.

        # @ONLY LOOK FOR THE SUFFIX ID FROM THE CONFIG.
        os.chdir(internal_config[exp_name]['directories']['dirCalibrated'])
        contents = os.listdir()
        ncfiles = [file for file in contents if '.nc' in file and 'cal' in file]
        ncfiles.sort()
        ntot = np.size(ncfiles)

        for ncal, cal_nc in enumerate(ncfiles):
            # Name of the output file for this archive chunk
            name_components = cal_nc.split('.')[0].split('_')
            # As we allow '_' in the experiment name/suffix we can't reliably
            # count elemnts of the list. To get around this, we look for a
            # reliably known element (the 'cal' indicator) and base our
            # string indexing on the location of that element.
            cal_str_ind = name_components.index('cal')
            outname_date = name_components[cal_str_ind + 2]

            print('Physically labeling ' + cal_nc + ' (' + str(ncal + 1)
                  + ' of ' + str(ntot) + ')')

            # Open each calibrated file.
            os.chdir(internal_config[exp_name]['directories']['dirCalibrated'])
            dstemp = xr.open_dataset(cal_nc)

            # Resample to a common time stamp interval with the reference/bath
            # instruments. Do this using a linear interpolation.
            # Round to the nearest time interval
            dt_start = pd.Timestamp(dstemp.time.values[0]).round(delta_t)
            dt_end = pd.Timestamp(dstemp.time.values[-1]).round(delta_t)

            # Catch a weird edge case for a single time step.
            if not dt_start == dt_end:
                # Create a regular interval time stamp index
                reg_time_pd = pd.date_range(start=dt_start,
                                            end=dt_end,
                                            freq=delta_t)
                dstemp = dstemp.interp(time=reg_time_pd,
                                       kwargs={'fill_value': 'extrapolate'})

            # Add a delta t attribute
            dstemp.attrs['dt'] = delta_t

            # Clean up unused variables and labels.
            vars_to_drop = [v for v in dstemp.data_vars
                            if v not in vars_to_keep]
            coords_to_drop = [c for c in dstemp.coords
                              if c not in coords_to_keep]

            # Clean up attributes and dropped the unused ones.
            dt = dstemp.attrs['dt']
            dLAF = dstemp.attrs['dLAF']
            dstemp = dstemp.drop(vars_to_drop).drop(coords_to_drop)
            dstemp.attrs = []
            dstemp.attrs['dt']  = dt
            dstemp.attrs['dLAF'] = dLAF

            for ploc in phys_locs:
                # Relabel the locations. This allows locations to
                # change after calibrating, as the calibration only
                # cares about the location of the reference baths.
                dstemp_ploc = pyfocs.labelLoc_additional(dstemp,
                                                        lib[ploc],
                                                        ploc)

                # Assign physical labels
                dstemp_ploc = pyfocs.labeler.dtsPhysicalCoords_3d(dstemp_ploc,
                                                                  lib[ploc])

                # Output each location type as a separate final file.
                outname = '_'.join(filter(None, [exp_name, 'final',
                                                 outname_date,
                                                 outname_suffix,
                                                 ploc])) + '.nc'
                os.chdir(internal_config[exp_name]['directories']['dirFinal'])
                dstemp_ploc.to_netcdf(outname, mode='w')

            # Make sure we don't reprocess files.
            finished_files.extend(nc_cal)
