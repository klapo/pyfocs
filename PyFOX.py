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

# Cores (mostly relevant to multicore fibers)
coretype = internal_config['coretype']
cores = internal_config['cores']

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

        # Time step for resampling to a uniform time step
        delta_t = internal_config['resampling_time']

        # Get the external data to add if we are not using the instrument PT100s
        if internal_config['flags']['ref_temp_option'] == 'external':
            # Get the metdata
            ref_data = xr.open_dataset(internal_config['external_data'])
            ref_data = ref_data.resample(time=delta_t).interpolate('linear')

            probe1_name = internal_config['probe1']
            probe2_name = internal_config['probe2']

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
                # Output file name creation.
                if coretype == 'multicore':
                    for c in cores:
                        outname_core = c
                        outname = '_'.join(filter(None, [exp_name, 'cal',
                                                         outname_channel,
                                                         outname_date,
                                                         outname_suffix,
                                                         outname_core])) + '.nc'
                        # If the file exists indicate that we are to skip
                        # this file.
                        if os.path.exists(os.path.join(dir_cal, outname)):
                            skip_flag = True
                            # Determine if we are preserving or overwriting the
                            # existing files.
                            print('Preserving existing netcdf ' + outname
                                  + ' (' + str(nraw + 1)
                                  + ' of ' + str(ntot) + ')')
                        # If one of the core files does not exist, exit the
                        # loop and continue.
                        else:
                            skip_flag = False
                            continue

                elif coretype == 'singlecore':
                    outname_core = None
                    outname = '_'.join(filter(None, [exp_name, 'cal',
                                                     outname_channel,
                                                     outname_date,
                                                     outname_suffix,
                                                     outname_core])) + '.nc'
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

            # Add in external reference data here
            if internal_config['flags']['ref_temp_option'] == 'external':
                temp_ref_data = ref_data.reindex_like(dstemp.time,
                                                      method='nearest',
                                                      tolerance=1)

                dstemp[probe1_name] = temp_ref_data[probe1_name]
                dstemp[probe2_name] = temp_ref_data[probe2_name]

                # Add additional external data for this data stream.
                if internal_config['external_fields']:
                    for ext_dat in internal_config['external_fields']:
                        dstemp[ext_dat] = temp_ref_data[ext_dat]

                # If the bath pt100s and dts do not line up in time,
                # notify the user.
                if not (np.size(np.flatnonzero(~np.isnan(dstemp.temp.values))) > 0):
                    print('PT100 and DTS data do not line up in time for ' + raw_nc)

            # Step loss corrections if they are provided.
            if step_loss_flag:
                # Calculate the log power of the stokes/anti-stokes scattering
                dstemp['logPsPas'] = np.log(dstemp.Ps / dstemp.Pas)
                # Correct for any step-losses due to splicing.
                for spl_num, spl_LAF in enumerate(splice_LAF):
                    dstemp['logPsPas'] = dstemp.logPsPas.where((dstemp.LAF < spl_LAF),
                                                               dstemp.logPsPas + step_loss_corrections[spl_num])

            # Split up multicore data
            if coretype == 'multicore':
                for c in cores:
                    dstemp_core = dstemp

                    # Drop the unnecessary negative LAF indices
                    core_LAF_start = cores[c]['LAF'][0]
                    core_LAF_end = cores[c]['LAF'][-1]
                    dstemp_core = dstemp_core.sel(LAF=(slice(core_LAF_start, core_LAF_end)))
                    dstemp_core.attrs['LAF_beg'] = core_LAF_start
                    dstemp_core.attrs['LAF_end'] = core_LAF_end

                    # Location labels
                    for loc_type_cur in lib[c]:
                        dstemp_core = pyfocs.labelLoc_additional(dstemp_core,
                                                                 lib[c][loc_type_cur],
                                                                 loc_type_cur)
                    # Converting to a netcdf ruins this step unfortunately.
                    # dstemp_core.attrs['loc_general_long'] = dict((l, config_user['loc_general'][l]['long name'])

                    dstemp_core, _, _, _ = pyfocs.matrixInversion(dstemp_core,
                                                                  internal_config)

                    # Rename the instrument reported temperature field
                    dstemp_core = dstemp_core.rename({'temp': 'instr_temp'})
                    dstemp_core.coords['core'] = c

                    # Output the calibrated dataset
                    outname_core = c
                    outname = '_'.join(filter(None, [exp_name, 'cal',
                                                     outname_channel,
                                                     outname_date,
                                                     outname_suffix,
                                                     outname_core])) + '.nc'
                    os.chdir(internal_config[exp_name]['directories']['dirCalibrated'])
                    dstemp_core.to_netcdf(outname, engine='netcdf4')

            elif coretype == 'singlecore':
                # Drop the unnecessary negative LAF indices
                dstemp = dstemp.sel(LAF=slice(0, None))
                dstemp.attrs['LAF_beg'] = dstemp.LAF.values[0]

                # Location labels
                for loc_type_cur in lib:
                    # lib_laf = {l: lib[loc_type_cur][l]['LAF'] for l in lib[loc_type_cur]}
                    dstemp = pyfocs.labelLoc_additional(dstemp,
                                                        lib[loc_type_cur],
                                                        loc_type_cur)

                # Calibrate the temperatures.
                dstemp, _, _, _ = pyfocs.matrixInversion(dstemp, internal_config)
                dstemp = dstemp.rename({'temp': 'instr_temp'})

                # Output the calibrated dataset
                outname_core = None
                outname = '_'.join(filter(None, [exp_name, 'cal',
                                                 outname_channel,
                                                 outname_date,
                                                 outname_suffix,
                                                 outname_core])) + '.nc'

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

    # List of finished files so we can skip files when handling multicore fibers.
    finished_files = []

    # When finalizing the dataset all extraneous coordinates and data
    # is dropped, leaving behind these variables.
    coords_to_keep = ['xyz', 'time', 'x', 'y', 'z', 'core', 'LAF']
    vars_to_keep = ['cal_temp']
    if coretype == 'multicore':
        cores_to_proc = list(cores.keys())

    for exp_name in experiment_names:
        # Find all 'calibrated' netcdfs within the calibrated directory,
        # sort them by date and core name, and process each individually.
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

            # Deal with multicore fibers differently than single core fibers.
            if coretype == 'multicore':
                os.chdir(internal_config[exp_name]['directories']['dirCalibrated'])
                # Check if we already handled this core.
                if cal_nc in finished_files:
                    continue
                # Cores are separate files but share all other name components.
                # For a single datetime find the cores available.
                nc_cal_core_name = '_'.join(name_components[:-1])

                # Finds all core files for this time step.
                nc_cal_core = [file for file in ncfiles
                               if nc_cal_core_name in file and 'cal' in file]
                dstemp_out = {}

                # Assign physical coordinates. Each core is appended to a
                # list to be merged later.
                for c in nc_cal_core:
                    dstemp = xr.open_dataset(c)
                    dstemp.load()

                    # Drop everything except the calibrated temperature.
                    # Clean up unused variables and labels.
                    vars_to_drop = [v for v in dstemp.data_vars
                                    if v not in vars_to_keep]

                    coords_to_drop = [c for c in dstemp.coords
                                      if c not in coords_to_keep]
                    try:
                        coords_to_drop.remove(phys_locs)
                    except ValueError:
                        # Just continue onwards. The physical location was
                        # not in the list of coordinates.
                        coords_to_drop = coords_to_drop
                    dstemp = dstemp.drop(vars_to_drop).drop(coords_to_drop)

                    # Reformat the config locations to specify just a
                    # single core
                    core = str(dstemp.core.values)
                    # Make sure we intend to process this core
                    if core not in cores_to_proc:
                        continue

                    for ploc in phys_locs:
                        if ploc not in dstemp_out:
                            dstemp_out[ploc] = []
                        # Relabel the locations. This allows locations to
                        # change after calibrating, as the calibration only
                        # cares about the location of the reference baths.
                        dstemp = pyfocs.labelLoc_additional(dstemp,
                                                            lib[core][ploc],
                                                            ploc)

                        # Assign physical labels
                        dstemp_out[ploc].append(
                            pyfocs.labeler.dtsPhysicalCoords_3d(dstemp,
                                                                lib[core][ploc]))

                # Merge the cores
                for ploc in phys_locs:
                    interp_to = None
                    for nc, c in enumerate(dstemp_out[ploc]):
                        # Make the assumption that the specific core is
                        # not important as they should be identical within
                        # the physical resolution of the DTS.
                        if not interp_to and c.cal_temp.any():
                            interp_to = c
                        elif c.cal_temp.any() and interp_to:
                            dstemp_out[ploc][nc] = dstemp_out[ploc][nc].interp(xyz=interp_to.xyz)

                    dstemp_out[ploc] = xr.concat(dstemp_out[ploc], dim='core')

                    # Clean up attributes and dropped the unused ones.
                    dt = dstemp_out[ploc].attrs['dt']
                    dLAF = dstemp_out[ploc].attrs['dLAF']
                    dstemp_out[ploc].attrs = []
                    dstemp_out[ploc].attrs['dt'] = dt
                    dstemp_out[ploc].attrs['dLAF'] = dLAF

                    # Output each location type as a separate final file.
                    outname = '_'.join(filter(None, [exp_name, 'final',
                                                     outname_date,
                                                     outname_suffix,
                                                     ploc])) + '.nc'
                    os.chdir(internal_config[exp_name]['directories']['dirFinal'])
                    dstemp_out[ploc].to_netcdf(outname, mode='w')

                # Make sure we don't reprocess files.
                finished_files.extend(nc_cal_core)

            # Deal with the single core fiber.
            else:
                # Open each calibrated file.
                os.chdir(internal_config[exp_name]['directories']['dirCalibrated'])
                dstemp = xr.open_dataset(cal_nc)

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
