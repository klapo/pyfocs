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
# import matplotlib.pyplot as plt

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
archiving_flag = internal_config['flags']['archiving_flag']
archive_read_flag = internal_config['flags']['archive_read_flag']
calibrate_flag = internal_config['flags']['calibrate_flag']
final_flag = internal_config['flags']['final_flag']

# Step loss corrections
step_loss_flag = internal_config['step_loss']['flag']
if step_loss_flag:
    splice_LAF = internal_config['step_loss']['LAF']
    step_loss_corr = internal_config['step_loss']['correction']

# Channels/experiments/output file names
channelNames = internal_config['channelNames']
experiment_names = internal_config['experiment_names']
outname_suffix = internal_config['outname_suffix']

# -----------------------------------------------------------------------------
# Archive and create raw netcdfs
# -----------------------------------------------------------------------------
# %% Archive/read the raw xml files
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

        # Get the external data if it was indicated.
        if cal['external_flag']:
            # Get the metdata
            ref_data = xr.open_dataset(internal_config['external_data'])
            # ref_data = ref_data.resample(time=delta_t).interpolate('linear')
        elif not cal['external_flag']:
            ref_data = None

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
            if cal['double_ended']:
                outname_channel = 'merged'
                fw_channel = raw_nc.split('.')[0].split('_')[-2]
            else:
                outname_channel = raw_nc.split('.')[0].split('_')[-2]
            outname_date = raw_nc.split('.')[0].split('_')[-1]

            # The file name for the calibrated data
            outname = '_'.join(filter(None, [exp_name, 'cal',
                                             outname_channel,
                                             outname_date,
                                             outname_suffix,
                                             ])) + '.nc'
            # Check of the files to create and determine if we are overwriting.
            if write_mode == 'preserve':
                dir_cal = internal_config[exp_name]['directories']['dirCalibrated']

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
            print('Calibrating ' + outname_date + ' (' + str(nraw + 1)
                  + ' of ' + str(ntot) + ')')

            # Open up raw file for this interval.
            os.chdir(internal_config[exp_name]['directories']['dirRawNetcdf'])
            dstemp = xr.open_dataset(raw_nc)

            # Assign reference data
            dstemp, probe_names = pyfocs.assign_ref_data(dstemp,
                                                         cal,
                                                         ref_data=ref_data)

            # Drop the unnecessary negative LAF indices
            LAFmin = internal_config['min_fiber_limit']
            LAFmax = internal_config['max_fiber_limit']
            if LAFmax == -1:
                LAFmax = dstemp.LAF.values[-1]

            # Execute single-ended methods
            if not cal['double_ended']:
                # Apply fiber limits
                dstemp = dstemp.sel(LAF=slice(LAFmin, LAFmax))
                dstemp.attrs['LAF_beg'] = dstemp.LAF.values[0]
                dstemp.attrs['LAF_end'] = dstemp.LAF.values[-1]
                dstemp.attrs['calibration_method'] = cal['method']

                # Step loss corrections if they are provided.
                # This step should only be executed for single-ended methods.
                if step_loss_flag:
                    # Calculate the log power of the stokes/anti-stokes scattering
                    dstemp['logPsPas'] = np.log(dstemp.Ps / dstemp.Pas)
                    # Correct for any step-losses due to splicing.
                    for spl_num, spl_LAF in enumerate(splice_LAF):
                        dstemp['logPsPas'] = dstemp.logPsPas.where((dstemp.LAF < spl_LAF),
                                                                   dstemp.logPsPas + step_loss_corr[spl_num])

                # Calibrate the temperatures using the
                # explicit matrix inversion.
                if cal['method'] == 'matrix':
                    # Label the calibration locations
                    dstemp = pyfocs.labelLoc_additional(dstemp,
                                                        cal['library'],
                                                        'calibration')
                    cal_baths = [c for c in cal['library']
                                 if cal['library'][c]['type'] == 'calibration']
                    # Build a temporary configuration dictionary
                    # for this method.
                    temp_cfg = {}
                    for ncb, cb in enumerate(cal_baths):
                        s_ncb = str(ncb + 1)
                        temp_cfg['refField' +
                            s_ncb] = cal['library'][cb]['ref_sensor']
                        temp_cfg['refLoc' + s_ncb] = cb

                    dstemp, _, _, _ = pyfocs.matrixInversion(dstemp, temp_cfg)
                    # Drop the instrument reported temperature
                    dstemp = dstemp.drop('temp')

                if (cal['method'] == 'ols single'
                        or cal['method'] == 'wls single'):
                    dstore = pyfocs.data.to_datastore(dstemp,
                                                      internal_config,
                                                      False)
                    method = cal['method'].split()[0]
                    dstore = pyfocs.single_calibrate(dstore, method=method)
                    dstemp = pyfocs.data.from_datastore(
                        dstore,
                        datavars=probe_names
                    )

            elif cal['double_ended']:
                method = cal['method'].split()[0]
                dt = internal_config['resampling_time']

                # Get the backwards direction.
                bw_channel = cal['bw_channel']
                bw_raw_nc = raw_nc.replace(fw_channel,
                                           bw_channel)
                bw_dstemp = xr.open_dataset(bw_raw_nc)
                bw_dstemp, _ = pyfocs.assign_ref_data(bw_dstemp,
                                                      cal,
                                                      ref_data=ref_data)

                # Load both into memory to avoid annoying overhead.
                dstemp.load()
                bw_dstemp.load()

                # Line up the time stamps approximately.
                # This step means that the forward channel is
                # defined as the channel that is sampled first.
                bw_dstemp['time'] = bw_dstemp['time'] - pd.Timedelta(dt)

                len_fw = len(dstemp.time)
                len_bw = len(bw_dstemp.time)
                if not len_fw == len_bw:
                    len_time_off = len_fw - len_bw

                    if abs(len_time_off) > 2:
                        # @ Fix this error message to be less bad.
                        mess = ('The time dimension between the forward and '
                                'channels is off by more than a small amount.')
                        raise ValueError(mess)

                    # The logic here is that the forward channel is defined
                    # by being the channel sampled first. Therefor if the
                    # forward channel has 1 more time step we go searching in
                    # the next interval for the missing backwards measurement.
                    # If the backwards channels has one more measurement, we
                    # reindex to the forward channel to discard the extra
                    # observation that belongs to the previous interval.

                    # The backwards channel has more time steps or we can't
                    # grab more data.
                    if len_time_off < 0 or (nraw + 1) == ntot:
                        bw_dstemp = bw_dstemp.reindex_like(
                            dstemp,
                            method='nearest')

                    # The forward channel has more time steps.
                    elif len_time_off > 0 and not (nraw + 1) == ntot:
                        file_list = ncfiles[nraw:nraw + 1]
                        file_list = [fl.replace(fw_channel, bw_channel)
                                     for fl in file_list]
                        bw_dstemp = xr.open_mfdataset(
                            file_list,
                            combine='by_coords')
                        bw_dstemp.load()
                        bw_dstemp['time'] = bw_dstemp['time'] - \
                            pd.Timedelta(dt)
                        bw_dstemp = bw_dstemp.reindex_like(
                            dstemp,
                            method='nearest')
                        bw_dstemp, _ = pyfocs.assign_ref_data(bw_dstemp,
                                                              cal,
                                                              ref_data=ref_data)

                bw_dstemp = pyfocs.to_datastore(
                    bw_dstemp,
                    internal_config,
                    True)
                dstemp = pyfocs.to_datastore(
                    dstemp,
                    internal_config,
                    True)
                dstemp = pyfocs.merge_single(
                    dstemp,
                    bw_dstemp,
                    fixed_shift=cal['fixed_shift'])

                # Apply fiber limits after merging.
                dstemp = dstemp.sel(x=slice(LAFmin, LAFmax))
                dstemp.attrs['LAF_beg'] = dstemp.x.values[0]
                dstemp.attrs['LAF_end'] = dstemp.x.values[-1]
                dstemp.attrs['calibration_method'] = cal['method']

                # If the LAF is to be adjusted we must address the sections.
                sections = copy.deepcopy(dstemp.sections)
                for probe in sections:
                    delete_index = []
                    for nsect, probe_s in enumerate(sections[probe]):
                        if probe_s.start > LAFmax or probe_s.stop > LAFmax:
                            delete_index.append(nsect)
                        if probe_s.start < LAFmin or probe_s.stop < LAFmin:
                            delete_index.append(nsect)

                    # Delete from the last index to the beginning so that the index always works.
                    delete_index.sort()
                    delete_index.reverse()
                    for d in delete_index:
                        del sections[probe][d]

                # Rebuild the sections
                dstemp.sections = sections

                # Finally, we can calibrate.
                dstemp = pyfocs.double_calibrate(dstemp, method=method)
                dstemp = pyfocs.from_datastore(
                    dstemp,
                    datavars=probe_names,
                )

            dstemp.attrs['calibration_method'] = cal['method']

            # Assign the calibration and any available locations
            dstemp = pyfocs.labelLoc_additional(dstemp,
                                                cal['library'],
                                                'calibration')
            if 'phys_locs' in internal_config:
                for ploc in internal_config['phys_locs']:
                    # Label all physical locations provided
                    dstemp = pyfocs.labelLoc_additional(dstemp,
                                                        lib[ploc],
                                                        ploc)

            # @ Store the calibration decisions here.
            os.chdir(internal_config[exp_name]['directories']['dirCalibrated'])
            dstemp.to_netcdf(outname, engine='netcdf4')
            dstemp.close()

            print('')

# -----------------------------------------------------------------------------
# Finalize with physical coordinates
# -----------------------------------------------------------------------------

if final_flag:
    print('-------------')
    print('Final preparation: assigning physical coordinates, dropping unused fields and locations')
    print(' ')

    phys_locs = internal_config['phys_locs']
    align_locations = internal_config['align_locations']

    # Time step for resampling to a uniform time step
    delta_t = internal_config['resampling_time']

    # Time limits for processing only a subsection
    if internal_config['time_limits']['flag']:
        tstart = pd.Timestamp(internal_config['tstart'])
        tstop = pd.Timestamp(internal_config['tstop'])

    # When finalizing the dataset all extraneous coordinates and data
    # is dropped, leaving behind these variables.
    coords_to_keep = ['xyz', 'time', 'x', 'y', 'z', 'LAF']
    vars_to_keep = ['cal_temp']

    finished_files = []

    for exp_name in experiment_names:
        # Find all 'calibrated' netcdfs within the calibrated directory,
        # sort them by date, suffix, and process each individually.

        os.chdir(internal_config[exp_name]['directories']['dirCalibrated'])
        contents = os.listdir()
        if outname_suffix:
            ncfiles = [file for file in contents
                if '.nc' in file
                and 'cal' in file
                and outname_suffix in file]
        else:
            ncfiles = [file for file in contents
                if '.nc' in file
                and 'cal' in file]
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
            if (internal_config['time_limits']['flag'] and
                    (pd.Timestamp(dstemp.time.min().values) < tstart
                    or pd.Timestamp(dstemp.time.max().values) > tstop)):
                print('Outside time bounds. skipping...')
                dstemp.close()
                continue

            # Resample to a common time stamp interval with the reference/bath
            # instruments. Do this using a linear interpolation.
            # Round to the nearest time interval
            dt_start=pd.Timestamp(dstemp.time.values[0]).round(delta_t)
            dt_end=pd.Timestamp(dstemp.time.values[-1]).round(delta_t)

            # Catch a weird edge case for a single time step.
            if not dt_start == dt_end:
                # Create a regular interval time stamp index
                reg_time_pd=pd.date_range(start=dt_start,
                                            end=dt_end,
                                            freq=delta_t)
                dstemp=dstemp.interp(time=reg_time_pd,
                                       kwargs={'fill_value': 'extrapolate'})

            # Add a delta t attribute
            dstemp.attrs['dt']=delta_t

            # Clean up unused variables and labels.
            vars_to_drop=[v for v in dstemp.data_vars
                            if v not in vars_to_keep]
            coords_to_drop=[c for c in dstemp.coords
                              if c not in coords_to_keep]

            # Clean up attributes and dropped the unused ones.
            dt=dstemp.attrs['dt']
            dLAF=dstemp.attrs['dLAF']
            dstemp=dstemp.drop(vars_to_drop).drop(coords_to_drop)
            dstemp.attrs=[]
            dstemp.attrs['dt']=dt
            dstemp.attrs['dLAF']=dLAF

            # Use the automatic alignment between sections to map one
            # location onto another.
            if align_locations:
                common_sections=internal_config['common_sections']
                unique_sections=internal_config['unique_sections']
                locs_to_match=internal_config['location_matching']

                # @ Rename s1, s2 etc to `to` and `from` for consistent naming.
                for map_from, map_to in locs_to_match.items():
                    s1_list=[]
                    s2_list=[]
                    for section, shift in common_sections[map_from].items():
                        s1, s2, lib = pyfocs.interp_section(
                            dstemp, lib, map_to, map_from, section,
                            fixed_shift=shift,
                            dl=10, plot_results=True)

                        # If they are not the same size then this step failed.
                        if not s1.LAF.size == s2.LAF.size:
                            mess = (
                                '-'.join(map_from, section)
                                + ' was not successfully mapped to '
                                + '-'.join(map_to, section))
                            print(mess)
                            print('==================')
                            print(map_to)
                            print(s1)
                            print('')
                            print('==================')
                            print(map_from)
                            print(s2)
                            raise ValueError

                        s1.coords[map_to]=section
                        s2.coords[map_from]=section

                        s1_list.append(s1)
                        s2_list.append(s2)

                    ds_ploc1=xr.concat(s1_list, dim='LAF')
                    ds_ploc1=ds_ploc1.drop('x')
                    ds_ploc1=pyfocs.labeler.dtsPhysicalCoords_3d(
                        ds_ploc1,
                        lib[map_to])

                    ds_ploc2=xr.concat(s2_list, dim='LAF')
                    ds_ploc2=ds_ploc2.drop('x')
                    ds_ploc2=pyfocs.labeler.dtsPhysicalCoords_3d(
                        ds_ploc2,
                        lib[map_from])

                    # Output each location type as a separate final file.
                    os.chdir(internal_config[exp_name]
                             ['directories']['dirFinal'])
                    outname_ploc1='_'.join(filter(None, [exp_name, 'final',
                                                         outname_date,
                                                         outname_suffix,
                                                         map_to])) + '.nc'
                    outname_ploc2='_'.join(filter(None, [exp_name, 'final',
                                                         outname_date,
                                                         outname_suffix,
                                                         map_from])) + '.nc'

                    # @ Convert boolean attributes to 0/1
                    del ds_ploc1.attrs['reverse']
                    del ds_ploc2.attrs['reverse']
                    ds_ploc1.to_netcdf(outname_ploc1, mode='w')
                    ds_ploc2.to_netcdf(outname_ploc2, mode='w')

                    # @ Label unique locations

            else:
                for ploc in phys_locs:
                    # Relabel the locations. This allows locations to
                    # change after calibrating, as the calibration only
                    # cares about the location of the reference baths.
                    dstemp_ploc=pyfocs.labelLoc_additional(
                        dstemp,
                        lib[ploc],
                        ploc)

                    # Assign physical labels
                    dstemp_ploc=pyfocs.labeler.dtsPhysicalCoords_3d(dstemp_ploc,
                                                                    lib[ploc])

                    # Output each location type as a separate final file.
                    outname='_'.join(filter(None, [exp_name, 'final',
                                                   outname_date,
                                                   outname_suffix,
                                                   ploc])) + '.nc'
                    os.chdir(internal_config[exp_name]
                             ['directories']['dirFinal'])
                    dstemp_ploc.to_netcdf(outname, mode='w')

            dstemp.close()
            # Make sure we don't reprocess files.
            finished_files.extend(cal_nc)
