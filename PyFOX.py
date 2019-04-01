#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 10:37:49 2018
@author: anita
"""

# from IPython import get_ipython #clearing variables that python still knows from before
# get_ipython().magic('reset -sf')

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

# UBT's package for handling dts data
import btmm_process

# Ignore the future compatibility warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)

filename_configfile = '/Users/karllapo/Desktop/proj/DarkMix_SOPHABS/data/SOPHABS_config.yml'

with open(filename_configfile, 'r') as stream:
    config_user = yaml.load(stream)

#%% create directories

dir_original = os.path.join(config_user['directories']['dir_pre'],
                            config_user['directories']['folder_raw_data'])
# give an error if raw data folder can't be found
if not os.path.exists(dir_original):
    print('Raw data folder ' + dir_original + ' does not exist.')

dir_graphics = os.path.join(config_user['directories']['dir_pre'],
                            config_user['directories']['folder_graphics'],
                            config_user['directories']['folder_raw_data'])
# create the folder for the graphics if it doesn't already exist
# if not os.path.exists(dir_graphics):
#     os.makedirs(dir_graphics)

dir_processed = os.path.join(config_user['directories']['dir_pre'],
                             config_user['directories']['folder_processed'],
                             config_user['directories']['folder_raw_data'])
# create the super-folder for the processed files if it doesn't already exist
if not os.path.exists(dir_processed):
    os.makedirs(dir_processed)

# find all experiments to be processed/calibrated and make a list of them
try:
    os.chdir(dir_original)
except FileNotFoundError:
    os.chdir(dir_processed)
contents = os.listdir()
DTS_folders = [c for c in contents if 'external' not in c and not c[0] == '.']

# Now find the external data streams (if any)
dir_ext = os.path.join(config_user['directories']['dir_pre'], 'external')
external_data_content = os.listdir(dir_ext)
external_data_folders = {}
for c in external_data_content:
    if os.path.isdir(os.path.join(dir_ext, c)):
        temp_files = os.listdir(os.path.join(dir_ext, c))
    else:
        continue
    # Look for key phrases used in indicating an external data file.
    # Check for the RBR solo PT100s
    if 'RBR' in c.casefold() and not c[0] == '.':
        external_data_folders['RBR'] = [os.path.join(dir_ext, c, tf)
                                        for tf in temp_files
                                        if 'RBR' in tf and not c[0] == '.']
    # Check for a datalogger or multiplexer TOA5 file.
    if ('logger' in c.casefold() or 'ts_data' in c.casefold()) and not c[0] == '.':
        external_data_folders['logger'] = [os.path.join(dir_ext, c, tf)
                                           for tf in temp_files
                                           if '.dat' in tf]

# Loop through all of the DTS data directories
# assemble internal config for each dts folder within the experiment folder
dir_data = {}
internal_config = {}
for dtsf in DTS_folders:
    # get all the directories within this experiment folder
    dir_data[dtsf] = os.path.join(dir_original, dtsf)

for dtsf in DTS_folders:
    internal_config[dtsf] = {}
    internal_config[dtsf] = copy.deepcopy(config_user)
    internal_config[dtsf]['archive']['channelName'] = config_user['directories']['channelName']
    internal_config[dtsf]['archive']['sourcePath'] = dir_data[dtsf]
    internal_config[dtsf]['archive']['targetPath'] = os.path.join(dir_processed, dtsf)
    internal_config[dtsf]['directories']['dirData'] = os.path.join(dir_processed, dtsf)
    internal_config[dtsf]['directories']['dirProcessed'] = os.path.join(dir_processed, dtsf)
    internal_config[dtsf]['directories']['fileName']['prefix'] = config_user['directories']['folder_raw_data'] + '_' + dtsf
    # create the folder for the processed files if it doesn't already exist
    if not os.path.exists(internal_config[dtsf]['archive']['targetPath']):
        os.makedirs(internal_config[dtsf]['archive']['targetPath'])

# Archive/read the raw xml files for each labeled experiment.
for dtsf in DTS_folders:
    print('-------------')
    print(dtsf)
    print('-------------')
    # archiving
    if config_user['flags']['archiving_flag']:
            print('archiving ', dtsf)
            btmm_process.archiver(internal_config[dtsf])

    # write raw netCDF
    if config_user['flags']['raw_read_flag']:
            print('creating raw netcdf for experiment: ', dtsf)
            btmm_process.archive_read(internal_config[dtsf])

#%% Deal with an external datastream
# Go to the data logger directory and find the multiplexer files
if config_user['flags']['ref_temp_flag'] == 'external':
    print('-------------')
    print('Working on external data streams....')
    print('-------------')
    external_data_flag = True
    external_data = {}
    pt100s = None
    probe1Cols = config_user['dataProperties']['probe1_value']
    probe2Cols = config_user['dataProperties']['probe2_value']
    probe1Name = config_user['dataProperties']['probe1Temperature']
    probe2Name = config_user['dataProperties']['probe2Temperature']
    other_cols = config_user['dataProperties']['other_cols']

    # Construct the time zone argument
    offset = config_user['dataProperties']['UTC_offset']
    sign = offset * -1
    # Flip the sign for the POSIX convention.
    if sign > 0:
        sign = '+'
    elif sign < 0:
        sign = '-'
    tz = 'Etc/GMT' + sign + str(offset)

    for extdat_source in external_data_folders:
        # The external data sources are grouped according to type:
        # self-contained RBRs and data loggers. These data streams are dealt
        # with separately (indicated by extdat_source)
        print('Found external data stream: ' + extdat_source)
        temp_extdata = external_data_folders[extdat_source]

        # Read the PT100 data from the TOA5/csv file. These files have 4
        # header lines and quotation marks on every line.
        # Columns spread over multiple lines creates problems for pandas,
        # so read the column headers separately; they should be the same
        # for each TOA5 file.
        logger_col_names = pd.read_csv(temp_extdata[0], header=1,
                                       nrows=0, index_col=0)
        # Empty dictionary for the pandas Dataframe objects.
        pt100s = []

        # Loop over every TOA5 file found, load, format, resample, and concat
        for dat in temp_extdata:
            # Deal with PT100s on a multiplexer
            if extdat_source == 'logger':
                new_data = pd.read_csv(dat, header=None, skiprows=4,
                                       index_col=0, parse_dates=True,
                                       infer_datetime_format=True, sep=',')
                new_data.columns = logger_col_names.columns
            # Deal with any RBRs
            elif extdat_source == 'RBR':
                new_data = pd.read_csv(dat, header=0,
                                       index_col=0, parse_dates=True,
                                       infer_datetime_format=True, sep=',')
                new_data.columns = [os.path.split(dat)[1].split('.')[0]]

#             # Resample to the desired frequency.
            new_data.index.names = ['time']
            # Add the external data stream to a dictionary.
            pt100s.append(new_data)
            # Concatenate all of the external data files into a DataFrame

        if extdat_source == 'RBR':
            # Put together the RBR's into a single Dataframe
            pt100s = pd.concat(pt100s, axis='columns')
        elif extdat_source == 'logger':
            # Concatenate along the time dimension for the logger data.
            pt100s = pd.concat(pt100s, axis=0)

        # Need to force a monotically increasing time series here.

        # This check is specific to the DarkMix wind tunnel experiments.
        if hasattr(pt100s, 'experiment_flag_Max'):
            # Use this flag to determine if we have one-off experiments.
            experiment_flag = True
            for dtsf in dir_data:
                print('Constructing external data stream for ' + dtsf)
                print('')

                external_data[dtsf] = xr.Dataset()

                # Match DTS configure names to data logger experiment names.
                temp_pts = pt100s[(pt100s.experiment_flag_Max == -1)
                              & (pt100s.experiment_name == dtsf)]

                # Drop unneeded columns.
                drop_columns = [cols for cols in temp_pts.columns
                                if cols not in probe1Cols
                                and cols not in probe2Cols
                                and cols not in other_cols]
                temp_pts = temp_pts.drop(drop_columns, axis=1)

                # Conver the datetime to UTC
                time = temp_pts.index.tz_localize(tz).tz_convert('UTC')

                # Create a warm and cold bath average
                probe1 = xr.DataArray.from_series(temp_pts[probe1Cols].mean(axis=1))
                probe2 = xr.DataArray.from_series(temp_pts[probe2Cols].mean(axis=1))
                external_data[dtsf] = xr.Dataset({probe1Name: ('time', probe1),
                                                  probe2Name: ('time', probe2)},
                                                 coords={'time': time.values})

                os.chdir(internal_config[dtsf]['directories']['dirProcessed'])
                external_data[dtsf].to_netcdf('external_data.' + dtsf + '.nc')

        # External data streams when we do not have experiment_flags
        # (most other observations).
        else:
            experiment_flag = False
            # Conver the datetime to UTC
            time = pt100s.index.tz_localize(tz).tz_convert('UTC')

            # Create a warm and cold bath average
            if ((probe1Cols and probe1Cols in pt100s) and (probe2Cols and probe2Cols in pt100s)):
                probe1 = xr.DataArray.from_series(pt100s[probe1Cols].mean(axis=1))
                probe2 = xr.DataArray.from_series(pt100s[probe2Cols].mean(axis=1))
                external_data = xr.Dataset({probe1Name: ('time', probe1),
                                            probe2Name: ('time', probe2)},
                                           coords={'time': time.values})
            else:
                # No reference PT100 columns were found
                external_data = None

            if any(set(other_cols).intersection(set(pt100s))):
                # Search for the other data columns specified by the user
                others = []
                for oc in other_cols:
                    others.append(xr.DataArray.from_series(pt100s[oc]))
                others = xr.merge(others)

                # Merge the reference and other datastreams
                if external_data:
                    external_data = xr.concat([others,
                                               external_data],
                                              dim='time')
                else:
                    external_data = others

            # Sort by date to create a monotonic Dataset
            external_data = external_data.sortby(external_data.time)
            # Create a regular interval time stamp index
            reg_time_pd = pd.date_range(start=pd.to_datetime(external_data.time.values[0]),
                                        end=pd.to_datetime(external_data.time.values[0]),
                                        freq=delta_t)
            external_data = external_data.interp(time=reg_time_pd,
                                                 kwargs={'fill_value': 'extrapolate'})
            external_data = external_data.resample(config_user['dataProperties']['resampling_time']).mean()

else:
    external_data_flag = False

#%% Open the raw saved as netcdf
# Create an empty dictionary to hold the DTS data
allExperiments = {}
labels_general = config_user['loc_general']
if 'loc_ms' in config_user:
    labels_ms = config_user['loc_ms']
else:
    labels_ms = None
if 'loc_ref_instr' in config_user:
    labels_ref_instr = config_user['loc_ref_instr']
else:
    labels_ref_instr = None

print('-------------')
print('Processing and Calibrating DTS temperature data.')
print('-------------')
for dtsf in dir_data:
    # Find all the netcdfs for this experiment, open the raw data netcdfs,
    # append reference data (if applicable), calibrate (if applicable),
    # Save the output one-at-a-time per experiment label.

    # Provide a single message whether calibration will be performed.
    if not config_user['flags']['calibration_flag']:
        print('No calibration performed for: ' + dtsf)
    else:
        print('Calibrating: ' + dtsf)

    # Initialize the dts Dataset
    temp_dts = None

    # Find all 'raw' netcdfs within the processed directory,
    # sort them (by date), and process each individually.
    os.chdir(internal_config[dtsf]['directories']['dirProcessed'])
    contents = os.listdir()
    ncfiles = [file for file in contents if '.nc' in file and 'raw' in file]
    ncfiles.sort()
    ntot = np.size(ncfiles)

    for nraw, raw_nc in enumerate(ncfiles):
        # Name of the output file for this archive chunk
        outname_suffix = config_user['directories']['fileName']['suffix']
        outname_date = '_'.join(raw_nc.split('.')[0].split('_')[1:])
        outname = dtsf + '_processed_' + outname_date + outname_suffix + '.nc'
        print('Processing ' + outname_date + ' (' + str(nraw + 1)
              + ' of ' + str(ntot) + ')')

        # Open up raw file for this interval.
        dstemp = xr.open_dataset(raw_nc)

        # Assigning to the experiments dictionary
        dstemp = btmm_process.labelLoc_general(dstemp, labels_general)
        if labels_ms:
            dstemp = btmm_process.labelLoc_additional(dstemp,
                                                      labels_ms,
                                                      'loc_ms')
        if labels_ref_instr:
            dstemp = btmm_process.labelLoc_additional(dstemp,
                                                      labels_ref_instr,
                                                      'loc_ref_instr')

        # Drop the unnecessary negative LAF indices
        dstemp = dstemp.sel(LAF=dstemp['LAF'] > 0)
        dstemp.attrs['LAF_beg'] = dstemp.LAF.values[0]

        # Resample to a common time stamp interval with the reference/bath
        # instruments. Do this using a linear interpolation.

        # Round to the nearest time interval
        delta_t = config_user['dataProperties']['resampling_time']
        dt_start = pd.Timestamp(dstemp.time.values[0]).round(delta_t)
        dt_end = pd.Timestamp(dstemp.time.values[-1]).round(delta_t)

        # Create a regular interval time stamp index
        reg_time_pd = pd.date_range(start=dt_start, end=dt_end, freq=delta_t)
        dstemp = dstemp.interp(time=reg_time_pd,
                               kwargs={'fill_value': 'extrapolate'})

        # Add a delta t attribute
        dstemp.attrs['dt'] = config_user['dataProperties']['resampling_time']

        # Step loss corrections if they are provided.
        if 'step_loss_LAF' in config_user and 'step_loss_correction' in config_user:
            splice_LAF = np.atleast_1d(config_user['step_loss_LAF'])
            step_loss_corrections = np.atleast_1d(config_user['step_loss_correction'])

            # Make sure these are not NaNs
            if splice_LAF and step_loss_corrections:
                # Calculate the log power of the stokes/anti-stokes scattering
                dstemp['logPsPas'] = np.log(dstemp.Ps / dstemp.Pas)
                # Correct for any step-losses due to splicing.
                for spl_num, spl_LAF in enumerate(splice_LAF):
                    dstemp['logPsPas'] = dstemp.logPsPas.where((dstemp.LAF < spl_LAF),
                                                               dstemp.logPsPas + step_loss_corrections[spl_num])

#%% Construct dataset with all experiments/over all the measurement duration
        if external_data_flag:
            if experiment_flag:
                # Let's get a deep copy of the reference data to work with.
                ext_dat = copy.deepcopy(external_data[dtsf])
            else:
                ext_dat = copy.deepcopy(external_data)
            # Reindex the reference data to our chunk of DTS data.
            ext_dat = ext_dat.reindex_like(dstemp.time,
                                           method='nearest',
                                           tolerance=1)

            # Assign the external data to the DTS Dataset
            for ed in external_data:
                dstemp[ed] = ext_dat[ed]

       # Calibrate the temperatures. If the bath pt100s and dts do not line up
        # in time, do not calibrate.
        if (np.size(np.flatnonzero(~np.isnan(dstemp.temp.values))) > 0) and (config_user['flags']['calibration_flag']):
            dstemp, _, _, _ = btmm_process.matrixInversion(dstemp, internal_config[dtsf])
        elif not config_user['flags']['calibration_flag']:
            # Just remind the user.
            print('No calibration required, the cal_temp field will not be returned.')
        else:
            # Notify the user.
            print('PT100 and DTS data do not line up in time for ' + dtsf)
            print('The cal_temp field will not be returned.')

        # Rename the instrument reported temperature field
        dstemp.rename({'temp': 'instr_temp'}, inplace=True)

        # Output the calibrated data
        # if not os.path.exists(outname):
        dstemp.to_netcdf(outname, engine='netcdf4')
        print(dstemp)
        # else:
            # print('A netCDF with the name ' + dtsf + '_processed.nc already exists. The file was not overwritten.')

        print('')

# Write a csv of the locations because Matlab has problems reading them from the netcdf
header = ['section_name', 'beginning (LAF)', 'end (LAF)']

if (('loc_general' in config_user) and
    ('loc_ms' in config_user) and
    ('loc_ref_instr') and
    (config_user['loc_general'] is not None) and
    (config_user['loc_ms'] is not None) and
    (config_user['loc_ref_instr'] is not None)):
    filename_locations = os.path.join(dir_processed, dtsf, dtsf + '_locations.csv')
    locations_file = open(filename_locations, 'w')
    with locations_file:
        writer = csv.writer(locations_file, delimiter=',')
        writer.writerow(header)
        for key in config_user['loc_general']:
            writer.writerow([key, np.min(config_user['loc_general'][key]), np.max(config_user['loc_general'][key])])
        for key in config_user['loc_ms']:
            writer.writerow([key, np.min(config_user['loc_ms'][key]), np.max(config_user['loc_ms'][key])])
        for key in config_user['loc_ref_instr']:
            writer.writerow([key, config_user['loc_ref_instr'][key], config_user['loc_ref_instr'][key]])
