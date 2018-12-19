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

#%% open the config file
# root = tk.Tk()
# root.withdraw() # to close this stupid root window that doesn't let itself be closed anymore AT ALL otherwise
# filename_configfile = filedialog.askopenfilename()
filename_configfile = '/Users/karllapo/Desktop/proj/DarkMix_LEOPOLD/example_data/LEOPOLD_initial_config.yml'

with open(filename_configfile, 'r') as stream:
    config_user = yaml.load(stream)


#%% create directories

dir_original = os.path.join(config_user['directories']['dir_pre'], config_user['directories']['folder_raw_data'])
if not os.path.exists(dir_original): # give an error if raw data folder can't be found
    print('Raw data folder ' + dir_original + ' does not exist.')

dir_graphics = os.path.join(config_user['directories']['dir_pre'], config_user['directories']['folder_graphics'], config_user['directories']['folder_raw_data'])
if not os.path.exists(dir_graphics): # create the folder for the graphics if it doesn't already exist
    os.makedirs(dir_graphics)

dir_processed = os.path.join(config_user['directories']['dir_pre'], config_user['directories']['folder_processed'], config_user['directories']['folder_raw_data'])
if not os.path.exists(dir_processed): # create the super-folder for the processed files if it doesn't already exist
    os.makedirs(dir_processed)

# find all experiments to be processed and make a list of them
os.chdir(dir_original)
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
    if 'RBR' in c and not c[0] == '.':
        external_data_folders['RBR'] = [os.path.join(dir_ext, c, tf)
                                        for tf in temp_files
                                        if 'RBR' in tf and not c[0] == '.']
    if 'Logger' in c and not c[0] == '.':
        external_data_folders['multiplexer'] = [os.path.join(dir_ext, c, tf)
                                                for tf in temp_files
                                                if 'multiplexer_data' in tf
                                                and not c[0] == '.']

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
    internal_config[dtsf]['archive']['sourcePath'] = dir_data[dtsf]
    internal_config[dtsf]['archive']['targetPath'] = os.path.join(dir_processed, dtsf)
    internal_config[dtsf]['directories']['dirData'] = os.path.join(dir_processed, dtsf)
    internal_config[dtsf]['directories']['dirProcessed'] = os.path.join(dir_processed, dtsf)
    internal_config[dtsf]['fileName']['filePrefix'] = config_user['directories']['folder_raw_data'] + '_' + dtsf
    if not os.path.exists(internal_config[dtsf]['archive']['targetPath']): # create the folder for the processed files if it doesn't already exist
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
    external_data_flag = True
    external_data = {}
    pt100s = None
    probe1Cols = config_user['dataProperties']['probe1_value']
    probe2Cols = config_user['dataProperties']['probe2_value']
    probe1Name = config_user['dataProperties']['probe1Temperature']
    probe2Name = config_user['dataProperties']['probe2Temperature']

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
        print('-------------')
        print('Found external data stream: ' + extdat_source)
        print('-------------')
        temp_extdata = external_data_folders[extdat_source]

        # Deal with PT100s on a multiplexer
        if extdat_source == 'multiplexer':
            # Read the PT100 data from the multiplexer file
            # The columns with multiple lines creates some problems for pandas,
            # so read the column headers separately; they should be the same
            # for each multiplexer file.
            for multiplexer in temp_extdata:
                pt100s_col_names = pd.read_csv(temp_extdata[0], header=1,
                                               nrows=0, index_col=0)
                new_data = pd.read_csv(multiplexer, header=None, skiprows=4,
                                       index_col=0, parse_dates=True,
                                       infer_datetime_format=True, sep=',')
                new_data.columns = pt100s_col_names.columns
                new_data.index.names = ['time']
                new_data = new_data.resample(config_user['dataProperties']['resampling_time']).mean()
                # Read only the data; load the first one, append the others
                if pt100s is not None:
                    pt100s = pd.concat([pt100s, new_data], axis='columns')
                else:
                    pt100s = new_data

        # Deal with any RBRs
        if extdat_source == 'RBR':
            for solo in temp_extdata:
                new_data = pd.read_csv(solo, header=0,
                                       index_col=0, parse_dates=True,
                                       infer_datetime_format=True, sep=',')
                new_data.columns = [os.path.split(solo)[1].split('.')[0]]
                new_data = new_data.resample(config_user['dataProperties']['resampling_time']).mean()
                new_data.index.names = ['time']
                if pt100s is not None:
                    pt100s = pd.concat([pt100s, new_data], axis='columns')
                else:
                    pt100s = new_data

    # Now extract each individual experiment
    for dtsf in dir_data:
        print('Constructing external data stream for ' + dtsf)
        print('')

        if multiplexer:
            external_data[dtsf] = xr.Dataset()
            if hasattr(pt100s, 'experiment_flag_Max'):
                temp_pts = pt100s[(pt100s.experiment_flag_Max == -1) # This is a line specific to the DarkMix wind tunnel
                              & (pt100s.experiment_name == dtsf)]
            else: # for data without experiment_flags
                temp_pts = pt100s
            drop_columns = [cols for cols in temp_pts.columns
                            if cols not in probe1Cols
                            and cols not in probe2Cols]
            temp_pts = temp_pts.drop(drop_columns, axis=1)

            # Conver the datetime to UTC
            time = temp_pts.index.tz_localize(tz).tz_convert('UTC')

            # Create a warm and cold bath average
            probe1 = xr.DataArray.from_series(temp_pts[probe1Cols].mean(axis=1))
            probe2 = xr.DataArray.from_series(temp_pts[probe2Cols].mean(axis=1))
            external_data[dtsf] = xr.Dataset({probe1Name: ('time', probe1),
                                              probe2Name: ('time', probe2)},
                                             coords={'time': time.values})
            # external_data[dtsf] = external_data[dtsf].chunk({'time': config_user['dataProperties']['chunkSize']})
            # external_data[dtsf] = external_data[dtsf].resample(time=config_user['dataProperties']['resampling_time']).mean()

            os.chdir(internal_config[dtsf]['directories']['dirProcessed'])
            external_data[dtsf].to_netcdf('external_data.' + dtsf + '.nc')

else:
    external_data_flag = False

#%% Open the raw saved as netcdf
# Create an empty dictionary to hold the DTS data
allExperiments = {}
labels_general = config_user['loc_general']
labels_ms = config_user['loc_ms']
labels_ref_instr = config_user['loc_ref_instr']

print('-------------')
print('Calibrating and Processing DTS temperature data.')
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

    # We are going to loop over each raw dts netcdf file, performing the
    # calibration one-at-a-time.
    os.chdir(internal_config[dtsf]['directories']['dirProcessed'])
    contents = os.listdir()
    ncfiles = [file for file in contents if '.nc' in file and 'raw' in file]
    for raw_nc in ncfiles:
        dstemp = xr.open_dataset(raw_nc)

        # Assigning to the experiments dictionary
        dstemp = btmm_process.labelLoc_general(dstemp, labels_general)
        dstemp = btmm_process.labelLoc_additional(dstemp, labels_ms, 'loc_ms')
        dstemp = btmm_process.labelLoc_additional(dstemp,
                                                  labels_ref_instr,
                                                  'loc_ref_instr')

        # Drop the unnecessary negative LAF indices
        dstemp = dstemp.sel(LAF=dstemp['LAF'] > 0)

        # Resample to a common time stamp interval with the reference/bath
        # instruments. However, if we resample to 1 second for 1 minute
        # averages we will get out a matrix of mostly NaNs).
        dstemp = dstemp.resample(time=config_user['dataProperties']['resampling_time']).mean()

#%% Construct dataset with all experiments/over all the measurement duration
        if external_data_flag:
            # Let's get a deep copy of the reference data to work with.
            ext_dat = copy.deepcopy(external_data[dtsf])
            # Reindex the reference data to our chunk of DTS data.
            ext_dat = ext_dat.reindex_like(dstemp.time,
                                           method='nearest',
                                           tolerance=1)

            # Assign the pt100s to the xarray datasets
            dstemp[probe1Name] = ext_dat[probe1Name]
            dstemp[probe2Name] = ext_dat[probe2Name]

        # Calibrate the temperatures, if the bath pt100s and dts do not line up in time, do not calibrate
        if (np.size(np.flatnonzero(~np.isnan(dstemp.temp.values))) > 0) and (config_user['flags']['calibration_flag']):
            temp_array, _, _, _ = btmm_process.matrixInversion(dstemp, internal_config[dtsf])
        elif not config_user['flags']['calibration_flag']:
            # Just return some nans here, do not notify the user as they should expect this behavior.
            temp_array = xr.Dataset({'manualTemp':
                                    (['time', 'LAF'],
                                     np.ones_like(dstemp.temp.values)
                                     * np.nan)},
                                    coords={'LAF': dstemp.LAF,
                                            'time': dstemp.time})
        else:
            # Just return some nans here and notify the user.
            temp_array = xr.Dataset({'manualTemp':
                                    (['time', 'LAF'],
                                     np.ones_like(dstemp.temp.values)
                                     * np.nan)},
                                    coords={'LAF': dstemp.LAF,
                                            'time': dstemp.time})
            print('PT100 and DTS data do not line up in time for ' + dtsf)
            print('The cal_temp field will contain NaNs.')

        # Construct a new dataset with time as a timedelta object
        dstemp['cal_temp'] =temp_array.manualTemp

        # Concatenate into a single dataset
        if temp_dts:
            temp_dts = xr.concat([temp_dts, dstemp], dim='time')
        else:
            temp_dts = xr.Dataset(dstemp)

    # Force time to be monotonically increasing
    temp_dts = temp_dts.sortby('time')

#%% Data output
    # Output the calibrated and then processed data
    suffix_cal = config_user['fileName']['fileSuffix']

    if not os.path.exists(dtsf + '_processed.nc'):
        temp_dts.to_netcdf(dtsf + '_processed.nc', engine='netcdf4')
    else:
        print('A netCDF with the name ' + dtsf + '_processed.nc already exists. The file was not overwritten.')

    # Write a csv of the locations because Matlab has problems reading them from the netcdf
    header = ['section_name', 'beginning (LAF)', 'end (LAF)']

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

    print('')
