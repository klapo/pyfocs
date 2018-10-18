#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 10:37:49 2018

@author: anita
"""

import numpy as np
from datetime import timedelta
import pandas as pd
import xarray as xr

# OS interaction
import os
import yaml

# UBT's package for handling dts data
import btmm_process

# Ignore the future compatibility warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)


import tkinter as tk #for the dialog that lets you choose your config file
from tkinter import filedialog

import copy

import csv

#%% open the config file
root = tk.Tk()
root.withdraw() # to close this stupid root window that doesn't let itself be closed anymore AT ALL otherwise
filename_configfile = filedialog.askopenfilename()

with open(filename_configfile, 'r') as stream:
    config_user = yaml.load(stream)


#%% create directories

dir_original = os.path.join(config_user['directories']['dir_pre'], config_user['directories']['folder_raw_data'])
print(dir_original)
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
# The CR6 part will not be necessary for flyfox but leave it in here for later
DTS_folders = [c for c in contents if 'CR6' not in c and not c[0] == '.']

# Loop through all of the XT data directories
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


# archiving
if config_user['flags']['archiving_flag']:
    for dtsf in DTS_folders:
        print('archiving ', dtsf)
        btmm_process.archiver(internal_config[dtsf])

# write raw netCDF
if config_user['flags']['raw_read_flag']:
    for dtsf in DTS_folders:
        print('-------------')
        print(dtsf)
        print('-------------')
        print('creating raw netcdf for experiment: ', dtsf)
        btmm_process.archive_read(internal_config[dtsf])

#%% Open the raw saved as netcdf
# Create an empty dictionary to hold the DTS data

allExperiments = {}
labels = config_user['locations']

for dtsf in dir_data:
    print('-------------')
    print(dtsf)
    print('-------------')

    # Find all the netcdfs in this directory, open the raw data netcdfs, and put into the dictionary.
    os.chdir(internal_config[dtsf]['directories']['dirProcessed'])
    contents = os.listdir()
    ncfiles = [file for file in contents if '.nc' in file and 'raw' in file]

    allExperiments[dtsf] = {}
    ds = xr.open_mfdataset('*raw*.nc', concat_dim='time')
    ds = ds.sortby('time')
    # Assigning to the experiments dictionary
    ds = btmm_process.labelLocation(ds, labels)
    allExperiments[dtsf] = ds

    print('')

#%% Construct dataset with all experiments/over all the measurement duration

for dtsf in dir_data:

    # Initialize the tunnel_exp variable for each overarching experiment
    tunnel_exp = None

    # Provide a single message that no calibration will be performed.
    if not config_user['flags']['calibration_flag']:
        print('No calibration performed for: ' + dtsf)

    # Extract the experiment from the containing dictionary
    dstemp = allExperiments[dtsf]

    # Resample to a common time stamp interval with the reference/bath instruments
    # We need to make the resampling interval a variable from the configuration file (e.g., if we resample to 1 second for
    # 1 minute averages we will get out a matrix of mostly NaNs)
    dstemp = dstemp.resample(time=config_user['dataProperties']['resampling_time']).mean()

    # Calibrate the temperatures, if the bath pt100s and dts do not line up in time, do not calibrate
    if (np.size(dstemp.temp.where(~np.isnan(dstemp.temp), drop=True)) > 0) and (config_user['flags']['calibration_flag']):
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

    # Now construct a new array to store this data

    # New coordinate for time that is time since beginning of experiment.
    LAF = dstemp.LAF

    # Construct a new dataset with time as a timedelta object
    dstemp = xr.Dataset({'instr_temp': (['time', 'LAF'], dstemp.temp),
                         'cal_temp': (['time', 'LAF'], temp_array.manualTemp),
                         'warmProbe': (['time'], dstemp['warmProbe']),
                         'coldProbe': (['time'], dstemp['coldProbe']),
                         # 'location': (['LAF'], dstemp.location),
                         },
                        coords={'time': dstemp.time,
                                'LAF': LAF,
                                'location': (['LAF'], dstemp['location'])})

    # Drop the unnecessary negative LAF indices
    dstemp = dstemp.sel(LAF=dstemp['LAF'] > 0)
    # dstemp['location'] =  [0 if v is None else v for v in dstemp.location]
    # Concatenate into a single dataset
    if tunnel_exp:
        tunnel_exp = xr.concat([tunnel_exp, dstemp], coords='all')
    else:
        tunnel_exp = xr.Dataset(dstemp)

#%% Data output
    # Output the calibrated and then processed data

    os.chdir(internal_config[dtsf]['directories']['dirProcessed'])
    if not os.path.exists(dtsf + '_calibrated.nc'):
        # Save to a netcdf
        tunnel_exp.to_netcdf(dtsf + '_calibrated.nc', engine="netcdf4")
    else:
        print('A netCDF with the name ' + dtsf + '_calibrated.nc already exists. The file was not overwritten.')

    if not os.path.exists(dtsf + '_processed.nc'):
        ds.to_netcdf(dtsf + '_processed.nc', engine='netcdf4')
    else:
        print('A netCDF with the name ' + dtsf + '_processed.nc already exists. The file was not overwritten.')

    # Write a csv of the locations because Matlab has problems reading them from the netcdf
    header = ['section_name', 'beginning (LAF)', 'end (LAF)']

    filename_locations = os.path.join(dir_processed, dtsf, dtsf + '_locations.csv')
    locations_file = open(filename_locations, 'w')
    with locations_file:
        writer = csv.writer(locations_file, delimiter=',')
        writer.writerow(header)
        for key in config_user['locations']:
            writer.writerow([key, np.min(config_user['locations'][key]), np.max(config_user['locations'][key])])

#%%######
# # Plotting!
# import FlyFox_plotting_v1
#
# if config_user['flags']['plotting_flag']:
#     for dtsf in dir_data:
#         print('-------------')
#         print('Plotting experiment: ' + dtsf)
#         print('-------------')
#
#         # Find all the netcdfs in this directory, open the raw data netcdfs, and put into the dictionary.
#         os.chdir(internal_config[dtsf]['directories']['dirProcessed'])
#         ds = xr.open_dataset(dtsf + '_calibrated.nc')
#         FlyFox_plotting_v1.flyfox_heatmap(ds)
#
#         fig, ax = plt.subplots(1, 1, figsize=(12, 8))
#
#         ax.pcolor(ds.time, ds.LAF, ds.cal_temp, cmap='viridis')
#         ax.set_xlabel('time')
#         ax.set_ylabel('LAF (m)')
#         fig.colorbar()
#
#         # Save the figure
#         os.chdir(dir_graphics)
#         fig.savefig('caltemp_heatmap_' + dtsf + config_user['fileName']['fileSuffix'] + '.pdf')
#         fig.savefig('caltemp_heatmap_' + dtsf + config_user['fileName']['fileSuffix'] + '.png')
