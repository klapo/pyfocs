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
import btmm_process

# Ignore the future compatibility warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)

# -----------------------------------------------------------------------------
# Configuration file
# -----------------------------------------------------------------------------
# Quick fix for us to do local debugging.
try:
    filename_configfile_KL = 'bla'
    filename_configfile_AF = '/home/anita/Schreibtisch/python_programs/config_files/LOVE/LOVE_Simba_south_190711_config.yml'
    if os.path.exists(filename_configfile_KL):
        filename_configfile = filename_configfile_KL
    elif os.path.exists(filename_configfile_AF):
        filename_configfile = filename_configfile_AF
    else:
        raise NameError('Config file not found.')

# Expected, non-debugging behavior.
except NameError:

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

# This yaml syntax requires pyyaml v5.1
with open(filename_configfile, 'r') as stream:
    config_user = yaml.load(stream, Loader=yaml.FullLoader)

try:
    remote_dir = config_user['directories']['remote']['remote_directories']
except KeyError:
    remote_dir = []

# -----------------------------------------------------------------------------
# Look through the configuration file and prepare
# -----------------------------------------------------------------------------
#%% loop through the experiment_names
experiment_names = config_user['directories']['experiment_names']
internal_config = {}
for exp_name in experiment_names:
    #%% Get paths for the directories of the different processing steps.
    # Create directories if they don't exist and errors if needed data aren't
    # found.

    # --------
    # External data
    try:
        dir_ext = config_user['directories']['external']
        # Handling for relative paths
        if dir_ext == '.':
            dir_ext = os.getcwd()
        # Check for empty yaml lists (assigned NoneType)
        if dir_ext is not None:
            # check if the external data directory exists, error if not
            if (not os.path.exists(dir_ext)) and (config_user['flags']['ref_temp_option']):
                raise FileNotFoundError('External data directory not found at '
                                        + dir_ext)
    except KeyError:
        dir_ext = None

    # Handling for relative paths
    try:
        dir_pre_local = config_user['directories']['local']['dir_pre']
    except KeyError:
        dir_pre_local = ''

    try:
        dir_pre_remote = config_user['directories']['remote']['dir_pre']
    except KeyError:
        dir_pre_remote = ''

    if dir_pre_local == '.':
        dir_pre_local = os.getcwd()
    if dir_pre_remote == '.':
        dir_pre_remote = os.getcwd()

    #%%
    # --------
    # Raw xml files
    if 'raw_xml' in remote_dir:
        dir_raw_xml = os.path.join(dir_pre_remote,
                                   exp_name, config_user['directories']['raw'])
    else:
        dir_raw_xml = os.path.join(dir_pre_local,
                                   exp_name, config_user['directories']['raw_xml'])

    # Throw an error if raw files should be read but directory isn't found
    if (not os.path.exists(dir_raw_xml)) and (config_user['flags']['archiving_flag']):
        raise FileNotFoundError('Raw data directory not found: ' + dir_raw_xml)

    #%%
    # --------
    # Archived tar.gz of xmls
    if 'archive' in remote_dir:
        dir_archive = os.path.join(dir_pre_remote,
                                   exp_name, config_user['directories']['archive'])
    else:
        dir_archive = os.path.join(dir_pre_local,
                                   exp_name, config_user['directories']['archive'])

    if (not os.path.exists(dir_archive)):
        # create the folder for the archived files if it doesn't already exist and archived files will be created
        if config_user['flags']['archiving_flag']:
            os.makedirs(dir_archive)
            print('Archived directory not found. Created a new one at ' + dir_archive)
        # throw error if archived files are needed but neither found nor created
        elif config_user['flags']['archive_read_flag']:
            raise FileNotFoundError('Archived data directory not found at ' +
                                    dir_archive)

    #%%
    # --------
    # Raw xmls in netcdf format
    if 'raw_netcdf' in remote_dir:
        dir_raw_netcdf = os.path.join(dir_pre_remote,
                                      exp_name, config_user['directories']['raw_netcdf'])
    else:
        dir_raw_netcdf = os.path.join(dir_pre_local,
                                      exp_name, config_user['directories']['raw_netcdf'])

    if not os.path.exists(dir_raw_netcdf):
        # Create the folder for the raw netcdf files if it doesn't already
        # exist and raw netcdf files will be created.
        if config_user['flags']['archive_read_flag']:
            os.makedirs(dir_raw_netcdf)
            print('Calibrated directory not found. Created a new one at ' + dir_raw_netcdf)
        # throw error if raw netcdf files are needed but neither found nor created
        elif config_user['flags']['calibrate_flag']:
            raise FileNotFoundError('Raw netcdf data directory not found ' +
                                    dir_raw_netcdf)

    #%%
    # --------
    # Calibrated temperature (full matrix inversion)
    if 'calibrated' in remote_dir:
        dir_cal = os.path.join(dir_pre_remote,
                               exp_name, config_user['directories']['calibrated'])
    else:
        dir_cal = os.path.join(dir_pre_local,
                               exp_name, config_user['directories']['calibrated'])

    if not os.path.exists(dir_cal):
        # create the folder for the calibrated files if it doesn't already
        # exist and calibrated files will be created
        if config_user['flags']['calibrate_flag']:
            os.makedirs(dir_cal)
            print('Calibrated directory not found. Created a new one at '
                  + dir_cal)
        elif config_user['flags']['final_flag']:
            raise FileNotFoundError('Calibrated data directory not found ' +
                                    dir_cal)

    #%%
    # --------
    # Final data (converted into just the labeled segments with
    # a physical location)
    if 'final' in remote_dir:
        dir_final = os.path.join(dir_pre_remote,
                                 exp_name, config_user['directories']['final'])
    else:
        dir_final = os.path.join(dir_pre_local,
                                 exp_name, config_user['directories']['final'])

    # Create the folder for the final files if it doesn't already exist
    if (not os.path.exists(dir_final)) and (config_user['flags']['calibrate_flag']):
        os.makedirs(dir_final)
        print('Final directory not found. Created a new one at ' + dir_final)

    #%%
    # --------
    # Graphics folder (currently unused)
    if 'graphics' in remote_dir:
        dir_graphics = os.path.join(dir_pre_remote,
                                    exp_name, config_user['directories']['graphics'])
    else:
        dir_graphics = os.path.join(dir_pre_local,
                                    exp_name, config_user['directories']['graphics'])

    # assemble internal config for each dts folder within the experiment folder
    internal_config[exp_name] = {}
    internal_config[exp_name] = copy.deepcopy(config_user)
    internal_config[exp_name]['archive']['channelName'] = config_user['directories']['channelName']
    internal_config[exp_name]['archive']['sourcePath'] = dir_raw_xml
    internal_config[exp_name]['archive']['targetPath'] = dir_archive
    internal_config[exp_name]['directories']['dirArchive'] = dir_archive
    internal_config[exp_name]['directories']['dirRawNetcdf'] = dir_raw_netcdf
    internal_config[exp_name]['directories']['dirCalibrated'] = dir_cal
    internal_config[exp_name]['directories']['dirFinal'] = dir_final

    # Make sure the prefix/suffix fields exist
    try:
        outname_prefix = config_user['directories']['prefix']
    except KeyError:
        outname_prefix = None

    try:
        outname_suffix = config_user['directories']['suffix']
    except KeyError:
        outname_suffix = None

   # Determine fiber type
    coretype = internal_config[exp_name]['dataProperties']['fiber_type']
    if coretype == 'multicore':
        cores = internal_config[exp_name]['dataProperties']['cores']
    else:
        cores = None

# Location library and data labels
# Determine the list of all location types.
try:
    loc_type = config_user['dataProperties']['all_locs']
except KeyError:
    loc_type = []
    try:
        for l in config_user['location_library']:
            loc_type.append(config_user['location_library'][l]['loc_type'])
        loc_type = np.unique(loc_type).tolist()
    except KeyError:
        print('No location library found, all steps after creating raw' +
              ' netcdfs have been disabled.')
        config_user['flags']['calibrate_flag'] = False
        config_user['flags']['processed_flag'] = False
        config_user['flags']['final_flag'] = False
    except TypeError:
        print('No location library found, all steps after creating raw' +
              ' netcdfs have been disabled.')
        config_user['flags']['calibrate_flag'] = False
        config_user['flags']['processed_flag'] = False
        config_user['flags']['final_flag'] = False

# Step loss corrections
if ('step_loss_LAF' in config_user['dataProperties']
        and 'step_loss_correction' in config_user['dataProperties']):
    step_loss_flag = True
else:
    step_loss_flag = False

# Determine write mode:
try:
    write_mode = config_user['flags']['write_mode']
except KeyError:
    write_mode = 'overwrite'

# Determine calibration mode:
try:
    cal_mode = config_user['flags']['cal_mode']
except KeyError:
    cal_mode = 'instantaneous'

# -----------------------------------------------------------------------------
# Archive and create raw netcdfs
# -----------------------------------------------------------------------------
#%% Archive/read the raw xml files
for exp_name in experiment_names:
    print('-------------')
    print(exp_name)
    print('-------------')

    # archiving
    if config_user['flags']['archiving_flag']:
        print('-------------')
        print('Archiving raw xml files.')
        print(' ')
        print('archiving ', exp_name)
        btmm_process.archiver(internal_config[exp_name])

    # write raw netCDF
    if config_user['flags']['archive_read_flag']:
        print('-------------')
        print('Writing netcdfs from raw xml files.')
        print(' ')
        print('creating raw netcdf for experiment: ', exp_name)
        btmm_process.archive_read(internal_config[exp_name],
                                  write_mode=write_mode)

for exp_name in experiment_names:
    # --------------------------------------------------------------------------
    # Process DTS files
    # --------------------------------------------------------------------------
    if config_user['flags']['calibrate_flag']:

        # Get the external data to add if we are not using the instrument PT100s
        if internal_config[exp_name]['flags']['ref_temp_option'] == 'external':
            # Get the metdata
            os.chdir(dir_ext)
            ref_data = xr.open_dataset(config_user['directories']['filename_external'])
            ref_data = ref_data.resample(time = config_user['dataProperties']['resampling_time']).interpolate('linear')

            probe1_name = internal_config[exp_name]['dataProperties']['probe1Temperature']
            probe2_name = internal_config[exp_name]['dataProperties']['probe2Temperature']

        # Find all 'raw' netcdfs within the processed directory,
        # sort them (by date), and process each individually.
        os.chdir(internal_config[exp_name]['directories']['dirRawNetcdf'])
        contents = os.listdir()
        ncfiles = [file for file in contents
                   if '.nc' in file and 'raw' in file]
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
                        outname = '_'.join(filter(None, [outname_prefix,
                                                         exp_name, 'cal',
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
                    outname = '_'.join(filter(None, [outname_prefix,
                                                     exp_name, 'cal',
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
            delta_t = config_user['dataProperties']['resampling_time']
            dt_start = pd.Timestamp(dstemp.time.values[0]).round(delta_t)
            dt_end = pd.Timestamp(dstemp.time.values[-1]).round(delta_t)

            # Create a regular interval time stamp index
            reg_time_pd = pd.date_range(start=dt_start, end=dt_end, freq=delta_t)
            dstemp = dstemp.interp(time=reg_time_pd,
                                   kwargs={'fill_value': 'extrapolate'})

            # Add a delta t attribute
            dstemp.attrs['dt'] = config_user['dataProperties']['resampling_time']

            # Add in external reference data here
            if internal_config[exp_name]['flags']['ref_temp_option'] == 'external':
                temp_ref_data = ref_data.reindex_like(dstemp.time,
                                                      method='nearest',
                                                      tolerance=1)

                dstemp[probe1_name] = temp_ref_data[probe1_name]
                dstemp[probe2_name] = temp_ref_data[probe2_name]

                # Add additional external data for this data stream.
                if internal_config[exp_name]['dataProperties']['external_fields']:
                    for ext_dat in internal_config[exp_name]['dataProperties']['external_fields']:
                        dstemp[ext_dat] = temp_ref_data[ext_dat]

                # If the bath pt100s and dts do not line up in time,
                # notify the user.
                if not (np.size(np.flatnonzero(~np.isnan(dstemp.temp.values))) > 0):
                    print('PT100 and DTS data do not line up in time for ' + raw_nc)

            # Step loss corrections if they are provided.
            if step_loss_flag:
                splice_LAF = np.atleast_1d(config_user['dataProperties']['step_loss_LAF'])
                step_loss_corrections = np.atleast_1d(config_user['dataProperties']['step_loss_correction'])

                # Make sure these are not NaNs
                if splice_LAF and step_loss_corrections:
                    # Calculate the log power of the stokes/anti-stokes scattering
                    dstemp['logPsPas'] = np.log(dstemp.Ps / dstemp.Pas)
                    # Correct for any step-losses due to splicing.
                    for spl_num, spl_LAF in enumerate(splice_LAF):
                        dstemp['logPsPas'] = dstemp.logPsPas.where((dstemp.LAF < spl_LAF),
                                                                   dstemp.logPsPas + step_loss_corrections[spl_num])

            # Split up multicore data
            if coretype == 'multicore':
                for c in cores:
                    print(c)
                    dstemp_core = dstemp

                    # Drop the unnecessary negative LAF indices
                    core_LAF_start = config_user['dataProperties']['cores'][c]['LAF'][0]
                    core_LAF_end = config_user['dataProperties']['cores'][c]['LAF'][-1]
                    dstemp_core = dstemp_core.sel(LAF=((dstemp_core['LAF'] > core_LAF_start)
                                                        & (dstemp_core['LAF'] < core_LAF_end)))
                    dstemp_core.attrs['LAF_beg'] = core_LAF_start
                    dstemp_core.attrs['LAF_end'] = core_LAF_end

                    # Location labels
                    for loc_type_cur in loc_type:
                        location = {}
                        for l in config_user['location_library']:
                            if loc_type_cur == config_user['location_library'][l]['loc_type']:
                                location[l] = config_user['location_library'][l]['LAF'][c]
                        dstemp_core = btmm_process.labelLoc_additional(dstemp_core,
                                                                       location,
                                                                       loc_type_cur)
                    # Converting to a netcdf ruins this step unfortunately.
                    # dstemp_core.attrs['loc_general_long'] = dict((l, config_user['loc_general'][l]['long name'])

                    # Calibrate the temperatures
                    if cal_mode == 'smooth':
                        dstemp_core = btmm_process.timeAvgCalibrate(dstemp_core, internal_config[exp_name])
                    else:
                        dstemp_core, _, _, _ = btmm_process.matrixInversion(dstemp_core, internal_config[exp_name])

                    # Rename the instrument reported temperature field
                    dstemp_core = dstemp_core.rename({'temp': 'instr_temp'})
                    dstemp_core.coords['core'] = c

                    # Output the calibrated dataset
                    outname_core = c
                    outname = '_'.join(filter(None, [outname_prefix,
                                                     exp_name, 'cal',
                                                     outname_channel,
                                                     outname_date,
                                                     outname_suffix,
                                                     outname_core])) + '.nc'
                    os.chdir(internal_config[exp_name]['directories']['dirCalibrated'])
                    dstemp_core.to_netcdf(outname, engine='netcdf4')

            elif coretype == 'singlecore':
                # Drop the unnecessary negative LAF indices
                dstemp = dstemp.sel(LAF=dstemp['LAF'] > 0)
                dstemp.attrs['LAF_beg'] = dstemp.LAF.values[0]

                # Location labels
                for loc_type_cur in loc_type:
                    location = {}
                    for l in config_user['location_library']:
                        if loc_type_cur == config_user['location_library'][l]['loc_type']:
                            location[l] = config_user['location_library'][l]['LAF']
                    dstemp = btmm_process.labelLoc_additional(dstemp,
                                                              location,
                                                              loc_type_cur)

                # Calibrate the temperatures.
                dstemp, _, _, _ = btmm_process.matrixInversion(dstemp, internal_config[exp_name])
                dstemp = dstemp.rename({'temp': 'instr_temp'})

                # Output the calibrated dataset
                outname_core = None
                outname = '_'.join(filter(None, [outname_prefix,
                                                 exp_name, 'cal',
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

if config_user['flags']['final_flag']:
    print('-------------')
    print('Final preparation: assigning physical coordinates, dropping unused fields and locations')
    print(' ')

    # Locations to physically label
    phys_locs = config_user['dataProperties']['phys_locs']

    # List of finished files so we can skip files when handling multicore fibers.
    finished_files = []

    # When finalizing the dataset all extraneous coordinates and data
    # is dropped, leaving behind these variables.
    coords_to_keep = ['xyz', 'time', 'x', 'y', 'z', 'core']
    vars_to_keep = ['cal_temp']

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
            try:
                outname_suffix = config_user['directories']['fileName']['suffix']
            except KeyError:
                outname_suffix = ''
            name_components = cal_nc.split('.')[0].split('_')
            outname_channel = name_components[-3]
            outname_date = name_components[-2]

            print('Finalizing ' + cal_nc + ' (' + str(ncal + 1)
                  + ' of ' + str(ntot) + ')')

            # Location dictionary to pass to the physical coordinates function.
            location = {}

            # Deal with multicore fibers differently than single core fibers.
            if coretype == 'multicore':
                os.chdir(internal_config[exp_name]['directories']['dirCalibrated'])
                # Check if we already handled this core.
                if cal_nc in finished_files:
                    continue
                nc_cal_core_name = '_'.join(name_components[:-1])
                nc_cal_core = [file for file in ncfiles if nc_cal_core_name in file and 'cal' in file]
                dstemp_out = {}

                # Assign physical coordinates. Each core is appended to a list to be merged later.
                for c in nc_cal_core:
                    dstemp = xr.open_dataset(c)
                    dstemp.load()

                    # Reformat the config locations to specify just a single core
                    core = str(dstemp.core.values)
                    # Make sure we intend to process this core
                    if core not in config_user['dataProperties']['cores']:
                        continue

                    for l in config_user['location_library']:
                        location[l] = copy.deepcopy(config_user['location_library'][l])
                        location[l]['LAF'] = config_user['location_library'][l]['LAF'][core]

                    for ploc in config_user['dataProperties']['phys_locs']:
                        if not ploc in dstemp_out:
                            dstemp_out[ploc] = []
                        # Find all the locations to label in this location type
                        temp_loc = {loc:location[loc] for loc in location if ploc==location[loc]['loc_type']}
                        dstemp_out[ploc].append(btmm_process.labeler.dtsPhysicalCoords_3d(dstemp, temp_loc))

                # Merge the cores
                for ploc in config_user['dataProperties']['phys_locs']:
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

                    # Clean up unused variables and labels.
                    vars_to_drop = [v for v in dstemp_out[ploc].data_vars
                                    if v not in vars_to_keep]

                    coords_to_drop = [c for c in dstemp_out[ploc].coords
                                      if c not in coords_to_keep]
                    coords_to_drop.remove(ploc)

                    # Clean up attributes and dropped the unused ones.
                    dt = dstemp_out[ploc].attrs['dt']
                    dLAF = dstemp_out[ploc].attrs['dLAF']
                    dstemp_out[ploc] = dstemp_out[ploc].drop(vars_to_drop).drop(coords_to_drop)
                    dstemp_out[ploc].attrs = []
                    dstemp_out[ploc].attrs['dt'] = dt
                    dstemp_out[ploc].attrs['dLAF'] = dLAF

                    # Output each location type as a separate final file.
                    outname = '_'.join(name_components[:-1]).replace('cal', 'final') + '_' + ploc + '.nc'
                    os.chdir(internal_config[exp_name]['directories']['dirFinal'])
                    dstemp_out[ploc].to_netcdf(outname)

                # Make sure we don't reprocess files.
                finished_files.extend(nc_cal_core)

            # Deal with the single core fiber.
            else:
                # Open each calibrated file.
                os.chdir(internal_config[exp_name]['directories']['dirCalibrated'])
                dstemp = xr.open_dataset(cal_nc)
                dstemp = xr.open_dataset(c)

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

                location = {}
                for l in config_user['location_library']:
                    location[l] = copy.deepcopy(config_user['location_library'][l])
                    location[l]['LAF'] = config_user['location_library'][l]['LAF']

                for ploc in config_user['dataProperties']['phys_locs']:
                    temp_loc = {loc:location[loc] for loc in location if ploc in location[loc]['loc_type']}
                    dstemp_out = btmm_process.labeler.dtsPhysicalCoords_3d(dstemp, temp_loc)
                    os.chdir(internal_config[exp_name]['directories']['dirFinal'])
                    outname = cal_nc.split('.')[0].replace('cal', 'final') + '_' + ploc + '.nc'
                    dstemp_out.to_netcdf(outname)
