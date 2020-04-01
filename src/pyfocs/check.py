import os
import yaml
from pyfocs import yamlDict
import copy
import numpy as np


def config(fn_cfg):
    '''
    Function checks the basic integrity of a configuration file and returns the
    config dictionary, location library, and physical locations dictionary
    expected by pyfocs.
    '''

    # The internal config file and location library are handed back to pyfocs
    in_cfg = {}
    lib = {}

    # This yaml syntax requires pyyaml v5.1
    try:
        cfg = yamlDict(fn_cfg)
    # The parse error appears to be raised a couple different ways. Try to
    # catch both methods of reporting a ParserError.
    except yaml.parser.ParserError as exc:
        print('The config yaml file is not properly formatted.')
        print('Examine the following bad line(s):')
        print(exc)

        return None
        # Example output:
        # while parsing a block mapping
        # in "fine_name.yml", line 2, column 1

    try:
        remote_dir = cfg['directories']['remote']['remote_directories']
    except KeyError:
        remote_dir = []

    # -------------------------------------------------------------------------
    # Look through the configuration file and prepare
    experiment_names = cfg['directories']['experiment_names']
    mess = 'experiment_names should be a list (e.g., all items start with a dash)'
    assert isinstance(cfg['directories']['experiment_names'], list), experiment_names

    # Verify that the flags exist and are booleans.
    in_cfg['flags'] = cfg['flags']
    flags = ['ref_temp_option', 'write_mode', 'archiving_flag',
             'archive_read_flag', 'calibrate_flag', 'final_flag']

    assert all([fl in flags for fl in in_cfg['flags']]), 'Not all flags were found.'
    # assert all([isinstance(fl, bool) for fl in in_cfg['flags'] if fl in bool_flags]), 'Flags must be boolean.'

    for exp_name in experiment_names:
        # Create directories if they don't exist and
        # errors if needed data are not found.

        # --------
        # External data
        try:
            dir_ext = cfg['directories']['external']
            # Handling for relative paths
            if dir_ext == '.':
                dir_ext = os.getcwd()
            # Check for empty yaml lists (assigned NoneType)
            if dir_ext is not None:
                # check if the external data directory exists, error if not
                if (not os.path.exists(dir_ext)) and (cfg['flags']['ref_temp_option']):
                    raise FileNotFoundError('External data directory not found at '
                                            + dir_ext)
        except KeyError:
            dir_ext = None

        # Handling for relative paths
        try:
            dir_pre_local = cfg['directories']['local']['dir_pre']
        except KeyError:
            dir_pre_local = ''

        try:
            dir_pre_remote = cfg['directories']['remote']['dir_pre']
        except KeyError:
            dir_pre_remote = ''

        if dir_pre_local == '.':
            dir_pre_local = os.getcwd()
        if dir_pre_remote == '.':
            dir_pre_remote = os.getcwd()

        # --------
        # Raw xml files
        if 'raw_xml' in remote_dir:
            dir_raw_xml = os.path.join(dir_pre_remote,
                                       exp_name, cfg['directories']['raw'])
        else:
            dir_raw_xml = os.path.join(dir_pre_local,
                                       exp_name, cfg['directories']['raw_xml'])

        # Throw an error if raw files should be read but directory isn't found
        if (not os.path.exists(dir_raw_xml)) and (cfg['flags']['archiving_flag']):
            raise FileNotFoundError('Raw data directory not found: ' + dir_raw_xml)

        # --------
        # Archived tar.gz of xmls
        if 'archive' in remote_dir:
            dir_archive = os.path.join(dir_pre_remote,
                                       exp_name, cfg['directories']['archive'])
        else:
            dir_archive = os.path.join(dir_pre_local,
                                       exp_name, cfg['directories']['archive'])

        if (not os.path.exists(dir_archive)):
            # create the folder for the archived files if it doesn't already
            # exist and archived files will be created
            if cfg['flags']['archiving_flag']:
                os.makedirs(dir_archive)
                print('Archived directory not found. Created a new one at ' + dir_archive)
            # throw error if archived files are needed but neither found nor created
            elif cfg['flags']['archive_read_flag']:
                raise FileNotFoundError('Archived data directory not found at ' +
                                        dir_archive)

        # --------
        # Raw xmls in netcdf format
        if 'raw_netcdf' in remote_dir:
            dir_raw_netcdf = os.path.join(dir_pre_remote,
                                          exp_name, cfg['directories']['raw_netcdf'])
        else:
            dir_raw_netcdf = os.path.join(dir_pre_local,
                                          exp_name, cfg['directories']['raw_netcdf'])

        if not os.path.exists(dir_raw_netcdf):
            # Create the folder for the raw netcdf files if it doesn't already
            # exist and raw netcdf files will be created.
            if cfg['flags']['archive_read_flag']:
                os.makedirs(dir_raw_netcdf)
                print('Calibrated directory not found. Created a new one at ' + dir_raw_netcdf)
            # throw error if raw netcdf files are needed but neither found
            # nor created
            elif cfg['flags']['calibrate_flag']:
                raise FileNotFoundError('Raw netcdf data directory not found ' +
                                        dir_raw_netcdf)

        # --------
        # Calibrated temperature
        if 'calibrated' in remote_dir:
            dir_cal = os.path.join(dir_pre_remote,
                                   exp_name, cfg['directories']['calibrated'])
        else:
            dir_cal = os.path.join(dir_pre_local,
                                   exp_name, cfg['directories']['calibrated'])

        if not os.path.exists(dir_cal):
            # create the folder for the calibrated files if it doesn't already
            # exist and calibrated files will be created
            if cfg['flags']['calibrate_flag']:
                os.makedirs(dir_cal)
                print('Calibrated directory not found. Created a new one at '
                      + dir_cal)
            elif cfg['flags']['final_flag']:
                raise FileNotFoundError('Calibrated data directory not found ' +
                                        dir_cal)

        # --------
        # Final data (converted into just the labeled segments with
        # a physical location)
        if 'final' in remote_dir:
            dir_final = os.path.join(dir_pre_remote,
                                     exp_name, cfg['directories']['final'])
        else:
            dir_final = os.path.join(dir_pre_local,
                                     exp_name, cfg['directories']['final'])

        # Create the folder for the final files if it doesn't already exist
        if ((not os.path.exists(dir_final))
                and ((cfg['flags']['calibrate_flag'])
                     or (cfg['flags']['final_flag']))):
            os.makedirs(dir_final)
            print('Final directory not found. Created a new one at ' + dir_final)

        # --------
        # Graphics folder
        if 'graphics' in remote_dir:
            dir_graphics = os.path.join(dir_pre_remote,
                                        exp_name, cfg['directories']['graphics'])
        else:
            dir_graphics = os.path.join(dir_pre_local,
                                        exp_name, cfg['directories']['graphics'])

        # assemble internal config for each dts folder within the
        # experiment folder
        in_cfg[exp_name] = {}
        in_cfg[exp_name]['archive'] = copy.deepcopy(cfg['archive'])
        in_cfg[exp_name]['archive']['channelName'] = cfg['directories']['channelName']
        in_cfg[exp_name]['archive']['sourcePath'] = dir_raw_xml
        in_cfg[exp_name]['archive']['targetPath'] = dir_archive
        in_cfg[exp_name]['directories'] = copy.deepcopy(cfg['directories'])
        in_cfg[exp_name]['directories']['dirArchive'] = dir_archive
        in_cfg[exp_name]['directories']['dirRawNetcdf'] = dir_raw_netcdf
        in_cfg[exp_name]['directories']['dirCalibrated'] = dir_cal
        in_cfg[exp_name]['directories']['dirFinal'] = dir_final
        in_cfg[exp_name]['directories']['dirGraphics'] = dir_graphics

    # -------------------------------------------------------------------------
    # Integrity of the location library

    # Determine if a file suffix was provided.
    try:
        in_cfg['outname_suffix'] = cfg['directories']['suffix']
    except KeyError:
        in_cfg['outname_suffix'] = None

   # Determine fiber type
    try:
        in_cfg['coretype'] = in_cfg[exp_name]['dataProperties']['fiber_type']
    except KeyError:
        in_cfg['coretype'] = 'singlecore'
    assert (in_cfg['coretype'] == 'singlecore' or
            in_cfg['coretype'] == 'multicore')

    if in_cfg['coretype'] == 'multicore':
        in_cfg['cores'] = in_cfg[exp_name]['dataProperties']['cores']
    else:
        in_cfg['cores'] = None

    # Determine the list of location types. If none is provided, find all
    # unique location types listed in the location library.
    try:
        loc_type = cfg['dataProperties']['all_locs']
        check_loc_library = True
    except KeyError:
        loc_type = []
        try:
            for l in cfg['location_library']:
                loc_type.append(cfg['location_library'][l]['loc_type'])
            loc_type = np.unique(loc_type).tolist()
            check_loc_library = True
        except KeyError:
            print('No location library found, all steps after creating raw' +
                  ' netcdfs have been disabled.')
            check_loc_library = False
            cfg['flags']['calibrate_flag'] = False
            cfg['flags']['processed_flag'] = False
            cfg['flags']['final_flag'] = False
        except TypeError:
            print('No location library found, all steps after creating raw' +
                  ' netcdfs have been disabled.')
            check_loc_library = False
            cfg['flags']['calibrate_flag'] = False
            cfg['flags']['processed_flag'] = False
            cfg['flags']['final_flag'] = False

    # -------------------------------------------------------------------------
    # Labeling sections and physical coordinates
    missing_mess = 'Coordinates for {ploc} were not defined by pairs of LAF, x, y, and z in the location library.'

    if check_loc_library and in_cfg['flags']['final_flag']:
        try:
            phys_locs = config_user['dataProperties']['phys_locs']
        except KeyError as kerr:
            print('The phys_locs field in dataProperties is empty but labeling physical coordinates is turned on.')
            raise kerr

    if check_loc_library and in_cfg['coretype'] == 'singlecore':
        for loc_type_cur in loc_type:
            lib[loc_type_cur] = {}
            for l in cfg['location_library']:
                if loc_type_cur == cfg['location_library'][l]['loc_type']:
                    lib[loc_type_cur][l] = cfg['location_library'][l]

                    # All locations are defined by a pair of LAFs.
                    assert 'LAF' in lib[loc_type_cur][l], 'No LAFs provided for ' + l
                    assert len(lib[loc_type_cur][l]['LAF']) == 2, missing_mess.format(loc_type_cur)

                    # All physical coordinates are defined by pairs of x, y, z
                    if (cfg['flags']['final_flag']
                            and loc_type_cur in phys_locs):
                        assert 'x' in lib[loc_type_cur][l], missing_mess.format(loc_type_cur)
                        assert 'y' in lib[loc_type_cur][l], missing_mess.format(loc_type_cur)
                        assert 'z' in lib[loc_type_cur][l], missing_mess.format(loc_type_cur)

                        assert len(lib[loc_type_cur][l]['x']) == 2, missing_mess.format(loc_type_cur)
                        assert len(lib[loc_type_cur][l]['y']) == 2, missing_mess.format(loc_type_cur)
                        assert len(lib[loc_type_cur][l]['z']) == 2, missing_mess.format(loc_type_cur)

            if not bool(lib[loc_type_cur]):
                print(loc_type_cur + ' was not found in location library.')

        # Asserts that lib is an iterable list with stuff in it.
        mess = 'No items found in location library.'
        # Returns False if it is emtpy, True if is not.
        assert bool(lib), mess

    # Split up multicore data
    elif check_loc_library and in_cfg['coretype'] == 'multicore':
        for c in cores:
            if c in cfg['cores']:
                mess = ('core LAF limits for was poorly formatted. Expected'
                        ' a list of two LAF values for each core.')
                assert len(cfg['cores'][c]['LAF']) == 2, mess
            else:
                print('Expected to find field "LAF" for each core listed in ' +
                      'dataProperties of config file.')
                raise KeyError

            for loc_type_cur in loc_type:
                lib[loc_type_cur] = {}
                for l in cfg['location_library']:
                    if loc_type_cur == cfg['location_library'][l]['loc_type']:
                        lib[loc_type_cur][c][l] = cfg['location_library'][l]['LAF'][c]

    # -------------------------------------------------------------------------
    # Further variables to assess.

    # Step loss corrections
    if ('step_loss_LAF' in cfg['dataProperties']
            and 'step_loss_correction' in cfg['dataProperties']):
        step_loss_flag = True
        in_cfg['step_loss'] = {}
        in_cfg['step_loss']['flag'] = step_loss_flag
        in_cfg['step_loss']['LAF'] = cfg['dataProperties']['step_loss_LAF']
        in_cfg['step_loss']['correction'] = cfg['dataProperties']['step_loss_correction']
    else:
        step_loss_flag = False
        in_cfg['step_loss'] = {}
        in_cfg['step_loss']['flag'] = step_loss_flag

    # Determine write mode:
    try:
        write_mode = cfg['flags']['write_mode']
    except KeyError:
        write_mode = 'overwrite'
    # Add check here that write_mode has an expected value.
    in_cfg['write_mode'] = write_mode

    # The calibration mode will be changed in an upcoming update.
    # # Determine calibration mode:
    # try:
    #     cal_mode = cfg['flags']['cal_mode']
    # except KeyError:
    #     cal_mode = 'instantaneous'

    # Channels to process -- this should be on the outside of the script
    channelNames = cfg['directories']['channelName']

    return in_cfg, lib

    # Items to unpack in pyfocs
    # experiment_names
    # coretype
    # cores
    # core LAFs
    # outname_suffix
    # step_loss_flag and LAF
    # write_mode
    # channelNames
    # cal_mode
    # return of None -- custom error message.
