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

        return None, None
        # Example output:
        # while parsing a block mapping
        # in "fine_name.yml", line 2, column 1

    # -------------------------------------------------------------------------
    # Look through the configuration file and prepare
    experiment_names = cfg['directories']['experiment_names']
    mess = 'experiment_names should be a list (e.g., all items start with a dash)'
    assert isinstance(cfg['directories']['experiment_names'], list), experiment_names
    in_cfg['experiment_names'] = experiment_names

    # Verify that the flags exist and are booleans.

    # Remove deprecated flags.
    try:
        del cfg['flags']['cal_mode']
        in_cfg['flags'] = cfg['flags']
    except KeyError:
        in_cfg['flags'] = cfg['flags']

    # Catch an old name for the "final" flag.
    if 'processed_flag' in in_cfg['flags']:
        in_cfg['flags']['final_flag'] = in_cfg['flags']['processed_flag']
        del in_cfg['processed_flag']

    flags = ['ref_temp_option', 'write_mode', 'archiving_flag',
             'archive_read_flag', 'calibrate_flag', 'final_flag']
    if not all([fl in flags for fl in in_cfg['flags']]):
        missing_flags = [fl for fl in in_cfg['flags'] if fl not in flags]
        raise KeyError('Not all flags were found:\n' + '\n'.join(missing_flags))

    # --------
    # Check paths that do not vary between experiments.

    # Directories located locally.
    try:
        dir_pre_local = cfg['directories']['local']['dir_pre']
    except KeyError:
        dir_pre_local = ''

    # List of directories located remotely.
    try:
        remote_dir = cfg['directories']['remote']['remote_directories']
    except KeyError:
        remote_dir = []
    # Remote path for running on a server.
    try:
        dir_pre_remote = cfg['directories']['remote']['dir_pre']
    except KeyError:
        dir_pre_remote = ''

    # Handling for relative paths
    if dir_pre_local == '.':
        dir_pre_local = os.getcwd()
    if dir_pre_remote == '.':
        dir_pre_remote = os.getcwd()

    # External data directory
    try:
        dir_ext = cfg['directories']['external']
        # Handling for relative paths
        if dir_ext == '.':
            dir_ext = os.getcwd()
        # Check for empty yaml lists (assigned NoneType)
        if dir_ext is not None:
            # check if the external data directory exists, error if not
            if (not os.path.exists(dir_ext)) and (cfg['flags']['ref_temp_option']):
                raise FileNotFoundError('External data directory not found at\n'
                                        + dir_ext)
    except KeyError:
        dir_ext = None

    # External data file
    if in_cfg['flags']['ref_temp_option'] == 'external':
        try:
            ext_fname = os.path.join(dir_ext,
                                     cfg['directories']['filename_external'])
            os.path.isfile=(ext_fname)
        except KeyError:
            mess = ('Config file indicates to use external data, but no file '
                    'was found at:\n' + ext_fname)
            raise FileNotFoundError(mess)

        # Make sure the file is a netcdf
        if not ext_fname.endswith('.nc'):
            mess = ('Config file indicates to use external data, but a '
                    'netcdf was not found at:\n' + ext_fname)
            raise FileNotFoundError()
        # Build the path to the file.
        in_cfg['external_data'] = os.path.join(dir_ext, ext_fname)

    # Prepare calibration relevant parameters.
    if in_cfg['flags']['calibrate_flag']:
        # Resampling time
        try:
            dt = cfg['dataProperties']['resampling_time']
        except KeyError:
            mess = 'A resampling time must be provided when calibrating.'
            raise KeyError(mess)
        in_cfg['resampling_time'] = dt

        # Calibration probe names
        if in_cfg['flags']['ref_temp_option'] == 'external':
            try:
                probe1_name = cfg['dataProperties']['probe1Temperature']
                probe2_name = cfg['dataProperties']['probe2Temperature']
            except KeyError:
                mess = ('When using external reference temperature probes '
                        'the probe1Temperature and probe2Temperature fields '
                        'must be included.')
                raise KeyError(mess)

            # Determine if additional external fields are to be added.
            try:
                in_cfg['external_fields'] = cfg['dataProperties']['external_fields']
            except KeyError:
                in_cfg['external_fields'] = []

        else:
            probe1_name = cfg['dataProperties']['probe1Temperature']
            probe2_name = cfg['dataProperties']['probe2Temperature']
        # Calibration field
        in_cfg['calibration'] = copy.deepcopy(cfg['calibration'])
        # The multiiple labelings of the reference probe will be removed with
        # the coming update to the calibration.
        in_cfg['probe1'] = probe1_name
        in_cfg['probe2'] = probe2_name

    # --------
    # Prepare each requested experiment.
    for exp_name in experiment_names:
        # Create directories if they don't exist and
        # errors if needed data are not found.

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
            mess = ('Creating archived data is on but the raw xml data '
                    'directory was not found at ' + dir_raw_xml)
            raise FileNotFoundError(mess)

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
                mess = ('Reading archived data is on but the archived data '
                        'directory was not found at ' + dir_archive)
                raise FileNotFoundError(mess)

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
                mess = ('Calibration is on but the raw netcdf data '
                        'directory was not found at ' + dir_raw_netcdf)
                raise FileNotFoundError(mess)
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
                mess = ('Finalizing is on but the calibrated data '
                        'directory was not found at ' + dir_cal)
                raise FileNotFoundError(mess)

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
        in_cfg[exp_name]['probe1'] = probe1_name
        in_cfg[exp_name]['probe2'] = probe2_name


    # -------------------------------------------------------------------------
    # Integrity of the location library

    # Determine if a file suffix was provided.
    try:
        in_cfg['outname_suffix'] = cfg['directories']['suffix']
    except KeyError:
        in_cfg['outname_suffix'] = None

   # Determine fiber type
    try:
        in_cfg['coretype'] = cfg['dataProperties']['fiber_type']
    except KeyError:
        in_cfg['coretype'] = 'singlecore'
    assert (in_cfg['coretype'] == 'singlecore' or
            in_cfg['coretype'] == 'multicore')

    if in_cfg['coretype'] == 'multicore':
        in_cfg['cores'] = cfg['dataProperties']['cores']
        for c in in_cfg['cores']:
            assert len(in_cfg['cores'][c]['LAF']) == 2, 'Cores property must include LAF pair of core limits.'
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
            in_cfg['flags']['calibrate_flag'] = False
            in_cfg['flags']['final_flag'] = False
        except TypeError:
            print('No location library found, all steps after creating raw' +
                  ' netcdfs have been disabled.')
            check_loc_library = False
            in_cfg['flags']['calibrate_flag'] = False
            in_cfg['flags']['final_flag'] = False
    if 'calibration' not in loc_type and calibrate_flag:
        loc_type.append('calibration')

    # -------------------------------------------------------------------------
    # Labeling sections and physical coordinates
    if check_loc_library and in_cfg['flags']['final_flag']:
        try:
            phys_locs = cfg['dataProperties']['phys_locs']
            in_cfg['phys_locs'] = phys_locs
        except KeyError as kerr:
            mess = ('The phys_locs field in dataProperties is empty but '
                    'labeling physical coordinates is turned on.')
            print(mess)
            raise kerr

    # Error messages for the assert statements.
    missing_mess = ('Coordinates for {ploc} were not defined by pairs of LAF, '
                    'x, y, and z in the location library. Found intsead:\n'
                    '{lib}')
    missing_mess_mc = ('Coordinates for core {ploc} {co} were not defined '
                       'by a pairs of LAF in the location library. Found '
                       'instead:\n {lib}'
                       )

    # Check the integrity of a single core location library.
    if check_loc_library and in_cfg['coretype'] == 'singlecore':
        for loc_type_cur in loc_type:
            lib[loc_type_cur] = {}
            for l in cfg['location_library']:
                if loc_type_cur == cfg['location_library'][l]['loc_type']:
                    lib[loc_type_cur][l] = cfg['location_library'][l]

                    # All locations are defined by a pair of LAFs.
                    assert 'LAF' in lib[loc_type_cur][l], 'No LAFs provided for ' + l
                    assert len(lib[loc_type_cur][l]['LAF']) == 2, missing_mess.format(ploc=l, lib=lib[loc_type_cur][l])
                    if any(np.isnan(lib[loc_type_cur][l]['LAF'])):
                        del lib[loc_type_cur][l]
                        print(l + ' will not be labeled due to NaNs in LAF. Check library.')
                        continue

                    # All physical coordinates are defined by pairs of x, y, z
                    if (cfg['flags']['final_flag']
                            and loc_type_cur in phys_locs):
                        assert 'x_coord' in lib[loc_type_cur][l], missing_mess.format(ploc=l, lib=lib[loc_type_cur][l])
                        assert 'y_coord' in lib[loc_type_cur][l], missing_mess.format(ploc=l, lib=lib[loc_type_cur][l])
                        assert 'z_coord' in lib[loc_type_cur][l], missing_mess.format(ploc=l, lib=lib[loc_type_cur][l])

                        assert len(lib[loc_type_cur][l]['x_coord']) == 2, missing_mess.format(ploc=l, lib=lib[loc_type_cur][l])
                        assert len(lib[loc_type_cur][l]['y_coord']) == 2, missing_mess.format(ploc=l, lib=lib[loc_type_cur][l])
                        assert len(lib[loc_type_cur][l]['z_coord']) == 2, missing_mess.format(ploc=l, lib=lib[loc_type_cur][l])

            if not bool(lib[loc_type_cur]):
                print(loc_type_cur + ' was not found in location library.')

        # Asserts that lib is an iterable list with stuff in it.
        mess = 'No items found in location library.'
        # Returns False if it is emtpy, True if is not.
        assert bool(lib), mess
        # Make sure that calibration labels are available
        assert 'calibration' in lib, 'Calibration location type is missing'

    # Check the integrity of a multicore location library.
    elif check_loc_library and in_cfg['coretype'] == 'multicore':
        try:
            mess = ('core LAF limits for was poorly formatted. Expected'
                    ' a list of two LAF values for each core.')
            assert len(in_cfg['cores'][c]['LAF']) == 2, mess
        except KeyError:
            mess = ('Expected to find field "LAF" for {co} listed in '
                    'dataProperties of config file.')
            raise KeyError(mess.format(co=c))

        # There is a difference between what makes sense to pass to pyfocs
        # and what is human readable. To handle that discrepancy we
        # unfortunately need to repack the dictionary.
        for c in in_cfg['cores']:
            lib[c] = {}

            for loc_type_cur in loc_type:
                lib[c][loc_type_cur] = {}

                for l in cfg['location_library']:
                    if loc_type_cur == cfg['location_library'][l]['loc_type']:
                        lib_this_loc = copy.deepcopy(cfg['location_library'][l])
                        # Debugging lines, delete after more testing.
                        # print(l)
                        # print(lib_this_loc)
                        # print('---------------------------------------')
                        # Drop any other other cores from the LAF field.
                        lib[c][loc_type_cur][l] = lib_this_loc
                        lib[c][loc_type_cur][l]['LAF'] = lib[c][loc_type_cur][l]['LAF'][c]

                        # All locations are defined by a pair of LAFs.
                        assert 'LAF' in lib[c][loc_type_cur][l], 'No LAFs provided for ' + l
                        assert len(lib[c][loc_type_cur][l]['LAF']) == 2, missing_mess_mc.format(ploc=l, co=c, lib=lib)

                        # Remove labels that do not have an LAF
                        if any(np.isnan(lib[c][loc_type_cur][l]['LAF'])):
                            del lib[c][loc_type_cur][l]
                            print(l + ' will not be labeled due to NaNs in LAF.')
                            continue

                        # All physical coordinates are defined by pairs of x, y, z
                        if (cfg['flags']['final_flag']
                                and loc_type_cur in phys_locs):
                            assert 'x_coord' in lib[c][loc_type_cur][l], missing_mess.format(ploc=l, lib=lib[loc_type_cur][l])
                            assert 'y_coord' in lib[c][loc_type_cur][l], missing_mess.format(ploc=l, lib=lib[loc_type_cur][l])
                            assert 'z_coord' in lib[c][loc_type_cur][l], missing_mess.format(ploc=l, lib=lib[loc_type_cur][l])

                            assert len(lib[c][loc_type_cur][l]['x_coord']) == 2, missing_mess.format(ploc=l, lib=lib[loc_type_cur][l])
                            assert len(lib[c][loc_type_cur][l]['y_coord']) == 2, missing_mess.format(ploc=l, lib=lib[loc_type_cur][l])
                            assert len(lib[c][loc_type_cur][l]['z_coord']) == 2, missing_mess.format(ploc=l, lib=lib[loc_type_cur][l])

            assert 'calibration' in lib[c], 'Calibration location type is missing for core {c}'.format(c=c)

    # -------------------------------------------------------------------------
    # Further variables to assess.

    # Step loss corrections
    if ('step_loss_LAF' in cfg['dataProperties']
            and 'step_loss_correction' in cfg['dataProperties']):
        in_cfg['step_loss'] = {}
        steploss_laf = np.atleast_1d(cfg['dataProperties']['step_loss_LAF']).astype('float')
        steploss_corr = np.atleast_1d(cfg['dataProperties']['step_loss_correction']).astype('float')
        for sl_laf, sl_co in zip(steploss_laf, steploss_corr):
            # Make sure the step loss corrections are numeric
            if np.isnan(sl_laf) or np.isnan(sl_co):
                step_loss_flag = False
                in_cfg['step_loss'] = {}
                break

            else:
                step_loss_flag = True
        in_cfg['step_loss']['flag'] = step_loss_flag
        in_cfg['step_loss']['LAF'] = steploss_laf
        in_cfg['step_loss']['correction'] = steploss_corr
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

    # Channels to process -- this should be on the outside of the script
    in_cfg['channelNames'] = cfg['directories']['channelName']

    return in_cfg, lib

    # Add:
    # [x] Resampling time
    # [x] External data directory/filename
    # [x] Probe names for calibration

    # Items to unpack in pyfocs
    # [x] physical locations
    # [x] experiment_names
    # [x] coretype
    # [x] cores
    # [?] core LAFs
    # [x] outname_suffix
    # [x] step_loss_flag, LAF, and corrections
    # [x] write_mode
    # [x] channelNames
    # [ ] cal_mode -> ignore for now until the new calibration method works.
    # [x] return of None -- custom error message.
