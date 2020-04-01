import yaml


def config(cfg):
    '''
    Function checks the basic integrity of a configuration file and returns the
    config dictionary, location library, and physical locations dictionary
    expected by pyfocs.
    '''

    # The internal config file and location library are handed back to pyfocs
    internal_config = {}
    lib = {}

    # This yaml syntax requires pyyaml v5.1
    try:
        cfg = pyfocs.yamlDict(example_config)
    except yaml.parser.ParserError as exc:
        print('The config yaml file is not properly formatted.')
        print('Block mapping errors will return bad lines to examine:')
        print(exc)
        # Example output:
        # while parsing a block mapping
        # in "fine_name.yml", line 2, column 1

    try:
        remote_dir = config_user['directories']['remote']['remote_directories']
    except KeyError:
        remote_dir = []

    # -------------------------------------------------------------------------
    # Look through the configuration file and prepare
    # -------------------------------------------------------------------------
    #%% loop through the experiment_names
    experiment_names = config_user['directories']['experiment_names']
    assert isinstance(config_user['directories']['experiment_names'], str)

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

        # --------
        # Archived tar.gz of xmls
        if 'archive' in remote_dir:
            dir_archive = os.path.join(dir_pre_remote,
                                       exp_name, config_user['directories']['archive'])
        else:
            dir_archive = os.path.join(dir_pre_local,
                                       exp_name, config_user['directories']['archive'])

        if (not os.path.exists(dir_archive)):
            # create the folder for the archived files if it doesn't already
            # exist and archived files will be created
            if config_user['flags']['archiving_flag']:
                os.makedirs(dir_archive)
                print('Archived directory not found. Created a new one at ' + dir_archive)
            # throw error if archived files are needed but neither found nor created
            elif config_user['flags']['archive_read_flag']:
                raise FileNotFoundError('Archived data directory not found at ' +
                                        dir_archive)

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
            # throw error if raw netcdf files are needed but neither found
            # nor created
            elif config_user['flags']['calibrate_flag']:
                raise FileNotFoundError('Raw netcdf data directory not found ' +
                                        dir_raw_netcdf)

        # --------
        # Calibrated temperature
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
        if ((not os.path.exists(dir_final))
                and ((config_user['flags']['calibrate_flag'])
                     or (config_user['flags']['final_flag']))):
            os.makedirs(dir_final)
            print('Final directory not found. Created a new one at ' + dir_final)

        # --------
        # Graphics folder
        if 'graphics' in remote_dir:
            dir_graphics = os.path.join(dir_pre_remote,
                                        exp_name, config_user['directories']['graphics'])
        else:
            dir_graphics = os.path.join(dir_pre_local,
                                        exp_name, config_user['directories']['graphics'])

        # assemble internal config for each dts folder within the
        # experiment folder
        internal_config[exp_name] = {}
        internal_config[exp_name] = copy.deepcopy(config_user)
        internal_config[exp_name]['archive']['channelName'] = config_user['directories']['channelName']
        internal_config[exp_name]['archive']['sourcePath'] = dir_raw_xml
        internal_config[exp_name]['archive']['targetPath'] = dir_archive
        internal_config[exp_name]['directories']['dirArchive'] = dir_archive
        internal_config[exp_name]['directories']['dirRawNetcdf'] = dir_raw_netcdf
        internal_config[exp_name]['directories']['dirCalibrated'] = dir_cal
        internal_config[exp_name]['directories']['dirFinal'] = dir_final
        internal_config[exp_name]['directories']['dirGraphics'] = dir_graphics

    # -------------------------------------------------------------------------
    # Integrity of the location library
    # -------------------------------------------------------------------------
    # Determine if a file suffix was provided.
    try:
        internal_config['outname_suffix'] = config_user['directories']['suffix']
    except KeyError:
        internal_config['outname_suffix'] = None

   # Determine fiber type
    try:
        internal_config['coretype'] = internal_config[exp_name]['dataProperties']['fiber_type']
    except KeyError:
        internal_config['coretype'] = 'singlecore'
    assert (internal_config['coretype'] == 'singlecore' or
            internal_config['coretype'] == 'multicore')

    if coretype == 'multicore':
        internal_config['cores'] = internal_config[exp_name]['dataProperties']['cores']
    else:
        internal_config['cores'] = None

    # Determine the list of location types. If none is provided, find all
    # unique location types listed in the location library.
    try:
        loc_type = config_user['dataProperties']['all_locs']
        check_loc_library = True
    except KeyError:
        loc_type = []
        try:
            for l in config_user['location_library']:
                loc_type.append(config_user['location_library'][l]['loc_type'])
            loc_type = np.unique(loc_type).tolist()
            check_loc_library = True
        except KeyError:
            print('No location library found, all steps after creating raw' +
                  ' netcdfs have been disabled.')
            check_loc_library = False
            config_user['flags']['calibrate_flag'] = False
            config_user['flags']['processed_flag'] = False
            config_user['flags']['final_flag'] = False
        except TypeError:
            print('No location library found, all steps after creating raw' +
                  ' netcdfs have been disabled.')
            check_loc_library = False
            config_user['flags']['calibrate_flag'] = False
            config_user['flags']['processed_flag'] = False
            config_user['flags']['final_flag'] = False

    # The physical location stuff is broken here.
    # Make sure the location library is formatted correctly and has stuff in
    # itself.
    if coretype == 'singlecore':
        for loc_type_cur in loc_type:
            lib[loc_type_cur] = {}
            for l in config_user['location_library']:
                if loc_type_cur == config_user['location_library'][l]['loc_type']:
                    try:
                        lib[loc_type_cur][l] = config_user['location_library'][l]['LAF']
                    except KeyError:
                        print('No LAFs provided for ' + l)
                    assert np.size(lib[loc_type_cur][l]) == 2
            else:
                print(loc_type_cur ' was not found in location library.')
        # Asserts that lib is an iterable list with stuff in it.
        assert self.assertTrue(len(lib))

    # Split up multicore data
    elif coretype == 'multicore':
        for c in cores:
            if c in config['cores']:
                mess = ('core LAF limits for was poorly formatted. Expected'
                        ' a list of two LAF values for each core.')
                assert len(config['cores'][c]['LAF']) == 2, mess
            else:
                print('Expected to find field "LAF" for each core listed in ' +
                      'dataProperties of config file.')
                raise KeyError

            for loc_type_cur in loc_type:
                lib[loc_type_cur] = {}
                for l in config_user['location_library']:
                    if loc_type_cur == config_user['location_library'][l]['loc_type']:
                        lib[loc_type_cur][l] = config_user['location_library'][l]['LAF'][c]

    # -------------------------------------------------------------------------
    # Physical coordinates

    # These location types are to be labeled
    phys_locs = config_user['dataProperties']['phys_locs']

    if coretype == 'multicore':
        location[l]['LAF'] = config_user['location_library'][l]['LAF'][core]
    elif coretype == 'singlecore':
        location[l]['LAF'] = config_user['location_library'][l]['LAF']

    # -------------------------------------------------------------------------
    # Further variables to assess.

    # Step loss corrections
    if ('step_loss_LAF' in config_user['dataProperties']
            and 'step_loss_correction' in config_user['dataProperties']):
        step_loss_flag = True
        internal_config['step_loss']['flag'] = step_loss_flag
        internal_config['step_loss']['LAF'] = config_user['dataProperties']['step_loss_LAF']
        internal_config['step_loss']['correction'] = config_user['dataProperties']['step_loss_correction']
    else:
        step_loss_flag = False
        internal_config['step_loss']['flag'] = step_loss_flag

    # Determine write mode:
    try:
        write_mode = config_user['flags']['write_mode']
    except KeyError:
        write_mode = 'overwrite'
    # Add check here that write_mode has an expected value.
    internal_config['write_mode'] = write_mode

    # The calibration mode will be changed in an upcoming update.
    # # Determine calibration mode:
    # try:
    #     cal_mode = config_user['flags']['cal_mode']
    # except KeyError:
    #     cal_mode = 'instantaneous'

    # Channels to process
    channelNames = config_user['directories']['channelName']

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
