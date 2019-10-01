# This script is the light weight version of PyFOX intended to be used for
# archiving DTS data in real time on the Ultimas.

#from IPython import get_ipython #clearing variables that python still knows from before
#get_ipython().magic('reset -sf')

import numpy as np
from datetime import timedelta
import pandas as pd
import xarray as xr

# OS interaction
import os
import yaml
import copy

# UBT's package for handling dts data
import btmm_process

# Ignore the future compatibility warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)

#%% open the config file
# The filename is hard coded to make it easier to use the
# Windows Task Scheduler.
filename_configfile = '/Users/karllapo/Desktop/software/python/scripts/PyFOX/config_silixa_archiver.yml' 

with open(filename_configfile, 'r') as stream:
    config_user = yaml.load(stream)

#%% create directories
dir_original = config_user['directories']['dir_pre']
if not os.path.exists(dir_original): # give an error if raw data folder can't be found
    print('Raw data folder ' + dir_original + ' does not exist.')

# Archive location
dir_archive = config_user['directories']['dir_archive']

# find all experiments to be processed and make a list of them
os.chdir(dir_original)
contents = os.listdir()

# The CR6 part will not be necessary for flyfox but leave it in here for later
DTS_folders = config_user['directories']['dir_exps']

# Loop through all of the DTS data directories
# assemble internal config for each dts folder within the experiment folder
dir_raw_data = {}
dir_proc_data = {}
internal_config = {}
for dtsf in np.atleast_1d(DTS_folders):
    # get all the directories within this experiment folder
    dir_raw_data[dtsf] = os.path.join(dir_original, dtsf)
    dir_proc_data[dtsf] = os.path.join(dir_archive, dtsf)

for dtsf in np.atleast_1d(DTS_folders):
    internal_config[dtsf] = {}
    internal_config[dtsf] = copy.deepcopy(config_user)
    internal_config[dtsf]['archive']['sourcePath'] = dir_raw_data[dtsf]
    internal_config[dtsf]['archive']['targetPath'] = dir_proc_data[dtsf]
    internal_config[dtsf]['fileName']['filePrefix'] = dtsf
    if not os.path.exists(internal_config[dtsf]['archive']['targetPath']): # create the folder for the processed files if it doesn't already exist
        os.makedirs(internal_config[dtsf]['archive']['targetPath'])

# archiving
if config_user['flags']['archiving_flag']:
    for dtsf in np.atleast_1d(DTS_folders):
        print('archiving ', dtsf)
        btmm_process.archiver(internal_config[dtsf])
