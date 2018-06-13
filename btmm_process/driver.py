import yaml
import os
from .labeler import yamlDict
from .readDTS import tar_read
from .calibrate import prepCalibrate
from .dtsarch import dtsarch

def run_btmm_process(configFile, archive_flag=True, tar_read_flag=True, calibrateFlag=True):

    # Load the configuration file
    cfg = yamlDict(configFile)

    # Archive the raw xml files and move them to the desired locations
    if archive_flag:
        # Read the configure file and pass the arguments to the archiver
        mode = cfg['archive']['mode']
        channel_1 = cfg['archive']['channel_1']
        channel_2 = cfg['archive']['channel_2']
        channel_3 = cfg['archive']['channel_3']
        channel_4 = cfg['archive']['channel_4']
        channels = [channel_1, channel_2, channel_3, channel_4]
        sourcePath = cfg['archive']['sourcePath']
        targetPath = cfg['archive']['targetPath']
        
        # Determine if we are backing up to an external directory
        dirBackUp = cfg['archive']['dirBackUp']
        if os.path.isdir(dirMobile):
            externalBackUp = True
            logfile = cfg['archive']['logfile']
        else:
            externalBackUp = False
            logfile = []
        
        # Call the archiving function
        dtsarch.archiver(sourcePath, targetPath, channels, mode,
                         externalBackUp=externalBackUp,
                         dirBackUp=dirBackUp, logfile=logifle):

    # Call readDTS.tar_read
    if tar_read_flag:
        tar_read(cfg)

    # Call the calibration routine
    if calibrateFlag:
        prepCalibrate(cfg)
