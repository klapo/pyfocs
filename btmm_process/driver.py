import os
import dtsarch


def run_btmm_process(configFile, archive_flag=True,
                     tar_read_flag=True, calibrateFlag=True):

    # Load the configuration file
    cfg = yamlDict(configFile)

    # Archive the raw xml files and move them to the desired locations
    if archive_flag:
        # Call the archiving function
        dtsarch.archiver(cfg)

    # Call readDTS.tar_read
    if tar_read_flag:
        archive_read(cfg)

    # Call the calibration routine
    if calibrateFlag:
        prepCalibrate(cfg)
