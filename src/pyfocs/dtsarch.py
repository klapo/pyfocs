import os
import datetime
import tarfile
import numpy as np


def make_tarfile(tarName, filesToZip):
    '''
    This helper function creates the command to tar.gz the DTS XML files.
    INPUT:
        tarName = Name of the tar.gz file to create
        filesToZip = List with the name of the source files to compress.
    '''

    print('Archiving ' + tarName)

    # Open the tarball and prepare it for writing
    # and compress the xml files
    with tarfile.open(tarName, "w:gz") as tar:
        for f in filesToZip:
            tar.add(f, arcname=os.path.basename(f))


def round_dt(dt, delta, type='floor'):
    '''
    Rounds a datetime object to the nearest minute delta
    dt - datetime object
    delta - rounding interval in minutes (e.g., nearest 15, delta=15)
    type - specifies if we round up or down
    '''
    if type == 'floor':
        dt = dt - datetime.timedelta(minutes=dt.minute % delta,
                                     seconds=dt.second,
                                     microseconds=dt.microsecond)
    elif type == 'ceil':
        dt = dt + datetime.timedelta(minutes=delta - dt.minute % delta - 1,
                                     seconds=60 - dt.second - 1,
                                     microseconds=10**6 - dt.microsecond)
    else:
        raise ValueError('Unrecognized dt rounding type. Valid options are'
                         'floor and ceil.')
    return(dt)


# ------------------------------------------------------------------------------
# Formats a datetime object into its components pieces
def dt_strip(dt, str_convert=False):
    year = dt.year
    month = dt.month
    day = dt.day
    hour = dt.hour
    minute = dt.minute
    sec = dt.second
    msec = dt.microsecond

    if str_convert:
        year = str(dt.year)

        if month < 10:
            month = '0' + str(month)
        else:
            month = str(month)

        if day < 10:
            day = '0' + str(day)
        else:
            day = str(day)

        if hour < 10:
            hour = '0' + str(hour)
        else:
            hour = str(hour)

        if minute < 10:
            minute = '0' + str(minute)
        else:
            minute = str(minute)

        if sec < 10:
            sec = '0' + str(sec)
        else:
            sec = str(sec)

        if msec < 10:
            msec = '0' + str(msec)
        else:
            msec = str(msec)

    return(year, month, day, hour, minute, sec, msec)
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Formats a datetime object into its components pieces or a string
def dt_string_label(t):
    '''
    t - string of dts filename (with a date/time indicated) to be converted
        to an actual python datetime.
    '''
    if 'UTC' not in t:
        # This is assumed to be an XT file if it lacks 'UTC'.
        # The format looks like: channel 1_20171109203321986.xml
        t = t.split('_')[-1]
        t = t.split('.')[0]
        year = t[0:4]
        month = t[4:6]
        day = t[6:8]
        hour = t[8:10]
        minute = t[10:12]
        sec = t[12:14]
        if len(t) > 14:
            msec = t[14:16]
    else:
        t = t.split('_')[-2:]
        t_ymd = t[0]
        t = t[-1].split('.')[0:2]
        t_hms = t[0]
        t_us = t[1]
        year = t_ymd[0:4]
        month = t_ymd[4:6]
        day = t_ymd[6:8]
        hour = t_hms[0:2]
        minute = t_hms[2:4]
        sec = t_hms[4:6]
        msec = t_us
    out_dt = datetime.datetime(int(year), int(month),
                               int(day), int(hour),
                               int(minute), int(sec),
                               int(msec) * 1000)
    return(out_dt)
# ------------------------------------------------------------------------------


def archiver(cfg):
    '''
    Archive data
    Script to tar and gzip the Ultima data.
    Requires a configuration file to run.
    '''

    # Read the configure file for the archiver arguments
    mode = cfg['archive']['mode']
    channels = np.atleast_1d(cfg['archive']['channelName'])
    sourcePath = cfg['archive']['sourcePath']
    targetPath = cfg['archive']['targetPath']
    delta_minutes = cfg['archive']['archiveInterval']

    # Determine if the archiver should clean up the original data directory
    try:
        cleanup_flag = cfg['archive']['cleanup_flag']
    except KeyError:
        cleanup_flag = False

    for ch in channels:
        channelPath = os.path.join(sourcePath, ch)
        print(channelPath)

        ########
        # Archive data to .tar.gz files.
        # Check if the channel directory exists.
        if os.path.isdir(channelPath):
            contents = os.listdir(channelPath)
        # If it doesn't, move on to the next channel
        else:
            print('Did not find ' + ch)
            continue
        # Only select xml files
        contents = [c for c in contents if '.xml' in c]
        # Sort the file list alphabetically.
        contents.sort()

        # Verify that xml files exist before proceeding.
        if len(contents) == 0:
            print('No xml files found in ' + channelPath)
            break

        # First datetime of the raw data
        t = contents[0]
        dtInit = dt_string_label(t)

        # Last datetime of the raw data
        t = contents[-1]
        dtFinal = dt_string_label(t)

        # First time step interval
        dt1 = round_dt(dtInit, delta_minutes, type='floor')
        dt2 = round_dt(dtInit, delta_minutes, type='ceil')

        # Round up to the nearest interval for the end
        dtFinal = round_dt(dtFinal, delta_minutes, type='ceil')
        # Empty container for the files in this interval
        interval_contents = []
        xml_counts = 0
        # Iterate through the (date) sorted list of raw xml files
        while xml_counts <= len(contents):
            # Split the file name string into the datetime components
            c = contents[xml_counts]
            dt = dt_string_label(c)

            # Determine if we are still within this interval.
            if dt < dt2 and dt >= dt1:
                interval_contents.append(os.path.join(sourcePath, ch, c))
                xml_counts = xml_counts + 1

            # We have spanned this interval, save the data and move on
            elif dt >= dt2:
                # Create file names for this hour
                year, month, day, hour, minute, _, _ = dt_strip(dt1, str_convert=True)

                dateFileName = '_' + year + month + day + '-' + hour + minute
                outFile = os.path.join(targetPath, ch + dateFileName + '.tar.gz')

                # Check if any files fall within the interval
                if len(interval_contents) == 0:
                    print('No xml files in the interval: ' + str(dt1))
                else:
                    make_tarfile(outFile, interval_contents)
                dt1 = dt2
                dt2 = dt2 + datetime.timedelta(minutes=delta_minutes)

                # Determine if the xml files should be removed
                if cleanup_flag:
                    sourceFile = interval_contents
                    print('Cleaning up the raw xml files...')
                    for f in sourceFile:
                        os.remove(f)

                interval_contents = []

            # We reached the end of the file list, save and exit.
            if xml_counts == len(contents):

                if mode == 'archiving':
                    # Create file names for this hour
                    year, month, day, hour, minute, _, _ = dt_strip(dt1, str_convert=True)

                    dateFileName = '_' + year + month + day + '-' + hour + minute
                    outFile = os.path.join(targetPath, ch + dateFileName + '.tar.gz')
                    sourceFile = interval_contents
                    make_tarfile(outFile, sourceFile)

                    # Determine if the xml files should be removed
                    if cleanup_flag:
                        print('Cleaning up the raw xml files...')
                        for f in sourceFile:
                            os.remove(f)

                    print('Done with ' + ch + '.')
                    xml_counts = xml_counts + 1
                    break

                # If the archiving mode is active, do not process
                # the last archiving interval, as it may be incomplete. Instead
                # exit and handle this interval with the next scheduled call
                # to the archiver. This step just requires no action here.
                if mode == 'active':
                    print('Done with ' + ch + '.')
                    xml_counts = xml_counts + 1
                    break
# ------------------------------------------------------------------------------
