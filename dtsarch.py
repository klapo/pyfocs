import os
import subprocess
import datetime
import sys
import glob
import tarfile

# Script to tar and gzip the Ultima data

# ------------------------------------------------------------------------------
# File parameters including paths. Use single quotes as python expects strings.
# ------------------------------------------------------------------------------
# Mode to run the script in
mode = 'archiving'
# Path to search for local data.
sourcePath = '/Users/karllapo/Desktop/software/DTS_archive/'
# Path to where you want the data to end up
targetPath = '/Users/karllapo/Desktop/archive/'
# Names of the channels. The script will search for targetPath/channel_X.
channel_1 = 'channel 1'
channel_2 = 'channel 2'
channel_3 = 'channel 3'
channel_4 = 'channel 4'
channels = [channel_1, channel_2, channel_3, channel_4]

# Root variables for the unison command. rootLocal is where to find the
# local archive. rootBackup is where to find the backup drive. rootMobile is
# where to find removable disk. Both back up locations may not be needed. In
# that case leave the unneeded path as an empty string, e.g., ''. Logfiles are
# where you want the log of the synchronisation to be written.
rootLocal = targetPath
rootMobile = ''
logfileMobile = rootMobile + '_logfile.txt'

# ------------------------------------------------------------------------------
# Sub function -- actually zips and archives data
# ------------------------------------------------------------------------------
# Python libraries for controlling the tarballing
def make_tarfile(tarName, filesToZip):
    # Open the tarball and prepare it for writing
    with tarfile.open(tarName, "w:gz") as tar:
        # If a wildcard was passed in filesToZip...
        if '*' in filesToZip:
            # Glob the wildcard expressions together
            filesToZip = glob.glob(filesToZip)
            # Add each file to the tarball
            for f in filesToZip:
                tar.add(f, arcname=os.path.basename(f))
        # Else just tarball the whole thing
        else:
            tar.add(filesToZip, arcname=os.path.basename(filesToZip))

def archiveTool(outFile, sourceFile):
    '''
    This helper function creates the command to tar.gz the DTS XML files.
    INPUT:
        outFile = Name of the tar.gz file to create
        sourceFile = Name of the source files to compress. It MUST terminate with a
            wildcare (*) (I think).
    OUTPUT:
        flag = True if the archive was sucessfully created. False if an error was
            detected.
    '''
    # Check that the directory actually exists.
    if os.path.isdir(sourcePath):
        print('Source files: ' + sourceFile)
        print('Archiving to: ' + outFile)

        # Compress xml files
        make_tarfile(outFile, sourceFile)
        # Previous version of the code checked for errors in the archiving
        # and returned a boolean indicating if the process occurred successfully.
        return(True)

    # No files were found, exit and notify.
    else:
        print('Could not find files in specified paths. Please check sourcePath')
        return(False)

# ------------------------------------------------------------------------------
# Archive data
# ------------------------------------------------------------------------------
for ch in channels:
    channelPath = os.path.join(sourcePath, ch)

    ########
    # Active: Meant to be run with a cron job and uses the current time.
    if mode == 'active':
        now = datetime.datetime.now()
        yyyy = now.year
        mm = now.month
        dd = now.dd

        # Hours require special attention
        hh = now.hour
        if hh < 10:
            hh = '0' + str(hh)

        # Date to zip and archive
        dateFileName = '_' + yyyy + mm + dd + '-' + hh

        # Define file names using current time
        outFile = os.path.join(targetPath, channel_1 + '_' + dateFileName + '.tar.gz')
        sourceFile = os.path.join(sourcePath, channel_1, channel_1 + '_' + dateFileName + '*')

        # zip the data and move to archive
        archiveTool(outFile, sourceFile)

    ########
    # Archive: To archive previously aquired data.
    elif mode == 'archiving':
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

        # First datetimes
        t = contents[0]
        t = t.split('_')[-1]
        t = t.split('.')[0]
        year = t[0:4]
        month = t[4:6]
        day = t[6:8]
        hour = t[8:10]
        dtInit = datetime.datetime(int(year), int(month), int(day), int(hour), 0)

        # Last datetime
        t = contents[-1]
        t = t.split('_')[-1]
        t = t.split('.')[0]
        year = t[0:4]
        month = t[4:6]
        day = t[6:8]
        hour = t[8:10]
        dtFinal = datetime.datetime(int(year), int(month), int(day), int(hour), 0)

        # Span the time found in the specified directory
        dt = dtInit
        while dt <= dtFinal:
            yyyy = dt.year
            mm = dt.month
            dd = dt.day

            # Hours require special attention
            hh = dt.hour
            if hh < 10:
                hh = '0' + str(hh)
            else:
                hh = str(hh)

            # Create file names for this hour
            dateFileName = '_' + str(yyyy) + str(mm) + str(dd) + '-' + hh
            outFile = os.path.join(targetPath, ch + '_' + str(yyyy) +
                str(mm) + str(dd) + '-' + hh + '.tar.gz')
            sourceFile = os.path.join(sourcePath, ch,
                channel_1 + '_' + str(yyyy) + str(mm) + str(dd) + hh + '*')
            # Zip and archive this time period
            archiveTool(outFile, sourceFile)

            # Iterate the time
            dt = dt + datetime.timedelta(hours=1)

    print('Done with ' + ch + '. Backup files in: ' + targetPath)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Sync to the mobile backup drive.
if os.path.isdir(rootMobile):
    print('Backing up archives to mobile drive')
    print('Syncing ' + rootLocal + ' to ' + rootMobile)

    # Check the os. Unix = rsync; Windos = Unison.
    if os.name == 'nt':
        mergeCommand = [unisonPath, rootLocal, rootMobile,
                        '-logfile ' + logfileMobile, ' -force ' + rootLocal,
                        ' -batch -nodeletion ' + rootMobile]
        subprocess.check_output(mergeCommand)
        try:
            p = subprocess.check_output(mergeCommand)
        # The tar command indicated an error.
        except subprocess.CalledProcessError:
            print('Warning: syncing to the mobile backup failed.')

    elif os.name == 'posix':
        # Add a trailing backslash to make rsync behave as expected.
        if not rootMobile[-1] == '/':
            rootMobile = rootMobile + '/'
        if not rootLocal[-1] == '/':
            rootLocal = rootLocal + '/'
        mergeCommand = ['rsync', '-az', rootLocal, rootMobile]
        try:
            p = subprocess.check_output(mergeCommand)
        # The tar command indicated an error.
        except subprocess.CalledProcessError:
            print('Warning: syncing to the mobile backup failed.')

else:
    print('Warning: Mobile back-up was not found in the specified path.')
# ------------------------------------------------------------------------------
