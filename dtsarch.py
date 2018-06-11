import os
import subprocess
import datetime
import sys
import glob
import tarfile
import dirsync
import logging

# Script to tar and gzip the Ultima data

# ------------------------------------------------------------------------------
# Define paths here. Use single quotes as python expects strings.
# ------------------------------------------------------------------------------
# Mode to run the script in
mode = 'archiving'

# Path to search for local data.
sourcePath = '/Users/karllapo/Desktop/software/DTS_archive/'
# Path to where you want the data to end up
targetPath = '/Users/karllapo/Desktop/archive/'

# Paths to the source data (dirLocal) and backup drive (dirMobile)
dirLocal = targetPath
dirMobile = '/Volumes/NO NAME/test_archive'
logfileMobile = 'dtsarch_logfile.txt'

# Names of the channels. The script will search for targetPath/channel_X.
channel_1 = 'channel 1'
channel_2 = 'channel 2'
channel_3 = 'channel 3'
channel_4 = 'channel 4'
channels = [channel_1, channel_2, channel_3, channel_4]
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Sub function -- actually zips and archives data
# ------------------------------------------------------------------------------
# Python libraries for controlling the tarballing
def make_tarfile(sourceDir, tarName, filesToZip):
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
    sourceDir = os.path.dirname(filesToZip)
    
    if os.path.isdir(sourceDir):
        print('Source files: ' + sourceDir)
        print('Archiving to: ' + tarName)
    
        # Open the tarball and prepare it for writing and compress the xml files
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

        # Exit with a successful indicator
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
# ------------------------------------------------------------------------------
# Open up a file to log the syncing
os.chdir(dirMobile)
logger = logging.basicConfig(filename=logfileMobile, level='info',
                             format='%(asctime)s %(message)s',
                             datefmt='%m/%d/%Y %H:%M:%S')

if os.path.isdir(dirMobile):
    print('Backing up archives to mobile drive')
    print('Syncing ' + dirLocal + ' to ' + dirMobile)

    dirsync.sync(dirLocal, dirMobile, 'diff', logger=logger)
    dirsync.sync(dirLocal, dirMobile, 'sync', logger=logger)

else:
    print('Warning: Mobile back-up was not found in the specified path.')
# ------------------------------------------------------------------------------
