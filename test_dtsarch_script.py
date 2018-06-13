import dtsarch

# To test the archiver just 'python test_dtsarch_script' and all should be cool.

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

dtsarch.archiver(sourcePath, targetPath, channels, mode='archiving',  externalBackUp=False)