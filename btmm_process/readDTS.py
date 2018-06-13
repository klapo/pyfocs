import os
import xmltodict
import subprocess
import pandas as pd
import xarray as xr
import glob
import tarfile
import numpy as np
from .labeler import labelLocation, yamlDict

def xml_read(dumbXMLFile):
    '''
    Opens the given xml file and reads the dts data contained within.
    '''

    with open(dumbXMLFile) as dumb:
        doc = xmltodict.parse(dumb.read())

    # Remove all of the bullshit
    doc = doc['logs']['log']

    # Extract units/metadata info out of xml dictionary
    metaData = {'LAF_beg': float(doc['startIndex']['#text']),
                'LAF_end': float(doc['endIndex']['#text']),
                'dLAF': float(doc['stepIncrement']['#text']),
                'dt_start': pd.to_datetime(doc['startDateTimeIndex'],
                                           infer_datetime_format=True),
                'dt_end': pd.to_datetime(doc['endDateTimeIndex'],
                                         infer_datetime_format=True),
                'probe1Temperature': float(doc['customData']['probe1Temperature']['#text']),
                'probe2Temperature': float(doc['customData']['probe2Temperature']['#text']),
                'fiberOK': int(doc['customData']['fibreStatusOk']),
               }

    # Extract data
    data = doc['logData']['data']

    numEntries = np.size(data)
    LAF = np.empty(numEntries)
    Ps = np.empty_like(LAF)
    Pas = np.empty_like(LAF)
    temp = np.empty_like(LAF)

    # Check dts type based on the number of columns
    if len(data[0].split(',')) == 4:
        dtsType = 'single_ended'
    elif len(data[0].split(',')) == 6:
        dtsType = 'double_ended'
    else:
        raise IOError('Unrecognized xml format... dumping first row \n' + data[0])

    # Single ended data
    if 'single_ended' in dtsType:
        for dnum, dlist in enumerate(data):
            LAF[dnum], Ps[dnum], Pas[dnum], temp[dnum] = list(map(float,
                                                              dlist.split(',')))
        actualData = pd.DataFrame.from_dict({'LAF': LAF, 'Ps': Ps, 'Pas': Pas, 'temp': temp}).set_index('LAF')

    # Double ended data
    elif 'double_ended' in dtsType:
        rPs = np.empty_like(LAF)
        rPas = np.empty_like(LAF)

        for dnum, dlist in enumerate(data):
            LAF[dnum], Ps[dnum], Pas[dnum], rPs[dnum], rPas[dnum], temp[dnum], = list(map(float, dlist.split(',')))

        actualData = pd.DataFrame.from_dict({'LAF': LAF, 'Ps': Ps,
                                             'Pas': Pas, 'rPs': rPs,
                                             'rPas': rPas, 'temp': temp}).set_index('LAF')

    return(actualData, metaData)

def tar_read(cfg):
    '''
    Reads all xml files in the provided directory and turns them into netcdfs.
    '''
    # Assign values
    dirData = cfg['directories']['dirData']
    dirProcessed = cfg['directories']['dirProcessed']
    channelName = cfg['directories']['channelName']
    labelsFile = cfg['directories']['labelsFilePath']
    filePrefix = cfg['fileName']['filePrefix']
    fileSuffix = cfg['fileName']['fileSuffix']
    chunkSize = cfg['dataProperties']['chunkSize']

    # Read label configuration files
    labels = yamlDict(labelsFile)

    # Start keeping track of chunks
    prevNumChunk = 0

    # List of files to iterate over

    # Check directories
    if not os.path.isdir(dirData):
        raise IOError('Data directory was not found at ' + dirData)
    os.chdir(dirData)
    dirConTar = [dC for dC in os.listdir() if channelName in dC and '.tar.gz' in dC]
    dirConTar.sort()

    # Untar files
    for tFile in dirConTar:
        print(tFile)
        t = tarfile.open(tFile)
        t.extractall()
        t.close

        dirCon = [dC for dC in os.listdir() if channelName in dC and '.xml' in dC]
        dirCon.sort()
        nTotal = np.size(dirCon)
        ds = None

        # Read each xml file, assign to an xarray Dataset, concatenate along
        # the time dimension, and output data with a given chunk size to netcdf
        # format.
        for nDumb, someDumbFiles in enumerate(dirCon):
            if not '.xml' in someDumbFiles:
                continue
            print("\r", someDumbFiles + ' File ' + str(nDumb) + ' of ' + str(nTotal), end="")

            # Read the file
            df, meta = xml_read(someDumbFiles)

            # Create a temporary xarray Dataset
            temp_Dataset = xr.Dataset.from_dataframe(df)
            temp_Dataset.coords['time'] = meta['dt_start']
            temp_Dataset['probe1Temperature'] = meta['probe1Temperature']
            temp_Dataset['probe2Temperature'] = meta['probe2Temperature']
            temp_Dataset['fiberStatus'] = meta['fiberOK']

            if ds:
                ds = xr.concat([ds, temp_Dataset], dim='time')
            else:
                ds = temp_Dataset

            # Chunking/saving to avoid memory errors
            if np.mod(nDumb + 1, chunkSize) == 0 or nDumb == nTotal - 1:
                os.chdir(dirProcessed)
                numChunk = np.floor_divide(nDumb, chunkSize) + prevNumChunk
                ds.attrs = {'LAF_beg': meta['LAF_beg'],
                            'LAF_end': meta['LAF_end'],
                            'dLAF': meta['dLAF']}
                ds = labelLocation(ds, labels)
                ds.rename({'probe1Temperature': cfg['dataProperties']['probe1Temperature'],
                           'probe2Temperature': cfg['dataProperties']['probe2Temperature']},
                           inplace=True)
                ds.to_netcdf(filePrefix + '_raw' + str(numChunk) + '_'
                             + fileSuffix  + '.nc', 'w')
                ds.close()
                ds = None
                os.chdir(dirData)
        print('')
        # Remove the extracted xml files
        subprocess.Popen(['rm'] + glob.glob('*.xml'))
        # Preserve the chunk count across tar files
        prevNumChunk = numChunk + 1
