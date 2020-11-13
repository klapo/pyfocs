import os
import xmltodict
import subprocess
import pandas as pd
import xarray as xr
import glob
import tarfile
import numpy as np


# Error classes
class Error(Exception):
    """Base class for exceptions in this module."""
    pass


class CorruptedXMLError(Error):
    """Exception raised for errors in the input.

    Attributes:
    """


def xml_read(dumbXMLFile):
    '''
    Opens the given xml file and reads the dts data contained within.
    '''
    # The Ultimas/XT return incomplete files at times, which must be discarded.
    if 'incomplete' in dumbXMLFile:
        raise CorruptedXMLError
    # Continue to try and read the file.
    try:
        with open(dumbXMLFile) as dumb:
            doc = xmltodict.parse(dumb.read())
    except Exception:
        # The exception in the xmltodict code is poorly formatted. So we do a
        # general catch here and hope for the best.
        # Raising this error allows us to catch corrupted files.
        raise CorruptedXMLError
    # Remove all of the bullshit
    doc = doc['logs']['log']

    # Extract units/metadata info out of xml dictionary
    metaData = {'LAF_beg': float(doc['startIndex']['#text']),
                'LAF_end': float(doc['endIndex']['#text']),
                'dLAF': float(doc['stepIncrement']['#text']),
                'dt_start': pd.to_datetime(doc['startDateTimeIndex'],
                                           infer_datetime_format=True,
                                           utc=True),
                'dt_end': pd.to_datetime(doc['endDateTimeIndex'],
                                         infer_datetime_format=True,
                                         utc=True),
                'probe1Temperature': float(doc['customData']
                                           ['probe1Temperature']['#text']),
                'probe2Temperature': float(doc['customData']
                                           ['probe2Temperature']['#text']),
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
        raise IOError('Unrecognized xml format... dumping first row \n'
                      + data[0])

    # Single ended data
    if 'single_ended' in dtsType:
        for dnum, dlist in enumerate(data):
            LAF[dnum], Ps[dnum], Pas[dnum], temp[dnum] = list(map(float,
                                                              dlist.split(','))
                                                              )
        actualData = pd.DataFrame.from_dict({'LAF': LAF,
                                             'Ps': Ps,
                                             'Pas': Pas,
                                             'temp': temp}).set_index('LAF')

    # Double ended data
    elif 'double_ended' in dtsType:
        rPs = np.empty_like(LAF)
        rPas = np.empty_like(LAF)

        for dnum, dlist in enumerate(data):
            LAF[dnum], Ps[dnum], Pas[dnum], rPs[dnum], rPas[dnum], temp[dnum], = list(map(float, dlist.split(',')))

        actualData = pd.DataFrame.from_dict({'LAF': LAF,
                                             'Ps': Ps,
                                             'Pas': Pas,
                                             'rPs': rPs,
                                             'rPas': rPas,
                                             'temp': temp}).set_index('LAF')

    return(actualData, metaData)


def archive_read(cfg, write_mode='preserve', prevNumChunk=0):
    '''
    Reads all archived xml files in the provided directory
    and turns them into netcdfs.

    Optional arguments:
    write_mode    -  Determines if the function peaks to see if the file it
                     would create exists.
                        write_mode == 'overwrite', no peaking, creates new
                                      file.
                        write_mode == 'preserve', peaks, does not process files
                                      that already exist.
    prevNumChunk  -  The chunk number to assign the output netcdf name. When
                     running across multiple experiments/directories it can be
                     useful to specify your own value.
    '''

    # Assign values
    dirDataOriginal = cfg['directories']['dirArchive']
    dirProcessed = cfg['directories']['dirRawNetcdf']
    channelNames = cfg['directories']['channelName']

    # If corrupt files are found we need to keep track of them (and skip them).
    corrupt_file_count = 0
    corrupt_file_list = []

    # Loop through each channel provided
    for chan in np.atleast_1d(channelNames):
        # Check directories
        dirData = dirDataOriginal
        if not os.path.isdir(dirData):
            raise IOError('Data directory was not found at ' + dirData)
        os.chdir(dirData)

        # List of files to iterate over
        dirConTar = [dC for dC in os.listdir() if chan in dC
                     and '.tar.gz' in dC]
        dirConTar.sort()

        # Untar files
        for tFile in dirConTar:
            print(tFile)

            # Name of the resulting netcdf
            nc_out_name = 'raw_' + tFile.split('.')[0] + '.nc'

            # Skip this archive if the netCDF already exists and we are
            # not overwriting.
            if write_mode == 'preserve':
                if os.path.isfile(os.path.join(dirProcessed, nc_out_name)):
                    print('... exists. No overwriting.')
                    continue

            # Extract the archive
            t = tarfile.open(tFile)
            t.extractall()
            t.close

            # List of files to iterate over
            dirConXML = [dC for dC in os.listdir() if chan in dC
                         and '.xml' in dC]
            dirConXML.sort()
            nTotal = np.size(dirConXML)
            ds = None
            ds_list = []

            # Read each xml file, assign to an xarray Dataset, concatenate
            # along the time dimension, and output data with a given chunk size
            # to netcdf format.
            for nDumb, someDumbFiles in enumerate(dirConXML):
                if '.xml' not in someDumbFiles:
                    continue
                print("\r", someDumbFiles + ' File ' + str(nDumb + 1) + ' of '
                      + str(nTotal), end="")

                # Read the file
                try:
                    df, meta = xml_read(someDumbFiles)
                except CorruptedXMLError:
                    corrupt_file_count = corrupt_file_count + 1
                    corrupt_file_list.append(someDumbFiles)
                    continue

                # Create a temporary xarray Dataset
                temp_Dataset = xr.Dataset.from_dataframe(df)
                temp_Dataset.coords['time'] = meta['dt_start']

                # Assing the reference probes to the dataset. This is no
                # longer an option.
                temp_Dataset['probe1Temperature'] = meta['probe1Temperature']
                temp_Dataset['probe2Temperature'] = meta['probe2Temperature']

                # If the flag is 'external' than no probe temperature field is
                # returned as the external data stream must be handled
                # separately.

                # Create a list of xarray Datasets
                ds_list.append(temp_Dataset)
            print('\n Concatenating netcdfs within archive...')
            try:
                ds = xr.concat(ds_list, dim='time')
            except ValueError:
                raise ValueError('No xml files were found within: ' + tFile)

            # Create a raw netcdf file for each archive interval. This means
            # that the archive interval dicates the speed/efficiency of the
            # later calibration step.
            ds.attrs = {'LAF_beg': meta['LAF_beg'],
                        'LAF_end': meta['LAF_end'],
                        'dLAF': meta['dLAF']}

            # Delete this block once the new calibration routine is working.
            # # Label the Ultima PT100 data. These names are used in
            # # calibration and must match the 'refField' variables.
            # try:
            #     ds = ds.rename({'probe1Temperature': cfg['probe1'],
            #                     'probe2Temperature': cfg['probe2']})
            # # Passing no values for the probe names causes a ValueError
            # # since both will be `None`. In that case, skip over the
            # # probe naming.
            # except ValueError:
            #     print('No PT100 field names were passed.')
            # # No new names are handed to the probes.
            # except KeyError:
            #     print('No PT100 field names were passed.')


            # Save to netcdf
            os.chdir(dirProcessed)
            ds.to_netcdf(nc_out_name, 'w')

            # Close the netcdf and release the memory.
            ds.close()
            ds = None

            os.chdir(dirData)
            print('')
            # Remove the extracted xml files. The second wildcard character
            # catches the 'incomplete' xml files that occur with power outages.
            files_to_remove = glob.glob('*.xml*')
            for file in files_to_remove:
                os.remove(file)

    # Notify the user if corrupt data are found.
    if corrupt_file_count > 0:
        print('Corrupt files: ' + str(corrupt_file_count))
        print(*corrupt_file_list, sep='\n')
