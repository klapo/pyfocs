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
                'dt_start': pd.to_datetime(doc['startDateTimeIndex'], infer_datetime_format=True),
                'dt_end': pd.to_datetime(doc['endDateTimeIndex'], infer_datetime_format=True),
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

    for dnum, dlist in enumerate(data):
        LAF[dnum], Ps[dnum], Pas[dnum], temp[dnum] = list(map(float, dlist.split(',')))
    actualData = pd.DataFrame.from_dict({'LAF': LAF, 'Ps': Ps, 'Pas': Pas, 'temp': temp}).set_index('LAF')
    return(actualData, metaData)


def tar_read(dirData, filePrefix, fileSuffix='',
             channelName='channel 1', chunkSize=1000):
    '''
    Reads all xml files in the provided directory and turns them into netcdfs.
    '''

    prevNumChunk = 0

    # List of files to iterate over
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
            print("\r", someDumbFiles + 'File ' + str(nDumb) + ' of ' + str(nTotal), end="")

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
                ds = labelLocation(ds, location)
                ds.to_netcdf(filePrefix + '_' + str(numChunk) + '_'
                             + fileSuffix  + '.nc', 'w')
                ds.close()
                ds = None
                os.chdir(dirData)
        print('')
        # Remove the extracted xml files
        subprocess.Popen(['rm'] + glob.glob('*.xml'))
        # Preserve the chunk count across tar files
        prevNumChunk = numChunk + 1
