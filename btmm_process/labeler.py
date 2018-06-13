import os
import yaml

def labelLocation(ds, location):
    # Assign location tags
    ds.coords['location'] = (('LAF'), [None] * ds.LAF.size)
    ds.attrs['locations'] = ';'.join(list(location.keys()))
    for l in location:
        ds.coords['location'].loc[(ds.LAF > location[l][0]) & (ds.LAF < location[l][-1])] = l

    return(ds)

def yamlDict(yamlPath):
    if not os.path.exists(yamlPath):
        raise IOError('Could not find yml file at ' + yamlPath)

    with open(yamlPath, 'r') as ymlfile:
        locationLabels = yaml.load(ymlfile)
    return locationLabels
