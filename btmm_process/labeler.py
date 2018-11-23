import os
import yaml
import numpy as np
import xarray as xr


def labelLoc_general(ds, location):
    '''
    Assign location tags to an xarray Dataset containing DTS data.

    Input:
        ds       -  xarray dataset with DTS data. Expects to find a dimension
                    labeled 'LAF'
        location -  dictionary specifying location labels of the form:
                    dict['label'] = [float of section LAF start, float of
                    section LAF end]
                    If the section has a fiber that 'points'
    Output:
        ds       -  the same xarray that was passed to the function, but with
                    the 'location' and 'location_flip' coordinates.
    '''

    # Pre-alloate the new coordinates that are to be assigned.
    ds.coords['location'] = (('LAF'), [None] * ds.LAF.size)
    ds.coords['location_flip'] = (('LAF'), [False] * ds.LAF.size)
    ds.attrs['locations'] = ';'.join(list(location.keys()))

    # Loop over all labels and find where they exist in the LAF domain.
    for lc in location:
        shape = np.shape(location[lc])
        # For non-continuous label locations, loop through each labeled section
        if np.size(shape) > 1:
            for loc_num in np.arange(0, max(shape)):
                LAF1 = min(location[lc][loc_num])
                LAF2 = max(location[lc][loc_num])
                ds.coords['location'].loc[(ds.LAF > LAF1) & (ds.LAF < LAF2)] = lc

                # For locations where the relative start is at a larger LAF
                # than the relative end of the section indicate that we need
                # to flip the LAF
                if not location[lc][loc_num][0] > location[lc][loc_num][-1]:
                    ds.coords['location_flip'].loc[(ds.LAF > LAF1)
                                                   & (ds.LAF < LAF2)] = True

        # The label locations occur only once in the LAF domain.
        else:
            if np.size(location[lc]) > 1:
                LAF1 = min(location[lc])
                LAF2 = max(location[lc])
                ds.coords['location'].loc[(ds.LAF > LAF1) & (ds.LAF < LAF2)] = lc
    
                # For locations where the relative start is at a larger LAF than
                # the relative end of the section, indicate that we need to
                # flip the LAF.
                if not location[lc][0] > location[lc][-1]:
                    ds.coords['location_flip'].loc[(ds.LAF > LAF1)
                                                   & (ds.LAF < LAF2)] = True
            else:
                ds.coords['location'].loc[(ds.LAF == location[lc])] = lc 
                
    return ds

def labelLoc_additional(ds, location, loc_type):
    '''
    Assign location tags to an xarray Dataset containing DTS data.

    Input:
        ds       -  xarray dataset with DTS data. Expects to find a dimension
                    labeled 'LAF'
        location -  dictionary specifying location labels of the form:
                    dict['label'] = [float of section LAF start, float of
                    section LAF end]
                    If the section has a fiber that 'points'
    Output:
        ds       -  the same xarray that was passed to the function, but with
                    the new 'loc_type'.
    '''

    # Pre-alloate the new coordinates that are to be assigned.
    ds.coords[loc_type] = (('LAF'), [None] * ds.LAF.size)
    ds.attrs[loc_type] = ';'.join(list(location.keys()))

    # Loop over all labels and find where they exist in the LAF domain.
    for lc in location:
        shape = np.shape(location[lc])

        # For non-continuous label locations, loop through each labeled section
        if np.size(shape) > 1:
            for loc_num in np.arange(0, max(shape)):
                LAF1 = min(location[lc][loc_num])
                LAF2 = max(location[lc][loc_num])
                ds.coords[loc_type].loc[(ds.LAF > LAF1) & (ds.LAF < LAF2)] = lc

        # The label locations occur only once in the LAF domain.
        else:
            if np.size(location[lc]) > 1:
                LAF1 = min(location[lc])
                LAF2 = max(location[lc])
                ds.coords[loc_type].loc[(ds.LAF > LAF1) & (ds.LAF < LAF2)] = lc
            else:
                ds.coords[loc_type].loc[(ds.LAF == location[lc])] = lc 
                
    return ds

def dtsPhysicalCoords(ds, location, coord_opt='relative', align='right'):
    '''
    Assign a physical coordinate to the xarray Dataset containing DTS data.

    Input:
        ds            -  xarray dataset with DTS data. Expects to find a
                         dimension labeled 'LAF'
        location_list -  List containing location labels for convertitn LAF
                         to a physically coordinate
    Output:
        ds            -  xarray Dataset formatted to include just the
                         physically labeled coordinate system.
    '''

    all_sections = []

    if coord_opt == 'relative':
        x_max = 0

        ########
        # Extract out just the relative distance section.
        # Assumed that the first LAF section refers to the start of the
        # relative section and the last one is the end of the relative
        # distance section.
        for l in location:
            shape = np.shape(location[l])

            ########
            # For non-contiguous label locations, loop through each labeled
            # section.
            if np.size(shape) > 1:
                for loc_num in np.arange(0, max(shape)):
                    LAF1 = min(location[l][loc_num])
                    LAF2 = max(location[l][loc_num])

                    # Extract out just the section
                    section = ds.loc[dict(LAF=(ds.LAF > LAF1)
                                          & (ds.LAF < LAF2))]

                    # For locations where the relative start is at a larger
                    # LAF than the relative end of the section, indicate that
                    # we need to flip the LAF
                    if not location[l][loc_num][0] > location[l][loc_num][-1]:
                        section['LAF'] = np.flip(section.LAF.values, 0)
                        section.coords['x'] = section.LAF - section.LAF[-1]

                    # LAF increases with test section distance
                    else:
                        section.coords['x'] = section.LAF - section.LAF[0]

                    # Keep only the longest x dimension.
                    if section['x'].max() > x_max:
                        x_max = section['x'].max()

                    # Assign to a dictionary for further processing
                    all_sections.append(section)

            ########
            # The labeled location occurs only once in the LAF domain.
            else:
                LAF1 = min(location[l])
                LAF2 = max(location[l])

                # Extract out just the section
                section = ds.loc[dict(LAF=(ds.LAF > LAF1) & (ds.LAF < LAF2))]

                # For locations where the relative start is at a larger LAF
                # than the relative end of the section, indicate that we need
                # to flip the LAF.
                if not location[l][0] > location[l][-1]:
                    section['LAF'] = np.flip(section.LAF.values, 0)
                    section.coords['x'] = section.LAF - section.LAF[-1]

                # LAF increases with test section distance
                else:
                    section.coords['x'] = section.LAF - section.LAF[0]

                # Keep only the longest x dimension.
                if section['x'].max() > x_max:
                    x_max = section['x'].max()
                    x = section['x']

                # Assign to a dictionary for further processing
                all_sections.append(section)

        for section_num, section in enumerate(all_sections):

            # Find the difference in the length of each section
            delta_x = x_max - section['x'].max()

            # Should the distances line up on the left or right?
            if 'right' in align:
                section.coords['x'] = section['x'] + delta_x

            if 'left' in align:
                section.coords['x'] = section['x'] - delta_x

            # Reindex to the longest x distance
            section = section.swap_dims({'LAF': 'x'})
            section = section.reindex(x=x, method='nearest')

            # Create the dimension we will concatenate along (don't worry,
            # section_num will disappear after swapping dimensions)
            section.coords['location'] = (('section_num'),
                                          np.unique(section.location.values))
            all_sections[section_num] = section

        # Concatenate into a single xarray Dataset
        ds_out = xr.concat(all_sections, 'section_num')
        ds_out = ds_out.swap_dims({'section_num': 'location'})

        return ds_out


def yamlDict(yamlPath):
    '''
    Reads the .yml label location file and returns the location dictionary
    for labeling the DTS array.
    INPUT:
        yamlPath - path to the .yml file to read
    OUTPUT:
        locationLabels - python dictionary of the locations/LAF pairs.
    '''
    if not os.path.exists(yamlPath):
        raise IOError('Could not find yml file at ' + yamlPath)

    with open(yamlPath, 'r') as ymlfile:
        locationLabels = yaml.load(ymlfile)
    return locationLabels
