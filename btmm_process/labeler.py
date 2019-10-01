import os
import yaml
import numpy as np
import xarray as xr
import pandas as pd
from scipy import stats


# ------------------------------------------------------------------------------
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

    # If no location is passed, return the file to the user and notify.
    if location is None:
        return ds

    # Pre-alloate the new coordinates that are to be assigned.
    ds.coords['loc_general'] = (('LAF'), [None] * ds.LAF.size)
    ds.coords['location_flip'] = (('LAF'), [False] * ds.LAF.size)
    ds.attrs['loc_general'] = ';'.join(list(location.keys()))

    # Loop over all labels and find where they exist in the LAF domain.
    for lc in location:
        shape = np.shape(location[lc])
        # For non-continuous label locations, loop through each labeled section
        if np.size(shape) > 1:
            for loc_num in np.arange(0, max(shape)):
                LAF1 = min(location[lc][loc_num])
                LAF2 = max(location[lc][loc_num])
                ds.coords['loc_general'].loc[(ds.LAF > LAF1) & (ds.LAF < LAF2)] = lc

                # For locations where the relative start is at a larger LAF
                # than the relative end of the section indicate that we need
                # to flip the LAF
                if not location[lc][loc_num][0] > location[lc][loc_num][-1]:
                    ds.coords['location_flip'].loc[(ds.LAF > LAF1)
                                                   & (ds.LAF < LAF2)] = True

        # The label locations occur only once in the LAF domain.
        elif np.size(shape) == 1:
            if np.size(location[lc]) > 1:
                LAF1 = min(location[lc])
                LAF2 = max(location[lc])
                ds.coords['loc_general'].loc[(ds.LAF > LAF1) & (ds.LAF < LAF2)] = lc

                # For locations where the relative start is at a larger LAF than
                # the relative end of the section, indicate that we need to
                # flip the LAF.
                if not location[lc][0] > location[lc][-1]:
                    ds.coords['location_flip'].loc[(ds.LAF > LAF1)
                                                   & (ds.LAF < LAF2)] = True

        # It is a single item element (i.e., a point) to label. Find the
        # nearest point to label.
        else:
            LAF_single_loc = ds.sel(LAF=location[lc], method='nearest').LAF
            ds.coords['loc_general'].loc[(ds.LAF == LAF_single_loc)] = lc

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

    if location is None:
        return ds

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
        elif np.size(shape) == 1:
            if np.size(location[lc]) > 1:
                LAF1 = min(location[lc])
                LAF2 = max(location[lc])
                ds.coords[loc_type].loc[(ds.LAF > LAF1) & (ds.LAF < LAF2)] = lc

        # It is a single item element (i.e., a point) to label. Find the
        # nearest point to label.
        else:
            LAF_single_loc = ds.sel(LAF=location[lc], method='nearest').LAF
            ds.coords[loc_type].loc[(ds.LAF == LAF_single_loc)] = lc

    return ds


# ------------------------------------------------------------------------------
def dtsPhysicalCoords(ds, location, loc_field='loc_general',
                      coord_opt='relative', align='right'):
    '''
    Assign a physical coordinate to the xarray Dataset containing DTS data.

    Input:
        ds            -  xarray dataset with DTS data. Expects to find a
                         dimension labeled 'LAF'
        location_list -  List containing location labels for converting LAF
                         to a physical coordinate
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
            section.coords[loc_field] = (('section_num'),
                                         np.unique(section[loc_field].values))
            all_sections[section_num] = section

        # Concatenate into a single xarray Dataset
        ds_out = xr.concat(all_sections, 'section_num')
        ds_out = ds_out.swap_dims({'section_num': loc_field})

        return ds_out


# ------------------------------------------------------------------------------
def dtsPhysicalCoords_3d(ds, location):
    '''
    Assign 3D physical coordinates to the xarray Dataset containing DTS data
    converting the 1d LAF dimension into a 3D location.

    Input:
        ds            -  xarray dataset with DTS data. Expects to find a
                         dimension labeled 'LAF'
        location_list -  Dictionary containing location labels for converting
                         LAF to a physical coordinate. Must contain the fields
                         'x_coord', 'y_coord', 'z_coord', each containing a
                         list or numpy array with two numbers, and 'LAF' with
                         two elements specifying the LAF value corresponding to
                         the two elements in each coord field.
    Output:
        ds            -  xarray Dataset formatted to include just the
                         physically labeled coordinate system.
    '''

    all_sections = []

    # Extract out just the relative distance section.
    # Assumed that the first LAF section refers to the start of the
    # relative section and the last one is the end of the relative
    # distance section.
    for l in location:

        # Get the xyz coordinates for this section.
        x = location[l]['x_coord']
        y = location[l]['y_coord']
        z = location[l]['z_coord']

        # Determine the LAF for this section.
        LAF = location[l]['LAF']
        LAF1 = min(LAF)
        LAF2 = max(LAF)
        dLAF, _ = stats.mode(np.diff(ds.LAF.values))

        # Extract out just the section in question
        section = ds.loc[dict(LAF=(ds.LAF > LAF1) & (ds.LAF < LAF2))]
        num_LAF = np.size(section.LAF.values)

        # Interpolate each coordinate into a line
        (x_int, dx) = np.linspace(x[0], x[1], num=num_LAF, retstep=True)
        (y_int, dy) = np.linspace(y[0], y[1], num=num_LAF, retstep=True)
        (z_int, dz) = np.linspace(z[0], z[1], num=num_LAF, retstep=True)

        # Crude check for the mapping's consistency.
        d = (dx**2 + dy**2 + dz**2)**(0.5)
        if np.abs(d - dLAF) > 0.5 * dLAF:
            delta_str = ('\ndx = ' + str(dx)
                         + '\n' + 'dy = ' + str(dy)
                         + '\n' + 'dz = ' + str(dz)
                         + '\n' + 'total = ' + str(dz)
                         + '\n' + 'dLAF = ' + str(dLAF))
            raise ValueError('Mapping problem detected for '
                             + location[l]['long name']
                             + '. Inferred data spacing is '
                             + delta_str)

        # Determine the orientation of the fiber. If
        # LAF decreases with distance along the section
        # flip the array.
        if not LAF1 > LAF2:
            section['LAF'] = np.flip(section.LAF.values, 0)

        # Assign the physical coordinates using a pandas multiindex
        midx = pd.MultiIndex.from_arrays([x_int, y_int, z_int],
                                         names=('x', 'y', 'z'))
        section.coords['xyz'] = ('LAF', midx)
        section = section.swap_dims({'LAF': 'xyz'})
        section = section.drop('LAF')
        all_sections.append(section)

    # Concatenate along the physical coordinate MultiIndex
    ds_out = xr.concat(all_sections, 'xyz')
    # xarray does not yet support writing a MultiIndex to netcdf format.
    ds_out = ds_out.reset_index('xyz')

    # Here is how to recreate the MultiIndex after resetting it like above.
    # midx = pd.MultiIndex.from_arrays([ds_out.x, ds_out.y, ds_out.z], names=('x', 'y', 'z'))
    # ds_out = ds_out.drop(['x', 'y', 'z'])
    # ds_out = ds_out.assign_coords(xyz = midx)
    return ds_out


# ------------------------------------------------------------------------------
def create_multiindex(ds, coords=['x', 'y', 'z']):
    '''
    xarray does not (yet?) support writing a MultiIndex to netcdf format. To
    get around this unsupported behavior, this function recreates a MultiIndex
    using the saved 'x', 'y', and 'z' coordinates. Supports the ability to use
    arbitrary coordinates with the default matching the behavior expected by
    PyFOX.
    '''
    # Recreate the multiindex for an xarray dataset saved as a netcdf
    coord_list = [ds[c] for c in coords]

    # Create a pandas multiindex
    midx = pd.MultiIndex.from_arrays(coord_list, names=coords)

    # Drop the old, uncombined coordinates since they conflic with assigning
    # the MultiIndex
    ds = ds.drop(coords)
    print(ds)
    # Assing the MultiIndex, this recreates the individual coordinates that
    # were just dropped.
    ds = ds.assign_coords(xyz = midx)

    return ds


# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
def yamlWrite(yamlPath):
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
