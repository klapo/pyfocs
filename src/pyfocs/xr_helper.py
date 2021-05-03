def swap_sel(ds, coord1, coord2, selection):
    '''
    Selects data using a non-dimension coordinate and then returns an xarray
    Dataset with the original dimensions.
    '''
    try:
        ret_ds = ds.swap_dims({coord1: coord2}).loc[{coord2: selection}].\
            swap_dims({coord2: coord1})
    except ValueError:
        ret_ds = ds.swap_dims({coord1: coord2}).loc[{coord2: selection}]

    return(ret_ds)


def xr_unique_index(ds, dim):
    '''
    It is not uncommon for DTS data to have shared xyz labels, e.g. at vertices
    between line segments, creating non-unique indices. This problem makes it
    impossible to use many xarray functions that rely on unique indices. It
    can be fixed by re-indexing using just unique values. Nominally, this
    approach throws away a point, but testing has revealed this to be
    acceptable.

    Parameters
    ----------
        ds : xarray dataset, contains dimension 'dim'
        dim : string, corresponding to dimension in ds

    Returns
    ----------
        ds : xarray with non-unique labels along 'dim' dropped
    '''


    _, index = np.unique(ds[dim], return_index=True)
    ds = ds.isel({dim: index})
    return ds
