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
