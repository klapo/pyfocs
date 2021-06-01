import numpy as np
import xarray as xr


def cc_step_detection(da, dim, window_size, resolution, add_padding=None):
    '''Author: Bart Schilperoort

    Computes the cross-correlation between a step function and the input
    DataArray, over the specified dimension.

    Returns an xarray DataArray with the coordinate of the maximum
    cross-correlation.
    '''

    if add_padding:
        # Bottom data:
        bottom_data = da.isel(
            {dim: slice(None, add_padding)}).mean(dim=dim).values

        bottom_height = np.arange(-window_size, 0)*resolution +\
            da.isel({dim: 0})[dim].values

        da_bottom = da.interp(
            {dim: bottom_height}, kwargs={'fill_value': 0}
        )
        da_bottom = da_bottom + np.vstack(bottom_data)

        # Top data:
        top_data = da.isel(
            {dim: slice(-add_padding, None)}).mean(dim=dim).values

        top_height = np.arange(0, window_size)*resolution +\
            da.isel({dim: -1})[dim].values + resolution

        da_top = da.interp(
            {dim: top_height}, kwargs={'fill_value': 0}
        )
        da_top = da_top + np.vstack(top_data)

        da = xr.concat([da_bottom, da, da_top], dim=dim)

    step = np.hstack((-np.ones(window_size//2),
                      np.ones(window_size//2)))

    da_step = xr.DataArray(step, dims=['window'])

    da = da.rolling(
        {dim: len(step)},
        center=True
        ).construct('window').dot(da_step)

    return da[dim][da.argmax(dim=dim)]
