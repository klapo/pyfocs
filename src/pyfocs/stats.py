import numpy as np
import xarray as xr


def noisymoments(x, maxlag=10):
    '''
    Function for estimating statistical moments from noisy time series
    Utilises the methods originally proposed in
    Lenschow, D. H., V. Wulfmeyer, and C. Senff (2000),
    Measuring second- through fourth-order moments in noisy data,
    J. Atmos. Oceanic Technol., 17, 1330� 1347.
    https://doi.org/10.1175/1520-0426(2000)017<1330:MSTFOM>2.0.CO;2
    INPUT:
     x = Numpy 1d array, should contain data with white noise (e.g.,
     observations from a time series)
    OUTPUT:
     var = variance estimates, different variables at different rows
     Sk = skewness estimates, different variables at different rows
     noisevar = noise variance estimates, different variables at different rows

     Adapted from code written by O. Peltola to python, Feb 2020
     '''

    lgs = np.arange(-maxlag, maxlag)
    # Select the lag indices for 1, 2, 3, 4
    f = slice(maxlag + 1, maxlag + 5)

    xp = x - np.nanmean(x)

    # variance
    c = norm_xcorr(xp, xp, lag=maxlag, scaleopt='biased')

    M = np.tile(lgs[f], [2, 1]).T ** (0, 1)

    # Linear algebra solution to a x = b for x
    p, _, _, _ = np.linalg.lstsq(M, c[f], rcond=None)
    p = np.flip(p)

    # Easy peasy conversion with a bonus reduction of redundant operations
    x_pval = np.polyval(p, 0)
    noisevar = c[maxlag] - x_pval
    var = x_pval
    if var < 0:
        var = np.nan

    # skewness
    # Unnormalized autocorrelation
    c = norm_xcorr(xp, xp * xp, lag=maxlag, scaleopt='biased')
    p, _, _, _ = np.linalg.lstsq(M, c[f], rcond=None)
    p = np.flip(p)
    Sk = np.polyval(p, 0) / var**(3 / 2)

    return var, Sk, noisevar


def norm_xcorr(x1, x2, lag=None, remove_mean=False, scaleopt='none'):
    '''
    Normalized cross correlations
    '''
    if remove_mean or scaleopt=='coef':
        x1 = x1 - np.nanmean(x1)
        x2 = x2 - np.nanmean(x2)

    if scaleopt == 'coef':
        # Returns normalized correlation coefficients.
        if not len(x1) == len(x2):
            raise ValueError('x1 and x2 must be the same length for scaleopt ' + scaleopt)
        norm_zeropad2 = np.hstack((np.arange(1, np.size(x1), 1),
                                   np.arange(np.size(x2), 0, -1))) / np.size(x2)
        full_norm = norm_zeropad2 * np.sqrt(np.sum(x1**2) * np.sum(x2**2))
    elif scaleopt == 'biased':
        # Equivalent to the matlab 'biased' option in xcorr
        if not len(x1) == len(x2):
            raise ValueError('x1 and x2 must be the same length for scaleopt ' + scaleopt)
        full_norm = len(x1)
    elif scaleopt == 'none':
        # Return the unscaled, raw cross correlations.
        full_norm = 1
    else:
        print('Unrecognized scaleopt')
        raise ValueError

    norm_xcorr = np.correlate(x1, x2, "full") / full_norm

    if not lag:
        return norm_xcorr
    elif type(lag) == int:
        norm_xcorr = norm_xcorr[norm_xcorr.size // 2 - lag: norm_xcorr.size//2 + lag + 1]
        return norm_xcorr
    elif not type(lag) == int:
        raise TypeError('lag must be an int')

    return norm_xcorr


def block_diff(da, indexer, window_size, step_size):
    '''
    Computes the finite difference derivative over an arbitrary window and dimension.
    Returns an xarray DataArray with the same coordinates as the original DataArray.
    '''

    # Rolling average
    da_roll = da.rolling({indexer: window_size}, center=True, min_periods=window_size // 2).mean()

    # Dictionaries of our indexer for this window
    sel_var_dict1 = {indexer: slice(window_size, None)}
    sel_var_dict2 = {indexer: slice(None, -window_size)}
    # And the dictionary of the window mid-point.
    sel_var_dict_mid = {indexer: slice(window_size // 2, -window_size // 2)}

    # Block difference using the numpy arrays to avoid xarray trying to automatically align
    # the two objects.
    bdiff_np = da_roll.isel(sel_var_dict1).values - da_roll.isel(sel_var_dict2).values

    # Convert to physical units
    bdiff_np = bdiff_np / (window_size * step_size)

    # Re-construct the xarray object. Here I make some assumptions just to make the example work.
    dim_names = da.coords.dims
    indexer_midpoint = da.isel(sel_var_dict_mid)[indexer]

    dim_list = [da[d] for d in dim_names if indexer not in d]
    dim_list.append(indexer_midpoint[indexer])

    bdiff_xr = xr.DataArray(bdiff_np, coords=dim_list, dims=dim_names)

    # Re-index to the original DataArray
    bdiff_xr = bdiff_xr.reindex_like(da)

    return bdiff_xr
