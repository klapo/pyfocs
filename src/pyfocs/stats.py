import numpy as np
import xarray as xr
from scipy.stats import f as f_dist
from scipy.stats import binned_statistic
from scipy import signal


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


def coherent_scales(
    ds1,
    ds2,
    window_middle,
    window_step,
    tcoord='time',
    scoord='LAF',
    thresh_type="sig",
    freq_sub_sample=1,
    level=5,
    **kwargs
):
    '''
    Calculates the time scale at which a specified coherence level is
    achieved between data1 (1d reference data, e.g. from sonic) and data2 (2d
    data, e.g. from DTS) at a variety of spatial averaging sizes.

    Input Parameters:
        ds1 : 1d xarray DataArray, e.g. sonic data. Must share a coordinate with
            ds2 (i.e. a time dimension, tcoord)
        ds2 : 2d xarray DataArray e.g. from DTS.
            Must share a coordinate with ds1 (i.e. a time dimension, tcoord)
        window_middle : the center of the spatial averaging window
        window_step : the size of steps over which to average spatially
        tcoord : string, name of the shared (i.e. time) coordinate in ds1 and ds2.
            Default coord name is 'time'
        scoord : string, name of the non-shared (i.e. space) coordinate in ds1 and ds2.
            Default coord name is 'LAF'
        thresh_type : Can be 'sig' or 'coh'. Decides whether time scales are determined
            based on a fixed minimum coherence or based on significance testing
            according to Biltoft and Pardyjak (2009). Defaults to 'sig'
        level : Defines either the required coherence (number between 0 and 1) or
            the required significance level for signals to be classified as coherent
            at the given time scale. Defaults to 5.
        kwargs : Optional kwargs to pass to the signal.coherence function.


    Further optional input arguments from scipy signal.coherence:
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.coherence.html

    Return:
        tau : Array of the averaging times necessary for achieving the specified
            coherence level for all spatial scales.

    '''

    # Checking input variables
    if tcoord not in ds1.dims or tcoord not in ds2.dims:
        raise ValueError('Expect 1d data for data1 (e.g. point observation time series).')

    if (ds1[tcoord] == ds2[tcoord]).sum(dim=tcoord) < len(ds1[tcoord]):
        raise ValueError('Lengths of compared time series aren`t matching.')

    if scoord not in ds2.dims:
        raise ValueError('Expected to find {} in ds2'.format(scoord))

    if window_middle < ds2[scoord].min() or window_middle > ds2[scoord].max():
        raise ValueError('Selected edges don`t exist.')

    # Determine the number of windows to search over
    num_windows = np.ceil(np.abs((ds2[scoord] - window_middle) / window_step)).max().values * 2
    # preallocate array for the resolvable time scales
    tau = np.full(num_windows.astype(int), np.nan)

    # Adding default details on the method of calculating the coherence.
    kwargs.setdefault('fs', 1) # assumes 1 Hz as default sampling frequency
    kwargs.setdefault('nperseg', 3600 * kwargs['fs']) # default length for fft: 1 hour

    # loop through different spatial averaging lengths
    for index, wnum in enumerate(np.arange(1, num_windows + 1)):
        #select locations and average spatially
        ds2_sel = ds2.sel(
            {
                scoord: slice(
                    window_middle - (window_step / 2) * wnum,
                    window_middle + (window_step / 2) * wnum
                )
            }
        ).mean(dim=scoord)
        #calculate coherence between signals
        freq, coh = signal.coherence(ds1.values, ds2_sel.values, **kwargs)

        # Spectral binning (necessary for larger nperseg when trying to
        # resolve high frequency and small values of coherence)
        if freq_sub_sample > 1:
            coh, freq, _ = binned_statistic(
                freq,
                coh,
                bins=freq[0:-1:freq_sub_sample],
                statistic='mean'
            )
            freq = freq[0:-1]

        #calculate threshold for threshold type 'significance'
        if thresh_type == 'sig':
            #number of independent spectra
            n_ind = np.size(ds1) / kwargs['nperseg']
            # Upper tail critical value of the F distribution
            F_p = f_dist.ppf(1 - level * 0.01, 2, n_ind)
            # threshold for the p% significance level according to Biltoft and Pardyjak (2009)
            thresh =  (2 * F_p) / (n_ind - 2 + 2 * F_p)

        # if a fixed threshold was chosen for the coherence
        elif thresh_type == 'coh':
            thresh = level

        else:
            raise ValueError('Threshold type can only be either \"sig\" or \"coh\"')

        # Find frequency where the coherence drops underneath the threshold for the first time.

        try:
            res = next(x for x, val in enumerate(coh[1:]) if val < thresh)
            freq_thresh = freq[res + 1]
            # turn this into a time scale
            tau[index] = 1 / freq_thresh
        except StopIteration:
            tau[index] = np.nan

    tau_spatial_scls = window_step * np.arange(1, num_windows + 1)

    return tau, tau_spatial_scls
