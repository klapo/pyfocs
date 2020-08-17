import numpy as np
from scipy import stats
import xarray as xr
import matplotlib.pyplot as plt
import copy
import pyfocs


def section_shift_x(
    ds,
    fn_lib,
    ploc1,
    ploc2,
    label,
    fixed_shift=None,
    dl=4,
    lag=50,
    temp_field='cal_temp',
    plot_results=False,
    restep=None):
    '''
    Aligns two sections of fiber based on the maximum cross-correlation.
    Inputs:
        ds - xarray Dataset in the pyfocs format
        fn_lib - location library as formatted by pyfocs.check.config()
        ploc1 - string, name of the physical location to map ploc2 onto.
        ploc2 - string, name of the physical location to be mapped onto ploc1.
        label - string, name of the section to be aligned.
        fixed_shift - float, a fixed distance to shift ploc2.
        dl - float, estimate of artifact distance to use when aligning
        lag - float, number of shift indices, modified by restep.
        temp_field - string, name of the temperature field in ds.
        plot_results - boolean
        restep - int, Multiplicative factor for the number of steps per LAF to
            subsample both sections. This parameter allows finer shifting of
            ploc1 and ploc2 instead of by integer indices.
    Outputs:
        s1 - xarray Dataset, data mapped by ploc1 and label
        s2_shift - xarray Dataset, data aligned to ploc1 for ploc2 and label
        shift_x - float, meters that s2_shift was shifted by. Alignment at
            sub-dLAF steps is possible with the restep option.
        c - numpy array, cross-correlation of s1 and s2
        lags - numpy array, lags corresponding to the cross correlation
    '''
    lims_s1 = fn_lib[ploc1][label]['LAF']
    lims_s2 = fn_lib[ploc2][label]['LAF']

    fn_lib = copy.deepcopy(fn_lib)

    # Determine which sections need to be flipped
    if lims_s1[0] > lims_s1[1]:
        s1_reverse = True
        s1_slice = slice(lims_s1[1] - dl, lims_s1[0] + dl)
    else:
        s1_reverse = False
        s1_slice = slice(lims_s1[0] - dl, lims_s1[1] + dl)

    if lims_s2[0] > lims_s2[1]:
        s2_reverse = True
        s2_slice = slice(lims_s2[1] - dl, lims_s2[0] + dl)
    else:
        s2_reverse = False
        s2_slice = slice(lims_s2[0] - dl, lims_s2[1] + dl)

    # Get each section
    s1 = ds.sel(LAF=s1_slice)
    s2 = ds.sel(LAF=s2_slice)

    # Assign a relative 'x' index.
    if s1_reverse:
        s1.coords['x'] = (('LAF'), s1.LAF.max() - s1.LAF)
    elif not s1_reverse:
        s1.coords['x'] = (('LAF'), s1.LAF - s1.LAF.min())

    if s2_reverse:
        s2.coords['x'] = (('LAF'), s2.LAF.max() - s2.LAF)
    elif not s2_reverse:
        s2.coords['x'] = (('LAF'), s2.LAF - s2.LAF.min())

    # Force the two sections to be the same length in LAF
    if len(s1.LAF.values) > len(s2.LAF.values):
        s1 = s1.isel(LAF=slice(0, len(s2.LAF.values)))
    if len(s1.LAF.values) < len(s2.LAF.values):
        s2 = s2.isel(LAF=slice(0, len(s1.LAF.values)))

    # No fixed_shift was provided so we will determine the optimal shift
    # using the cross-correlation.
    if fixed_shift is None:
        # Interpolate to a fine 'x' scale to allow non-integer dLAF adjustments
        if restep:
            num_steps = len(s1.LAF.values) * restep
            x_new = np.linspace(s1.x[0].values, s1.x[-1].values, num_steps)
            s1_vals = s1.swap_dims({'LAF': 'x'})
            s1_vals = s1_vals.interp(x=x_new, method='linear')

            x_new = np.linspace(s2.x[0].values, s2.x[-1].values, num_steps)
            s2_vals = s2.swap_dims({'LAF': 'x'})
            s2_vals = s2_vals.interp(x=x_new, method='linear')

            # Use the time mean to remove the effect of instrument noise on the lag analysis.
            s1_vals = s1_vals[temp_field].mean(dim='time').values
            s2_vals = s2_vals[temp_field].mean(dim='time').values

        else:
            # Use the time mean to remove the effect of instrument noise on the lag analysis.
            s1_vals = s1[temp_field].mean(dim='time').values
            s2_vals = s2[temp_field].mean(dim='time').values
            restep = 1

        if s1_reverse:
            s1_vals = np.flip(s1_vals)
        if s2_reverse:
            s2_vals = np.flip(s2_vals)

        # Find the index with the maximum correlation
        if not lag:
            lag = int(np.min([len(s1_vals) * restep, len(s2_vals) * restep]) // 2)
        lags = np.arange(-lag, lag + 1)
        if s2_reverse:
            lags = np.flip(lags)
        c = pyfocs.stats.norm_xcorr(
            s1_vals,
            s2_vals,
            scaleopt='coef',
            remove_mean=True,
            lag=lag)
        max_corr = np.argmax(c)
        shift = lags[max_corr]
        dlaf, _ = stats.mode(np.diff(ds.LAF.values))
        shift_x = shift * dlaf / restep

    # The fixed_shift was provided as a set value.
    elif fixed_shift is not None:
        dlaf, _ = stats.mode(np.diff(ds.LAF.values))
        shift_x = fixed_shift
        restep = 1
        c = 0
        lags = np.arange(-lag, lag + 1)

    # Convert lags from indices to meters
    lags = lags * dlaf / restep

    # Shift the relative index to align the two sections
    s2_shift = copy.deepcopy(s2)
    s2_shift['x'] = s2_shift.x + shift_x

    s1.attrs['reverse'] = s1_reverse
    s2_shift.attrs['reverse'] = s2_reverse

    if plot_results:
        fig, axes = plt.subplots(3, 1, figsize=(8, 8))
        ax = axes[0]
        ax.plot(lags, c)

        ax.set_ylabel('Correlation [-]')
        ax.set_xlabel('Shift [m]')
        ax.plot([shift_x, shift_x], [c.min(), c.max()], color='0.5')
        ax.set_ylim(c.min(), c.max())
        ax.set_xlim(np.min(lags), np.max(lags))
        ax.text(0.15, 0.85, 'max correlation at shift ' + str(shift_x), transform=ax.transAxes)
        ax.set_title('a) Cross-correlation')

        ax = axes[1]
        ax.plot(s1.x, s1[temp_field].mean(dim='time'), marker='.', label=ploc1)
        ax.plot(s2.x, s2[temp_field].mean(dim='time'), marker='.', label=ploc2)
        ax.set_xlabel('Relative LAF [m]')
        ax.set_ylabel('Temperature [C]')
        ax.legend(loc='upper right')
        ax.set_ylim(s1[temp_field].min(), s1[temp_field].max())
        ax.set_title('b) Original temperature')
        ax.set_xlim(0, np.max([s1.x.max(), s2.x.max()]))

        ax = axes[2]
        ax.plot(s1.x, s1[temp_field].mean(dim='time'), marker='.', label=ploc1)
        ax.plot(s2_shift.x, s2_shift[temp_field].mean(dim='time'), marker='.', label=ploc2)
        ax.set_xlabel('Relative LAF [m]')
        ax.set_ylabel('Temperature [C]')
        ax.legend(loc='upper right')
        ax.set_ylim(s1[temp_field].min(), s1[temp_field].max())
        ax.set_title('c) Shifted temperature')
        ax.set_xlim(0, np.max([s1.x.max(), s2.x.max()]))

        fig.tight_layout()

    return s1, s2_shift, shift_x, c, lags


def interp_section(
    ds,
    fn_lib,
    ploc1,
    ploc2,
    label,
    fixed_shift=None,
    dl=4,
    plot_results=False,
    temp_field='cal_temp',
    restep=5,
    fig_kwargs=None):
    '''
    Aligns two sections and interpolates them to a common index.
    Inputs:
        ds - xarray Dataset in the pyfocs format
        fn_lib - location library as formatted by pyfocs.check.config()
        ploc1 - string, name of the physical location to map ploc2 onto.
        ploc2 - string, name of the physical location to be mapped onto ploc1.
        label - string, name of the section to be aligned.
        dl - float, estimate of artifact distance to use when aligning
        temp_field - string, name of the temperature field in ds.
        plot_results - boolean
        restep - int, allows finer shifting of ploc1 and ploc2 instead of by
            integer indices.
        fig_kwards - dictionary, keywords passed to pyplot.subplots()
    Outputs:
        sub_s1 - xarray Dataset, data mapped by ploc1 and label
        sub_s2 - xarray Dataset, data mapped by ploc2 and label
            sub_s1 and sub_s2 will share a common index, 'x'
        fn_lib - pyfocs location library, LAF values adjusted to match aligned
            and interpolated sections. Note: These LAF values will not
            necessarily correspond to equal length sections.
    '''
    fn_lib = copy.deepcopy(fn_lib)

    # Determine the shift to apply
    s1, s2, shift_x, c, lags = section_shift_x(ds, fn_lib, ploc1, ploc2,
                                               label, restep=restep, dl=dl,
                                               temp_field=temp_field,
                                               plot_results=False,
                                               fixed_shift=fixed_shift)

    # In order to avoid flipping the orientation of s2 to match s1
    # through the interpolation step we use a "reverse flag".
    s1_reverse = s1.attrs['reverse']
    s2_reverse = s2.attrs['reverse']

    # The actual section for verifying a uniform length
    loc_s1 = fn_lib[ploc1][label]['LAF']
    sub_s1 = s1.sel(LAF=slice(np.min(loc_s1), np.max(loc_s1)))

    # Interpolate s2 to s1's x coordinates for the section
    s1 = s1.swap_dims({'LAF': 'x'})
    s2 = s2.swap_dims({'LAF': 'x'})
    sub_s2 = s2.interp(x=sub_s1.x.values).swap_dims({'x': 'LAF'})

    # The interp command causes s2 to take on s1's properties.
    sub_s2.attrs['reverse'] = s2.attrs['reverse']

    # Return the adjusted LAF values based on aligning s2 to s1
    LAFmin = float(sub_s2.LAF.min().values)
    LAFmax = float(sub_s2.LAF.max().values)
    if s2_reverse:
        fn_lib[ploc2][label]['LAF'] = [LAFmax, LAFmin]
    else:
        fn_lib[ploc2][label]['LAF'] = [LAFmin, LAFmax]

    if plot_results:
        if fig_kwargs is None:
            fig_kwargs = dict()
        if 'figsize' not in fig_kwargs:
            fig_kwargs['figsize'] = (10, 4)

        # No cross-correlation is perfomed when providing a fixed_shift
        if fixed_shift is not None:
            num_rows = 2
            fig, axes = plt.subplots(num_rows, 1, **fig_kwargs)
        else:
            num_rows = 3
            fig, axes = plt.subplots(num_rows, 1, **fig_kwargs)

        if 'time' in s1.coords:
            s1 = s1.mean(dim='time')
            sub_s1_mean = sub_s1.mean(dim='time')
        if 'time' in s2.coords:
            s2 = s2.mean(dim='time')
            sub_s2_mean = sub_s2.mean(dim='time')

        ax_ind = 0
        # A cross-correlation was performed, we will include it in the plot.
        if fixed_shift is None:
            ax = axes[ax_ind]
            ax_ind = ax_ind + 1
            ax.plot(lags, c)

            ax.set_ylabel('Correlation [-]')
            ax.set_xlabel('Shift [m]')
            ax.plot([shift_x, shift_x], [c.min(), c.max()], color='0.5')
            ax.set_ylim(c.min(), c.max())
            ax.set_xlim(np.min(lags), np.max(lags))
            ax.text(0.15, 0.85, 'max correlation at shift ' + str(shift_x) + 'm', transform=ax.transAxes)
            ax.set_title('a) Cross-correlation between ' + ploc1 + ' and ' + ploc2, loc='left')

        # Interpolated and shifted data
        ax = axes[0 + ax_ind]
        delta_T = (s1[temp_field].max() - s1[temp_field].min()) / 10

        colr1_interp = 'xkcd:blue'
        colr1_orig = 'xkcd:light blue'

        colr2_interp = 'xkcd:orange'
        colr2_orig = 'xkcd:light orange'

        ax.plot(s1.x, s1[temp_field].T, color=colr1_orig)
        ax.plot(s2.x, s2[temp_field].T, color=colr2_orig)

        ax.plot(sub_s2_mean.x, sub_s1_mean[temp_field].T, marker='o', label=ploc1, color=colr1_interp)
        ax.plot(sub_s2_mean.x, sub_s2_mean[temp_field].T, marker='o', label=ploc2 + ' (interpolated)', color=colr2_interp)

        ax.plot([np.min(sub_s1.x.values),
                 np.min(sub_s1.x.values)],
                [s1[temp_field].min().values - delta_T,
                 s1[temp_field].max().values + delta_T],
                'k-')
        ax.plot([np.max(sub_s1.x.values),
                 np.max(sub_s1.x.values)],
                [s1[temp_field].min().values - delta_T,
                 s1[temp_field].max().values + delta_T],
                'k-')

        ax.set_title('b) ' + label, loc='left')
        ax.set_xlabel('Common x index [m]')
        ax.set_ylabel('Temperature [C]')
        ax.legend()
        text = '\n'.join([ploc1 + ': ' + str(loc_s1),
                          ploc2 + ': ' + str(fn_lib[ploc2][label]['LAF']),
                          ploc1 + ' length=' + str(len(sub_s1.LAF.values)),
                          ploc2 + 'length=' + str(len(sub_s2.LAF.values))])
        ax.text(0.15, 0.7, text, transform=ax.transAxes, bbox={'facecolor': 'white', 'alpha': 0.5})
        ax.set_xlim(0, s1.x.max())
        ax.set_ylim(s1[temp_field].min() - delta_T, s1[temp_field].max() + delta_T)

        # Spatial derivative for verifying section boundaries
        ax = axes[1 + ax_ind]
        ax.plot(s1.differentiate(coord='x').x,
                s1[temp_field].differentiate(coord='x').T,
                marker='o', label=ploc1, color=colr1_orig)
        ax.plot(s2.differentiate(coord='x').x,
                s2[temp_field].differentiate(coord='x').T,
                marker='o', label=ploc2, color=colr2_orig)

        ax.plot(sub_s1_mean.differentiate(coord='x').x,
                sub_s1_mean[temp_field].differentiate(coord='x').T,
                marker='o', label=ploc1, color=colr1_interp)
        ax.plot(sub_s2_mean.differentiate(coord='x').x,
                sub_s2_mean[temp_field].differentiate(coord='x').T,
                marker='o', label=ploc2, color=colr2_interp)

        ax.plot([np.min(sub_s1.x.values),
                 np.min(sub_s1.x.values)],
                [-1.5, 1.5],
                'k-')
        ax.plot([np.max(sub_s1.x.values),
                 np.max(sub_s1.x.values)],
                [-1.5, 1.5],
                'k-')

        ax.set_title('c) ' + label + ' derivative in LAF', loc='left')
        ax.set_xlabel('Common x index [m]')
        ax.set_ylabel(r'$\frac{dT}{dz}$ [C]')

        ax.set_xlim(0, s1.x.max())
        ax.set_ylim(-1.5, 1.5)

        fig.tight_layout()

    return sub_s1, sub_s2, fn_lib


def section_limits(
    s1,
    s2,
    fn_lib,
    ploc1,
    ploc2,
    label,
    plot_results=False,
    temp_field='cal_temp'):
    '''
    Similar to interp_section, but without interpolation. Requires output from
    `section_shift_x`. Useful only for the first view of mapping locations.
    '''

    loc_s1 = fn_lib[ploc1][label]['LAF']

    s1_reverse = s1.attrs['reverse']
    s2_reverse = s2.attrs['reverse']

    # The actual section for verifying a uniform length
    sub_s1 = s1.sel(LAF=slice(np.min(loc_s1), np.max(loc_s1)))

    # Get the values of the new section limits for s2 based on index x
    # x can be reversed so we need to account for that here.
    if s2_reverse:
        loc_s2 = [np.max(sub_s1.x.values), np.min(sub_s1.x.values)]
        sub_s2 = s2.swap_dims({'LAF': 'x'}).sel(x=slice(np.max(sub_s1.x.values), np.min(sub_s1.x.values)))
    else:
        loc_s2 = [np.min(sub_s1.x.values), np.max(sub_s1.x.values)]
        sub_s2 = s2.swap_dims({'LAF': 'x'}).sel(x=slice(np.min(sub_s1.x.values), np.max(sub_s1.x.values)))

    # Return the adjusted LAF values based on aligning s2 to s1
    fn_lib[ploc2][label]['LAF'] = [float(sub_s2.LAF.isel(x=0).values),
                                   float(sub_s2.LAF.isel(x=-1).values)]

    if plot_results:
        if 'time' in s1.coords:
            s1 = s1.mean(dim='time')
        if 'time' in s2.coords:
            s2 = s2.mean(dim='time')

        delta_T = (s1[temp_field].max() - s1[temp_field].min()) / 10

        fig, axes = plt.subplots(2, 1, figsize=(18, 6))
        ax = axes[0]
        ax.plot(s1.x, s1[temp_field].T, marker='o', label=ploc1)
        ax.plot(s2.x, s2[temp_field].T, marker='o', label=ploc2)

        ax.plot([np.min(sub_s1.x.values),
                 np.min(sub_s1.x.values)],
                [s1[temp_field].min().values - delta_T,
                 s1[temp_field].max().values + delta_T],
                'k-')
        ax.plot([np.max(sub_s1.x.values),
                 np.max(sub_s1.x.values)],
                [s1[temp_field].min().values - delta_T,
                 s1[temp_field].max().values + delta_T],
                'k-')

        ax.set_title('a) ' + label, loc='left')
        ax.set_xlabel('Common x index [m]')
        ax.set_ylabel('Temperature [C]')
        ax.legend()
        ax.text(0.15, 0.85, ploc1 + ': ' + str(loc_s1), transform=ax.transAxes)
        ax.text(0.15, 0.8, ploc2 + ': ' + str(fn_lib[ploc2][label]['LAF']), transform=ax.transAxes)
        ax.text(0.15, 0.75, ploc1 + ' length=' + str(len(sub_s1.LAF.values)), transform=ax.transAxes)
        ax.text(0.15, 0.7, ploc2 + 'length=' + str(len(sub_s2.LAF.values)), transform=ax.transAxes)
        ax.set_xlim(0, s1.x.max())
        ax.set_ylim(s1[temp_field].min() - delta_T, s1[temp_field].max() + delta_T)

        ax = axes[1]
        ax.plot(s1.differentiate(coord='x').x,
                s1[temp_field].differentiate(coord='x').T,
                marker='o', label=ploc1)
        ax.plot(s2.differentiate(coord='x').x,
                s2[temp_field].differentiate(coord='x').T,
                marker='o', label=ploc2)

        ax.plot([np.min(sub_s1.x.values),
                 np.min(sub_s1.x.values)],
                [-1.5, 1.5],
                'k-')
        ax.plot([np.max(sub_s1.x.values),
                 np.max(sub_s1.x.values)],
                [-1.5, 1.5],
                'k-')

        ax.set_title('b) ' + label + ' derivative in LAF', loc='left')
        ax.set_xlabel('Common x index [m]')
        ax.set_ylabel(r'$\frac{dT}{dz}$ [C]')

        ax.set_xlim(0, s1.x.max())
        ax.set_ylim(-1.5, 1.5)
    return fn_lib
