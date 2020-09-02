import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator,
                               FormatStrFormatter)
from matplotlib.patches import Patch
import matplotlib.colors as colors
import numpy as np
from .xr_helper import swap_sel
import scipy
import importlib
import pyfocs
import copy

from matplotlib.lines import Line2D


def bias_violin(ds,
                in_callib,
                plot_var='bias',
                fig_kwargs=None,
                title=None,
                plot_lims=None,
                temp_field_name='cal_temp'):
    '''
    Violinplots of bath biases showing the distribution over both time and
    space.
    '''
    # Prep the calibration library by removing baths not in the data.
    callib = copy.deepcopy(in_callib)
    for bn_num, bn in enumerate(in_callib):
        try:
            _ = swap_sel(ds, 'LAF', 'calibration', bn)
        except KeyError:
            del callib[bn]

    if fig_kwargs is None:
        fig_kwargs = dict()
    if 'figsize' not in fig_kwargs:
        fig_kwargs['figsize'] = (6, 4)
    fig, ax = plt.subplots(1, 1, **fig_kwargs)

    if title:
        ax.set_title(title)

    colr_cal = 'xkcd:purple'
    colr_val = 'xkcd:orange'

    for bn_num, bn in enumerate(callib):
        try:
            bath = swap_sel(ds, 'LAF', 'calibration', bn)
        except KeyError:
            # @ Add re-labeling plus a warning here.
            continue

        ref = callib[bn]['ref_sensor']
        ref_type = callib[bn]['type']
        if ref_type == 'calibration':
            colr = colr_cal
        elif ref_type == 'validation':
            colr = colr_val

        if plot_var == 'bias':
            val = (bath[temp_field_name] - bath[ref])
            val_max = val.max().values
            val_min = val.min().values

        violin_parts = ax.violinplot(val.values.flatten(),
                                     [bn_num],
                                     showmeans=True)
        ax.text(bn_num - 0.25, -0.45,
                r'$\mu$={:01.2f} C'.format(val.mean(dim='time').mean(dim='LAF').values),
                bbox=dict(facecolor='white', alpha=0.5),
                )

        for pc in violin_parts['bodies']:
            pc.set_facecolor(colr)
            pc.set_edgecolor('black')
            pc.set_alpha(1)
        for comp in violin_parts:
            if comp == 'bodies':
                continue
            violin_parts[comp].set_color('black')

    if not plot_lims:
        ax.set_ylim(val_min, val_max)
    else:
        ax.set_ylim(np.min(plot_lims), np.max(plot_lims))
    ax.set_ylabel(r'Bias ($^{\circ}$C)')
    ax.set_xticks(np.arange(len(callib)))
    ax.set_xticklabels(callib.keys())

    # Custom legend for validation vs calibration bath types
    legend_elements = [Patch(facecolor=colr_cal,
                             edgecolor='k',
                             label='cal'),
                       Patch(facecolor=colr_val,
                             edgecolor='k',
                             label='val'),]
    ax.legend(handles=legend_elements, loc='best')

    return fig


def bath_validation(ds, in_callib, bath_lims=None, plot_var='bias',
                    fig_kwargs=None, title=None, temp_field_name='cal_temp'):
    '''
    Bart-like  plots of bath biases and power anomaly.
    '''
    callib = copy.deepcopy(in_callib)
    # Prep the calibration library by removing baths not in the data.
    for bn_num, bn in enumerate(in_callib):
        try:
            _ = swap_sel(ds, 'LAF', 'calibration', bn)
        except KeyError:
            del callib[bn]

    if plot_var == 'bias':
        label_text = 'Bias (K)'
        if not bath_lims:
            bath_lims = [-0.5, 0.5]
    if plot_var == 'power':
        label_text = r'$<log(\frac{Ps}{Pas})>$ [-]'
        if not bath_lims:
            bath_lims = [-0.005, 0.005]

    if fig_kwargs is None:
        fig_kwargs = dict()

    if 'figsize' not in fig_kwargs:
        fig_kwargs['figsize'] = (6, 4)
    fig = plt.figure(**fig_kwargs)

    if title:
        fig.suptitle(title, y=0.95)

    widths = [1., 5, 0.25]
    spec = fig.add_gridspec(ncols=3,
                            nrows=len(callib) + 1,
                            width_ratios=widths,
                            hspace=0.18, wspace=0.15,
                            )
    divnorm = colors.DivergingNorm(vmin=np.min(bath_lims),
                                   vcenter=0,
                                   vmax=np.max(bath_lims))

    # First, initiate all subplots
    ax_colorbar = fig.add_subplot(spec[1:, 2])
    ax_time = fig.add_subplot(spec[0, 1])
    ax_LAFs = []
    ax_maps = []
    for bn_num, bn in enumerate(callib):
        ax_LAFs.append(fig.add_subplot(spec[bn_num + 1, 0]))
        ax_maps.append(fig.add_subplot(spec[bn_num + 1, 1]))
    ax_maps[-1].get_shared_x_axes().join(ax_maps[-1], ax_time, )

    # Characterize each bath's calibration
    for bn_num, bn in enumerate(callib):
        ref = callib[bn]['ref_sensor']
        ref_type = callib[bn]['type']

        # Handle the axes for this bath
        ax_LAF = ax_LAFs[bn_num]
        ax_map = ax_maps[bn_num]

        # Create the biases
        bath = swap_sel(ds, 'LAF', 'calibration', bn)
        bath_start = bath.LAF.min()
        bath_end = bath.LAF.max()
        if plot_var == 'bias':
            val = (bath[temp_field_name] - bath[ref])
        # Create the power anomaly
        if plot_var == 'power':
            val = bath.logPsPas - bath.logPsPas.mean(dim='LAF')

        # Time-averaged variable, explicit over LAF
        ax_LAF.plot(val.mean(dim='time'), val.LAF)
        ax_LAF.plot([0, 0], [bath_start, bath_end], 'k')
        ax_LAF.set_ylim(bath_start, bath_end)
        ax_LAF.set_xlim(np.min(bath_lims), np.max(bath_lims))
        ax_LAF.set_ylabel('LAF (m)')

        if ax_LAF.is_last_row():
            ax_LAF.set_xlabel(label_text)
        else:
            ax_LAF.set_xticklabels([])

        # Map of biases
        im = ax_map.pcolormesh(val.time, val.LAF, val.T,
                               cmap=plt.get_cmap('RdBu'),
                               norm=divnorm)
        ax_map.set_yticklabels([])

        if not ax_map.is_last_row():
            ax_map.set_xticklabels([])
        if not ax_LAF.is_last_row():
            ax_LAF.set_xticklabels([])
        ax_map.get_shared_y_axes().join(ax_map, ax_LAF)

        # Annotation
        ax_map.text(0.0, 0.8, bn, transform=ax_map.transAxes,
                    bbox=dict(facecolor='white', alpha=0.5))

        # LAF-averaged, explicit over time
        ax_time.plot(val.time, val.mean(dim='LAF'), label=bn)
        ax_time.set_ylabel(label_text)
        ax_time.set_ylim(np.min(bath_lims), np.max(bath_lims))

    # Colorbar
    plt.colorbar(im, cax=ax_colorbar, extend='both')
    ax_colorbar.set_ylabel(label_text)

    # Join the time axis
    ax_time.legend(bbox_to_anchor=(1.0, 1.08))
    ax_time.set_xticklabels([])
    ax_time.set_xlim(val.time.min().values, val.time.max().values)

    fig.autofmt_xdate()

    return fig


def dts_loc_plot(ploc,
                 phys_locs,
                 ds,
                 lin_fit=False,
                 offset=50,
                 fig_kwargs=None):
    '''
    Plot the given location's mean and standard deviation of log(Ps/Pas) and
    stokes and anti-stokes intensities. This is used to verify that a given
    section behaves correctly (e.g., checking for bend-dependent power losses)
    and for verifying section limits.

    INPUTS:
        ploc - string that is a key within the phys_locs location dictionary
               returned by the check.config() function. This should indicate
               section for plotting.
        phys_locs - location library dictionary returned by the check.config()
        ds - xarray Dataset as formatted by pyfocs.
        offset - distance in meters added to the section LAF limits.
        lin_fig - indicates whether a linear fit should be made between
                  adjacent sections when testing for step- or bend-losses.
    OUTPUTS:

    '''

    # Limits of the section to be plotted
    start = np.min(phys_locs[ploc]['LAF'])
    end = np.max(phys_locs[ploc]['LAF'])
    # Limits with an offset to examine edge effects and bend losses at holders.
    s_start = start - offset
    s_end = end + offset

    if not fig_kwargs:
        fig_kwargs = dict()

    if 'figsize' not in fig_kwargs:
        fig_kwargs['figsize'] = (6, 12)

    # Only continue if an LAF is found
    if np.isnan(s_start) or np.isnan(s_end):
        print('LAFs for ' + ploc + ' returned a NaN value.')
        return
    fig, axes = plt.subplots(3, 1, sharex=True, **fig_kwargs)

    # Derive the useful quantities
    sect = ds.sel(LAF=slice(s_start, s_end))
    sect['logPsPas'] = np.log(sect.Ps / sect.Pas)
    sect_mean = sect.mean(dim='time')
    sect_ps_var = sect_mean['Ps'].rolling(LAF=10, center=True).std()
    sect_pas_var = sect_mean['Pas'].rolling(LAF=50, center=True).std()

    # Just the data within the section-of-interest for linear fitting.
    fit_sect = sect.sel(LAF=slice(start, end))
    fit_sect = fit_sect.mean(dim='time')

    # Mean(Log(Ps/Pas))
    ax = axes[0]
    ax.plot(sect_mean.LAF,
            sect_mean['logPsPas'].values)

    if lin_fit:
        logPsPas_m, logPsPas_b, _, _, _ = scipy.stats.linregress(fit_sect.LAF,
                                                                 fit_sect['logPsPas'].values)
        ax.plot(sect.LAF, logPsPas_b + logPsPas_m * sect.LAF.values, '--',
                label='linear fit', color='0.5')
        ax.legend()

    # Axis limits
    y_mean = sect_mean['logPsPas'].mean(dim='LAF').values
    y_min = sect_mean['logPsPas'].min().values - 0.004
    y_max = sect_mean['logPsPas'].max().values + 0.004

    # locations labels
    for ploc_labels in phys_locs:
        lims = phys_locs[ploc_labels]['LAF']
        # Only label points that are within the plotted limits.
        if (np.isnan(lims).any()
                or np.max(lims) < s_start
                or np.min(lims) > s_end):
            continue

        # Text locations for the label
        text_y_loc = y_mean + (y_max - y_mean) / 2
        if np.max(lims) > s_end:
            text_x_loc = np.min(lims)
        elif np.min(lims) < s_start:
            text_x_loc = s_start
        else:
            text_x_loc = np.mean(lims)

        # Label the location
        ax.text(text_x_loc, text_y_loc, ploc_labels)
        ax.fill_between(lims, 0, 4000,
                        edgecolor='k',
                        facecolor='0.9',
                        alpha=0.8)

    # Axis formatting
    ax.set_ylim(y_min, y_max)
    ax.set_ylabel('log(Ps/Pas)')
    ax.xaxis.grid(True, which='both')
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.xaxis.set_minor_locator(MultipleLocator(2))

    # Mean log(Stokes) and log(Anti-stokes)
    ax = axes[1]
    ax.plot(sect_mean.LAF,
            np.log(sect_mean['Ps'].values),
            label='Ps')
    ax.plot(sect_mean.LAF,
            np.log(sect_mean['Pas'].values),
            label='Pas')

    if lin_fit:
        Ps_m, Ps_b, _, _, _ = scipy.stats.linregress(fit_sect.LAF,
                                                     np.log(fit_sect['Ps'].values))
        Pas_m, Pas_b, _, _, _ = scipy.stats.linregress(fit_sect.LAF,
                                                       np.log(fit_sect['Pas'].values))
        ax.plot(sect.LAF, Ps_b + Ps_m * sect.LAF.values, '--',
                label=r'$log(P_s)$ fit')
        ax.plot(sect.LAF, Pas_b + Pas_m * sect.LAF.values, '--',
                label=r'$log(P_{as})$ fit')

    y_mean = np.log(sect_mean['Ps'].mean(dim='LAF').values)
    y_min = np.log(sect_mean['Pas'].min().values) - 0.1
    y_max = np.log(sect_mean['Ps'].max().values) + 0.1

    for ploc_labels in phys_locs:
        lims = phys_locs[ploc_labels]['LAF']
        if np.isnan(lims).any() or np.max(lims) < s_start or np.min(lims) > s_end:
            continue

        text_y_loc = y_mean + (y_max - y_mean) / 2

        if np.max(lims) > s_end:
            text_x_loc = np.min(lims)
        elif np.min(lims) < s_start:
            text_x_loc = s_start
        else:
            text_x_loc = np.mean(lims)
        ax.text(text_x_loc, text_y_loc, ploc_labels)
        ax.fill_between(lims, 0, 4000,
                        edgecolor='k',
                        facecolor='0.9',
                        alpha=0.8)
    ax.set_ylim(y_min, y_max)
    ax.set_ylabel(r'$log(P_s), log(P_{as})$ Backscatter Intensities [-]')
    ax.legend(loc='lower left')
    ax.xaxis.grid(True, which='both')
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.xaxis.set_minor_locator(MultipleLocator(2))

    # Spatial variances in stokes and anti-stokes intensities
    ax = axes[2]
    var_ratio = (sect_ps_var / sect_pas_var)
    ax.plot(sect_ps_var.LAF,
            var_ratio.values)

    y_mean = var_ratio.mean(dim='LAF').values
    y_min = np.max([var_ratio.min().values - 0.1, 0])
    y_max = var_ratio.max().values + 0.1
    for ploc_labels in phys_locs:
        lims = phys_locs[ploc_labels]['LAF']
        if (np.isnan(lims).any()
                or np.max(lims) < s_start
                or np.min(lims) > s_end):
            continue

        text_y_loc = y_mean + (y_max - y_mean) / 2

        if np.max(lims) > s_end:
            text_x_loc = np.min(lims)
        elif np.min(lims) < s_start:
            text_x_loc = s_start
        else:
            text_x_loc = np.mean(lims)
        ax.text(text_x_loc, text_y_loc, ploc_labels)
        ax.fill_between(lims, 0, 4000,
                        edgecolor='k', facecolor='0.9', alpha=0.8)

    ax.xaxis.grid(True, which='both')
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.xaxis.set_minor_locator(MultipleLocator(2))

    ax.set_ylim(y_min, y_max)
    ax.set_ylabel('$\frac{\sigma_{x}(P_s)}{\sigma_{x}(P_{as})}$ [-]')
    ax.set_xlim(s_start, s_end)
    ax.set_xlabel('LAF (m)')
    fig.tight_layout()

    return fig


def bath_check(ds,
               callib,
               lims_dict=None,
               fig_kwargs=None,
               title=None,
               include_temp=False,
               laf_lims=2,
               ):
    '''
    Helper function for determining bath limits and integrity.

    INPUTS:
        ds - xarray object following pyfocs protocol
        callib - The calibration library dictionary from pyfocs.config.check
        laf_lims - LAF padding to include for visual inspection, default is 2m.
        include_temp - a boolean flag indicating if temperature should be
            plottted.
        title - string for the suptitle
        fig_kwargs - dictionary of figure keywords
    RETURNS:
        figure handle for the created plot
    '''

    # 1st row = bias in instrument reported temperature (if requested)
    # 2nd row = power anomaly
    # 3rd row = standard deviation of power
    # 4th row = spatial derivative of power

    # I like purple and orange
    colr_cal = 'xkcd:purple'
    colr_val = 'xkcd:orange'

    # Determine if there are any figure keywords to pass
    if fig_kwargs is None:
        fig_kwargs = dict()
    if 'figsize' not in fig_kwargs:
        fig_kwargs['figsize'] = (2, 6)
    # Determine if we should include bias in
    # instrument reported temperature.
    if include_temp:
        numrows = 4
    else:
        numrows = 3

    numcols = len(callib)

    # Generate the figure
    fig, axes = plt.subplots(numrows,
                             numcols,
                             **fig_kwargs)
    axes = np.atleast_2d(axes)
    print(np.shape(axes))
    if title:
        fig.suptitle(title)

    # Limits for the different plots
    bath_lims = {'temp': [-2, 2],  # +/- from mean ref. temp
                 'power_anom': [-0.05, 0.05],
                 'power_std': [0, 0.015],
                 'power_deriv': [-0.05, 0.05],
                 }
    if lims_dict:
        # Sanitize the input dictionary
        subset_test = all([ld in bath_lims.keys() for ld in lims_dict])
        assert subset_test
        # Re-assign from lims_dict to bath_lims
        for ld in lims_dict:
            bath_lims[ld] = lims_dict[ld]

    # Determine dimension name, swap to calibration, and select bath.
    spacedim = [d for d in ds.dims if 'time' not in d]
    msg = 'Too many non-time dimensions detected. Expected only a 2D array.'
    assert len(spacedim) == 1, msg
    spacedim = spacedim[0]

    # Grab just this bath
    for bnum, bname in enumerate(callib):
        probename = callib[bname]['ref_sensor']
        ref_type = callib[bname]['type']
        if ref_type == 'calibration':
            colr = colr_cal
        elif ref_type == 'validation':
            colr = colr_val
        axcol = axes[:, bnum]

        ds = ds.swap_dims({spacedim: 'calibration'})
        bath = ds.loc[{'calibration': bname}]
        bath_start = bath[spacedim].min().values
        bath_end = bath[spacedim].max().values
        ds = ds.swap_dims({'calibration': spacedim})
        bathplus = ds.sel({spacedim: slice(bath_start - laf_lims,
                                           bath_end + laf_lims)})

        # How is temperature named?
        if 'cal_temp' in ds.data_vars:
            temp_var = 'cal_temp'
        elif 'temp' in ds.data_vars:
            temp_var = 'temp'
        elif 'instr_temp' in ds.data_vars:
            temp_var = 'instr_temp'
        else:
            print('Expected to find cal_temp, instr_temp, or temp data variable.')
            raise KeyError

        # Plot the temperature bias if it was requested
        axind = 0
        if include_temp:
            ax = axcol[0]
            axind = 1
            ax.fill_between([bath_start, bath_end],
                            -20, 60, edgecolor='k',
                            facecolor='0.9', alpha=0.8)
            ax.plot(bathplus.LAF,
                    bathplus.mean(dim='time')[temp_var],
                    'o-', color=colr, label=temp_var)

            probeval = bath[probename].mean(dim='time')
            ax.plot([bath_start, bath_end],
                    [probeval, probeval],
                    '--', color='0.6', label='Reference')

            ylims = bath[probename].mean(dim='time').values
            ylims = ylims + np.array(bath_lims['temp'])
            ax.set_ylim(ylims)
            if ax.is_first_col():
                ax.legend()
                ax.set_ylabel('Temperature (C)')
            if ax.is_first_row():
                ax.set_title(bname)

        # Power anomaly
        ax = axcol[0 + axind]
        power = np.log(bathplus.Ps / bathplus.Pas)
        mean_power = power.sel(LAF=slice(bath_start, bath_end)).mean(dim='LAF').mean(dim='time').values
        power_anom = power - mean_power
        power_anom_std = power_anom.std(dim='time')
        power_anom_mean = power_anom.mean(dim='time')

        ax.xaxis.grid(True, which='major')

        ax.fill_between([bath_start, bath_end], -20, 60,
                        edgecolor='k', facecolor='0.9', alpha=0.8)
        ax.plot(power_anom_mean.LAF, power_anom_mean, 'o-', color=colr)
        ax.fill_between(power_anom_mean.LAF,
                        power_anom_mean - power_anom_std,
                        power_anom_mean + power_anom_std,
                        color=colr, alpha=0.5)

        ax.set_xlim(bath_start - laf_lims, bath_end + laf_lims)
        ax.set_ylim(bath_lims['power_anom'])
        if ax.is_first_col():
            ax.set_ylabel(r'$log(\frac{Ps}{Pas}) - \overline{log(\frac{Ps}{Pas})}$')
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        if ax.is_first_row():
            ax.set_title(bname)

        # Standard deviation of power
        ax = axcol[1 + axind]
        ax.fill_between([bath_start, bath_end], -20, 60,
                        edgecolor='k', facecolor='0.9', alpha=0.8)

        power_std = power.std(dim='time')
        ax.plot(power.LAF, power_std, 'o-', color=colr)

        ax.set_xlim(bath_start - laf_lims, bath_end + laf_lims)
        ax.set_ylim(bath_lims['power_std'])
        if ax.is_first_col():
            ax.set_ylabel(r'$\sigma_{time}(log(\frac{Ps}{Pas}))$')
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))

        # Spatial derivative of power
        ax = axcol[2 + axind]
        dP_dLAF = power.diff(dim='LAF').mean(dim='time')
        dP_dLAF_sigma = power.diff(dim='LAF').std(dim='time')

        ax.fill_between([bath_start, bath_end], -20, 60,
                        edgecolor='k', facecolor='0.9', alpha=0.8)

        ax.plot(dP_dLAF.LAF, dP_dLAF, 'o-', color=colr)
        ax.fill_between(dP_dLAF.LAF,
                        dP_dLAF - dP_dLAF_sigma,
                        dP_dLAF + dP_dLAF_sigma,
                        color=colr, alpha=0.5)

        ax.set_xlim(bath_start - laf_lims, bath_end + laf_lims)
        ax.set_ylim(bath_lims['power_deriv'])
        if ax.is_first_col():
            ax.set_ylabel(r'$\frac{dlog(\frac{Ps}{Pas})}{dLAF}$')
        ax.set_xlabel('LAF (m)')
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))

    fig.autofmt_xdate()
    fig.tight_layout()

    return fig


def overview(ds,
             lims_dict=None,
             indexer='LAF',
             fig_kwargs=None,
             title=None,
             temp_field_name='cal_temp',
             plot_type='temperature'
            ):
    '''
    Helper function for plotting the time-averaged overview of the DTS data.
    INPUTS:
        ds - xarray object following pyfocs protocol
        plot_type [<'temperature'>, 'variance'] - specifies the type of plot
            'temperature' plot of instrument and calibrated temperature.
            'variance' plot running variance of the backscatter intensities.
        title - optional string, name for the plot
        indexer - name of the spatial indexer. If not provided assumes the
            the spatial indexer is called 'LAF'.
        temp_field_name - 'cal_temp' (default), 'instr_temp', or 'both'
        fig_kwards - keywords to pass to the figure creation command.
    '''

    # Determine if there are any figure keywords to pass
    if fig_kwargs is None:
        fig_kwargs = dict()
    if 'figsize' not in fig_kwargs:
        fig_kwargs['figsize'] = (6, 2)

    # Generate the figure
    fig, ax = plt.subplots(1, 1, **fig_kwargs)

    if title:
        fig.suptitle(title)

    # Limits for the different plots
    if lims_dict and indexer in lims_dict.keys():
        xlims = lims_dict[indexer]
        assert len(xlims) == 2, 'lims_dict was not properly formatted.'
        ds = ds.sel(ds=slice(xlims[0], xlims[-1]))
    else:
        xlims = [ds[indexer].min(), ds[indexer].max()]

    if lims_dict and plot_type in lims_dict:
        ylims = lims_dict[plot_type]
    elif plot_type == 'temperature':
        ylims = [ds[temp_field_name].min(),
                 ds[temp_field_name].max()]

    if plot_type == 'temperature':
        # Both instrument and calibrated temperatures
        if temp_field_name == 'both':
            ax.plot(ds.LAF, ds.instr_temp.mean(dim='time'),
                    label='Instrument reported temperature')
            ax.plot(ds.LAF, ds.cal_temp.mean(dim='time'),
                    label='Calibrated temperature')
        # Whatever temperature name the user passed
        else:
            ax.plot(ds.LAF, ds[temp_field_name].mean(dim='time'),
                    label=temp_field_name)
        ylabel = 'Temperature (C)'

    if plot_type == 'variance':
        ps_var = ds['Ps'].rolling(LAF=10, center=True).std().mean(dim='time')
        pas_var = ds['Pas'].rolling(LAF=10, center=True).std().mean(dim='time')
        if lims_dict and 'variance' in lims_dict:
            ylims = lims_dict['variance']
        else:
            ylims = [0, np.max([ps_var.max().values, pas_var.max().values])]

        ax.plot(ps_var[indexer], ps_var, label='Stokes')
        ax.plot(pas_var[indexer], pas_var, label='Anti-Stokes')
        ylabel = r'$\sigma_{LAF}$ (power) [-]'

    if plot_type == 'noise':
        print('Estimating instrument noise can take some time, especially for longer time slices')
        var = np.ones_like(ds.LAF) * np.nan
        noise = np.ones_like(ds.LAF) * np.nan

        for nlaf, laf in enumerate(ds.LAF):
            if np.isnan(ds.sel(LAF=laf)[temp_field_name].values).any():
                continue
            var[nlaf], _, noise[nlaf] = pyfocs.noisymoments(ds[temp_field_name].sel(LAF=laf).values)

        if lims_dict and 'noise' in lims_dict:
            ylims = lims_dict['noise']
        else:
            ylims = [0, 1]
        ax.plot(ds.LAF, noise, label='Estimated Instrument Noise')
        ylabel = 'Noise (K)'

    ax.legend()
    ax.set_xlabel(indexer + ' (m)')
    ax.set_xlim(xlims)
    ax.set_ylabel(ylabel)
    ax.set_ylim(ylims)

    return overview
