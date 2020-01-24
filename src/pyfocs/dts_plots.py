import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator,
                               FormatStrFormatter)
import matplotlib.colors as colors
import numpy as np
from .xr_helper import xr_swap_dims_sel
import scipy


def bath_check(ds,
               bath_define,
               plot_var='power',
               bath_lims=[-0.005, 0.005],
               fig_kwargs=None, title=None):
    '''
    Plots for checking bath boundaries with both power anomaly and bias (bias
    plots not implemented yet)
    '''

    if fig_kwargs is None:
        fig_kwargs = dict()
    fig, axes = plt.subplots(1, len(bath_define),
                             sharey=True,
                             figsize=(10, 6),
                             **fig_kwargs)
    if title:
        fig.suptitle(title)

    for bn_num, bn in enumerate(bath_define):
        temp_bath = xr_swap_dims_sel(ds, 'LAF', 'calibration', bn)
        bath_start = temp_bath.LAF.min()
        bath_end = temp_bath.LAF.max()
        bath = ds.sel(LAF=slice(bath_start - 2, bath_end + 2))

        if plot_var == 'power':
            power = np.log(bath.Ps / bath.Pas)
            mean_power = power.sel(LAF=slice(bath_start, bath_end)).mean(dim='LAF').mean(dim='time').values
            power = power.mean(dim='time')

        ax = axes[bn_num]
        ax.xaxis.grid(True, which='major')
        ax.set_title(bn)

        # Fill in the bath locations
        ax.fill_between([bath_start, bath_end], -20, 60,
                        edgecolor='k',
                        facecolor='0.9',
                        alpha=0.8)

        # Power perturbation
        ax.plot(power.LAF, power - mean_power)

        ax.set_xlim(bath_start - 1.5, bath_end + 1.5)
        ax.set_ylim(bath_lims)

        if ax.is_first_col():
            ax.set_ylabel(r'$log(\frac{Ps}{Pas}) - \overline{log(\frac{Ps}{Pas})}$')
        ax.set_xlabel('LAF (m)')
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    fig.autofmt_xdate()

    return fig


def bias_violin(ds, bath_define, plot_var='bias',
                fig_kwargs=None, title=None, plot_lims=[-0.5, 0.5]):
    '''
    Violinplots of bath biases showing the distribution over both time and
    space.
    '''
    if fig_kwargs is None:
        fig_kwargs = dict()
    fig, ax = plt.subplots(1, 1, figsize=(8, 6), **fig_kwargs)

    if title:
        ax.set_title(title)

    val_max = np.nanmax(plot_lims)
    val_min = np.nanmin(plot_lims)

    for bn_num, bn in enumerate(bath_define):
        bath = xr_swap_dims_sel(ds, 'LAF', 'calibration', bn)

        if plot_var == 'bias':
            val = (bath.cal_temp - bath[bath_define[bn]])
            val_max = np.max([val_max, val.max().values])
            val_min = np.min([val_min, val.min().values])

        violin_parts = ax.violinplot(val.values.flatten(),
                                     [bn_num],
                                     showmeans=True)
        ax.text(bn_num - 0.25, -0.45,
                '$\mu$={:01.2f} C'.format(val.mean(dim='time').mean(dim='LAF').values),
                bbox=dict(facecolor='white', alpha=0.5),
                )

        for pc in violin_parts['bodies']:
            pc.set_facecolor('0.5')
            pc.set_edgecolor('black')
            pc.set_alpha(1)
        for comp in violin_parts:
            if comp == 'bodies':
                continue
            violin_parts[comp].set_color('black')

    ax.set_ylim(val_min, val_max)
    ax.set_ylabel(r'Bias ($^{\circ}$C)')
    ax.set_xticks(np.arange(len(bath_define)))
    ax.set_xticklabels(bath_define.keys());

    return fig


def bath_validation(ds, bath_define, bath_lims=None, plot_var='bias',
                    fig_kwargs=None, title=None):
    '''
    Bart-like  plots of bath biases and power anomaly.
    '''
    if plot_var == 'bias':
        label_text = 'Bias (K)'
        if not bath_lims:
            bath_lims = [-0.5, 0.5]
    if plot_var == 'power':
        label_text = r'$log(\frac{Ps}{Pas}) - \overline{log(\frac{Ps}{Pas})}$ [-]'
        if not bath_lims:
            bath_lims = [-0.005, 0.005]

    if fig_kwargs is None:
        fig_kwargs = dict()
    fig = plt.figure(figsize=(12, 12), **fig_kwargs)

    if title:
        fig.suptitle(title, y=0.95)

    widths = [1., 5, 0.25]
    spec = fig.add_gridspec(ncols=3,
                            nrows=len(bath_define) + 1,
                            width_ratios=widths,
                            hspace=0.15, wspace=0.15,
                            )
    divnorm = colors.DivergingNorm(vmin=np.min(bath_lims),
                                   vcenter=0,
                                   vmax=np.max(bath_lims))

    # First, initiate all subplots
    ax_colorbar = fig.add_subplot(spec[1:, 2])
    ax_time = fig.add_subplot(spec[0, 1])
    ax_LAFs = []
    ax_maps = []
    for bn_num, bn in enumerate(bath_define):
        ax_LAFs.append(fig.add_subplot(spec[bn_num + 1, 0]))
        ax_maps.append(fig.add_subplot(spec[bn_num + 1, 1]))
    ax_maps[-1].get_shared_x_axes().join(ax_maps[-1], ax_time, )

    # Characterize each bath's calibration
    for bn_num, bn in enumerate(bath_define):

        # Handle the axes for this bath
        ax_LAF = ax_LAFs[bn_num]
        ax_map = ax_maps[bn_num]

        # Create the biases
        bath = xr_swap_dims_sel(ds, 'LAF', 'calibration', bn)
        bath_start = bath.LAF.min()
        bath_end = bath.LAF.max()
        if plot_var == 'bias':
            val = (bath.cal_temp - bath[bath_define[bn]])
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


def dts_loc_plot(ploc, phys_locs, ds, lin_fit=False, offset=50, c=None):
    '''
    Plot the given location's mean and standard deviation of log(Ps/Pas) and
    stokes and anti-stokes intensities. This is used to verify that a given
    section behaves correctly (e.g., checking for bend-dependent power losses)
    and for verifying section limits.

    INPUTS:

    OUTPUTS:

    '''

    # Limits of the section to be plotted, including an offset for
    # edge effects.
    if c:
        s_start = np.min(phys_locs[ploc]['LAF'][c]) - offset
        s_end = np.max(phys_locs[ploc]['LAF'][c]) + offset
    else:
        s_start = np.min(phys_locs[ploc]['LAF']) - offset
        s_end = np.max(phys_locs[ploc]['LAF']) + offset

    # Only continue if an LAF is found
    if np.isnan(s_start) or np.isnan(s_end):
        print('LAFs for ' + ploc + ' return a NaN value.')
        return
    fig, axes = plt.subplots(4, 1, figsize=(10, 15), sharex=True)

    # Derive the useful quantities
    sect = ds.sel(LAF=slice(s_start, s_end))
    sect['logPsPas'] = np.log(sect.Ps / sect.Pas)
    sect = sect.mean(dim='time')
    sect = sect.std(dim='time')

    # Just the data within the section-of-interest for linear fitting.
    fit_sect = sect.sel(LAF=slice(np.min(phys_locs[ploc]['LAF'][c]),
                                  np.max(phys_locs[ploc]['LAF'][c])))
    fit_sect = fit_sect.mean(dim='time')

    # Mean(Log(Ps/Pas))
    ax = axes[0]
    ax.plot(sect.LAF,
            sect['logPsPas'].values)

    if lin_fit:
        logPsPas_m, logPsPas_b, _, _, _ = scipy.stats.linregress(fit_sect.LAF,
                                                                 fit_sect['logPsPas'].values)
        ax.plot(sect.LAF, logPsPas_b + logPsPas_m * sect.LAF.values, '--',
                label='linear fit', color='0.5')
        ax.legend()

    # Axis limits
    y_mean = sect['logPsPas'].mean(dim='LAF').values
    y_min = sect['logPsPas'].min().values + 0.004
    y_max = sect['logPsPas'].max().values + 0.001

    # locations labels
    for ploc_labels in phys_locs:
        lims = phys_locs[ploc_labels]['LAF'][c]
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
    ax.xaxis.grid(True, which='minor')
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.xaxis.set_minor_locator(MultipleLocator(2))

    # Std(Log(Ps/Pas))
    ax = axes[1]
    ax.plot(sect.LAF,
            sect['logPsPas'].values)

    y_mean = sect['logPsPas'].mean(dim='LAF').values
    y_min = sect['logPsPas'].min().values
    y_max = sect['logPsPas'].max().values

    for ploc_labels in phys_locs:
        lims = phys_locs[ploc_labels]['LAF'][c]
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
                        edgecolor='k',
                        facecolor='0.9',
                        alpha=0.8)

    ax.set_ylim(y_min, y_max)
    ax.set_ylabel('std(log(Ps/Pas))')
    ax.xaxis.grid(True, which='minor')
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.xaxis.set_minor_locator(MultipleLocator(2))

    # Mean Stokes/Anti-stokes intensities
    ax = axes[2]
    ax.plot(sect.LAF,
            np.log(sect['Ps'].values),
            label='Ps')
    ax.plot(sect.LAF,
            np.log(sect['Pas'].values),
            label='Pas')

    if lin_fit:
        Ps_m, Ps_b, _, _, _ = scipy.stats.linregress(fit_sect.LAF,
                                                     np.log(fit_sect['Ps'].values))
        Pas_m, Pas_b, _, _, _ = scipy.stats.linregress(fit_sect.LAF,
                                                       np.log(fit_sect['Pas'].values))
        ax.plot(sect.LAF, Ps_b + Ps_m * sect.LAF.values, '--',
                label='Ps fit')
        ax.plot(sect.LAF, Pas_b + Pas_m * sect.LAF.values, '--',
                label='Pas fit')

    y_mean = np.log(sect['Ps'].mean(dim='LAF').values)
    y_min = np.log(sect['Pas'].min().values) - 0.1
    y_max = np.log(sect['Ps'].max().values) + 0.1

    for ploc_labels in phys_locs:
        lims = phys_locs[ploc_labels]['LAF'][c]
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
    ax.set_ylabel(r'$P_s, P_{as}$ Backscatter Intensities [-]')
    ax.legend(loc='lower left')
    ax.xaxis.grid(True, which='minor')
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.xaxis.set_minor_locator(MultipleLocator(2))

    # STD Stokes/Anti-stokes intensities
    ax = axes[3]
    ax.plot(sect.LAF,
            np.log(sect['Ps'].values),
            label='Ps')
    ax.plot(sect.LAF,
            np.log(sect['Pas'].values),
            label='Pas')

    y_mean = np.log(sect['Ps'].mean(dim='LAF').values)
    y_min = np.min([np.log(sect['Pas'].min().values),
                    np.log(sect['Ps'].max().values)]) - 0.1
    y_max = np.max([np.log(sect['Pas'].min().values),
                    np.log(sect['Ps'].max().values)]) + 0.1
    for ploc_labels in phys_locs:
        lims = phys_locs[ploc_labels]['LAF'][c]
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
        ax.fill_between(lims, 0, 4000, edgecolor='k', facecolor='0.9', alpha=0.8)

    ax.xaxis.grid(True, which='minor')
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.xaxis.set_minor_locator(MultipleLocator(2))

    ax.legend(loc='lower left')
    ax.set_ylim(y_min, y_max)
    ax.set_ylabel(r'$\sigma(P_s), \sigma(P_{as})$')
    ax.set_xlim(s_start, s_end)
    ax.set_xlabel('LAF (m)')
    fig.suptitle(c + ': ' + ploc)
    fig.tight_layout()

    return fig
