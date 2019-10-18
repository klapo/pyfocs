import matplotlib as plt
from matplotlib.ticker import (MultipleLocator,
                               FormatStrFormatter)
import matplotlib.colors as colors
import numpy as np
import .xr_helper


def bath_check(ds, bath_define, plot_var='power', bath_lims=[-0.005, 0.005],
               fig_kwargs=None, title=None):

    if fig_kwargs is None:
        fig_kwargs = dict()
    fig, axes = plt.subplots(1, len(bath_define),
                             sharey=True,
                             figsize=(10, 6),
                             **fig_kwargs)
    if title:
        fig.suptitle(title)

    for bn_num, bn in enumerate(bath_define):
        temp_bath = xr_helper.xr_swap_dims_sel(ds, 'LAF', 'calibration', bn)
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
        ax.fill_between([bath_start, bath_end], -20, 60, edgecolor='k', facecolor='0.9', alpha=0.8)

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


def bias_violin(ds, bath_define, plot_var='bias', fig_kwargs=None, title=None):

    if fig_kwargs is None:
        fig_kwargs = dict()
    fig, ax = plt.subplots(1, 1, figsize=(8, 6), **fig_kwargs)

    if title:
        ax.set_title(title)

    val_max = 0.5
    val_min = -0.5

    for bn_num, bn in enumerate(bath_define):
        bath = xr_swap_dims_sel(ds, 'LAF', 'calibration', bn)
        bath_start = bath.LAF.min()
        bath_end = bath.LAF.max()

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
    ax.set_ylabel('Bias ($^{\circ}$C)')
    ax.set_xticks(np.arange(len(bath_define)))
    ax.set_xticklabels(bath_define.keys());

    return fig


def bath_validation(ds, bath_define, bath_lims=None, plot_var='bias',
                    fig_kwargs=None, title=None):

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
                            left=0.08,
                            bottom=0.12,
                            right=0.9,
                            top=0.88)
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
    ax_time.legend()
    ax_time.set_xticklabels([])
    ax_time.set_xlim(val.time.min().values, val.time.max().values)

    fig.autofmt_xdate()

    return fig
