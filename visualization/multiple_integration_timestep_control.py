import settings
import matplotlib.pyplot as plt
import utils
from visualization import utils_visualization as utils_v
import numpy as np


def multiple_integration_timestep_control(beaufort, selection='w_10', y_label='Depth (m)', alpha=0.3,
                                x_label=r'Normalised Concentrations', fig_size=(8, 12),
                                ax_label_size=15, legend_size=10, single_select=0,
                                output_step=-1, boundary='Ceiling'):
    """
    figure comparing the influence of the integration timestep on the vertical concentration profiles. We do this for
    the three different rise velocities, and since these get mixed to varying depths, each subplot has its own y axis
    limits
    :param beaufort: sea state/wind speed on the Beaufort scale
    :param selection: selection criteria for loading parcels data
    :param y_label: label of y axis
    :param alpha: memory term for M-1 processes
    :param x_label: label of the x axis
    :param fig_size: size of the figure
    :param ax_label_size: fontsize of the axes labels
    :param legend_size: fontsize of the legend
    :param single_select: selection index related to 'selection'
    :param output_step: set at -1, so we plot the final concentration profile in the parcels simulation
    :param boundary: which boundary condition do we want to show
    :return:
    """
    # Getting the axis limits, but the y axis limits will be updated later
    ax_range = utils_v.get_axes_range(close_up=None, norm_depth=False)

    # Selecting which model data we want to plot based on the diffusion type
    swb, kpp, _ = utils_v.boolean_diff_type('all')

    # Selecting the timesteps for which we want to plot the data
    if boundary in ['Ceiling', 'Reflect']:
        dt_list = [30, 15, 10, 5, 1]
    elif boundary in ['Ceiling_Markov', 'Reflect_Markov']:
        dt_list = [30, 1]
    # Setting the linestyles and rise velocities
    line_style = ['-'] * len(dt_list)
    w_rise = [-0.03, -0.003, -0.0003]

    # Creating the axis, one for KPP and one for Kukulka
    shape = (3, 2)
    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size, shape=(3, 2), plot_num=6,
                             all_x_labels=True)

    # Getting the wind speed corresponding to the Beaufort scale
    mean_wind = np.mean(utils.utils_physics.beaufort_limits()[beaufort])

    # Setting the y limits in order of decreasing rise velocities
    ylim = [(-2, 0.1), (-15, 0.1), (-20, 0.1)]

    # The even column values correspond to the first column, the odd column values correspond to the second column
    for column in range(shape[0]):
        for count, dt in enumerate(dt_list):
            # Plotting the KPP profiles
            if kpp:
                profile_dict = utils_v.get_concentration_list([mean_wind], [w_rise[column]], selection, single_select,
                                                              output_step=output_step, diffusion_type='KPP',
                                                              boundary=boundary, alpha_list=[alpha], dt=dt)
                for counter in range(len(profile_dict['concentration_list'])):
                    ax[2 * column].plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                               label=label(dt),
                               linestyle=line_style[count], color=utils_v.discrete_color_from_cmap(count, len(dt_list)))
                    ax[2 * column].set_ylim(ylim[column])
            # Plotting the SWB profiles
            if swb:
                profile_dict = utils_v.get_concentration_list([mean_wind], [w_rise[column]], selection, single_select,
                                                              output_step=output_step, diffusion_type='Kukulka',
                                                              boundary=boundary, alpha_list=[alpha], dt=dt)
                for counter in range(len(profile_dict['concentration_list'])):
                    ax[2 * column + 1].plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                               label=label(dt),
                               linestyle=line_style[count], color=utils_v.discrete_color_from_cmap(count, len(dt_list)))
                    ax[2 * column + 1].set_ylim(ylim[column])

    lines, labels = ax[1].get_legend_handles_labels()
    # Adding the legend in the bottom right subfigure
    ax[-1].legend(lines, labels, fontsize=legend_size, loc='lower right')

    # Adding subplot titles
    ax[0].set_title(r'(a) KPP, $w_{r}=0.03$ m s$^{-1}$', fontsize=ax_label_size)
    ax[1].set_title(r'(b) SWB, $w_{r}=0.03$ m s$^{-1}$', fontsize=ax_label_size)
    ax[2].set_title(r'(c) KPP, $w_{r}=0.003$ m s$^{-1}$', fontsize=ax_label_size)
    ax[3].set_title(r'(d) SWB, $w_{r}=0.003$ m s$^{-1}$', fontsize=ax_label_size)
    ax[4].set_title(r'(e) KPP, $w_{r}=0.0003$ m s$^{-1}$', fontsize=ax_label_size)
    ax[5].set_title(r'(f) SWB, $w_{r}=0.0003$ m s$^{-1}$', fontsize=ax_label_size)

    # Saving the figure
    diff_dict = {'Ceiling': 'M0-Ceiling', 'Ceiling_Markov': 'M1_{}-Ceiling'.format(alpha),
                 'Reflect': 'M0-Reflect', 'Reflect_Markov': 'M1_{}-Reflect'.format(alpha)}
    plt.savefig(settings.figure_dir + 'Bft_{}_{}_dt_int_check_mld={}'.format(beaufort, diff_dict[boundary],
                                                                                settings.MLD) + '.png',
                bbox_inches='tight', dpi=600)


def label(dt):
    """ Setting the label for the integration timestep """
    return 'dt = {} s'.format(dt)