import settings
import matplotlib.pyplot as plt
from visualization import utils_visualization as utils_v
import utils
import numpy as np


def integration_timestep_control(w_rise_select, selection='w_10', close_up=None, y_label='Depth (m)', alpha=0.3,
                                x_label=r'Normalised Concentrations ($C/C_{max}$)', fig_size=(11, 6),
                                ax_label_size=16, legend_size=10, single_select=0,
                                output_step=-1, boundary='Reflect'):
    """
    Plotting the various concentration profiles with different integration timesteps, with the first subplot showing the
    effect with KPP diffusion and the second for SWB diffusion
    :param w_rise_select: the rise velocity which we want to use for the plot
    :param selection: we want to plot for a fixed wind condition
    :param close_up: setting the range of the x axis as (max min)
    :param y_label: the y axis label
    :param alpha: the memory term value
    :param x_label: the x axis label
    :param fig_size: the figure size
    :param ax_label_size: the fontsize of the axis labels
    :param legend_size: the fontsize of the legend
    :param single_select: the index for selecting the rise velocity
    :param output_step: we set this at -1, so we plot the last time step of the simulation
    :param boundary: setting the boundary condition
    :return:
    """
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=False)

    # Selecting which model data we want to plot based on the diffusion type
    swb, kpp, artificial = utils_v.boolean_diff_type('all')
    # Selecting the integration timesteps
    if boundary == 'Reflect':
        dt_list = [30, 15, 10, 5, 1]
    elif boundary == 'Reflect_Markov':
        dt_list = [30, 1]
    line_style = ['-'] * len(dt_list)
    w_rise = [w_rise_select] * len(dt_list)

    # Creating the axis, with one subplot for KPP and one for SWB
    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size, shape=(1, 2), plot_num=2,
                             all_x_labels=True)

    # Getting hte mean wind for Beaufort 4
    mean_wind = np.mean(utils.beaufort_limits()[4])

    # looping through all the timesteps we want to plot
    for count, dt in enumerate(dt_list):
        if kpp:
            profile_dict = utils_v.get_concentration_list([mean_wind], [w_rise[count]], selection, single_select,
                                                          output_step=output_step, diffusion_type='KPP',
                                                          boundary=boundary, alpha_list=[alpha], dt=dt)
            for counter in range(len(profile_dict['concentration_list'])):
                ax[0].plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                           label=label(w_rise[count], dt, boundary, alpha),
                           linestyle=line_style[count], color=utils_v.discrete_color_from_cmap(count, len(dt_list)))
        if swb:
            profile_dict = utils_v.get_concentration_list([mean_wind], [w_rise[count]], selection, single_select,
                                                          output_step=output_step, diffusion_type='SWB',
                                                          boundary=boundary, alpha_list=[alpha], dt=dt)
            for counter in range(len(profile_dict['concentration_list'])):
                ax[1].plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                           label=label(w_rise[count], dt, boundary, alpha),
                           linestyle=line_style[count], color=utils_v.discrete_color_from_cmap(count, len(dt_list)))

    lines, labels = ax[1].get_legend_handles_labels()
    # Adding the legends
    ax[1].legend(lines, labels, fontsize=legend_size, loc='lower right')

    # Adding subplot titles
    ax[0].set_title(r'(a) KPP', fontsize=ax_label_size)
    ax[1].set_title(r'(b) SWB', fontsize=ax_label_size)

    # Saving the figure
    diff_dict = {'Reflect': 'M0', 'Reflect_Markov': 'M1_{}'.format(alpha)}
    plt.savefig(settings.figure_dir + '{}_dt_int_check_w_rise={}_mld={}'.format(diff_dict[boundary], w_rise_select,
                                                                                settings.MLD) + '.png',
                bbox_inches='tight', dpi=600)


def label(w_rise, dt, boundary, alpha):
    if boundary == 'Reflect':
        return 'M-0, $w_r$ = ' + r'{}'.format(np.abs(w_rise)) + r' m s$^{-1}$, dt = ' + '{} s'.format(dt)
    else:
        return r'M-1, $\alpha=$ = ' + '{}, '.format(alpha) + '$w_r$ = ' + r'{}'.format(np.abs(w_rise)) + \
               r' m s$^{-1}$, dt = ' + '{} s'.format(dt)
