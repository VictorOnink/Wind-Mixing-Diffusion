import matplotlib.pyplot as plt
import utils, settings
from visualization import utils_visualization as utils_v
import numpy as np


def diffusion_markov_comparison(w_rise_list, selection='w_10', close_up=None, y_label='Depth (m)', alpha=0.15,
                                x_label=r'Normalised Concentrations', fig_size=(10, 6),
                                ax_label_size=16, legend_size=10, single_select=0,
                                output_step=-1, wave_roughness=False, add_variability=True, theta=1.0):
    """
    This creates a figure comparing the M-0 simulations with M-1 simulations with various values of alpha
    The first subfigure shows the results for KPP diffusion, while the second shows the results for SWB diffusion
    :param add_variability: if True, add variability shading
    :param w_rise_list: list of rise velocities
    :param selection: we are plotting for a fixed rise velocity
    :param close_up: defining the range of the x axis with (max, min)
    :param y_label: the label of y axis
    :param alpha: setting the opacity of the field observation markers
    :param x_label: the label of the x axis
    :param fig_size: the size of the figure
    :param ax_label_size: the fontsize of hte axis labels
    :param legend_size: the fontsize of hte lengend
    :param single_select: the index of the rise velocity list we pick for the plot
    :param output_step: selecting which time index we are plotting hte concentration for (default is the last)
    :param wave_roughness: if True, have surface roughness be wave height dependent
    :param theta: Langmuir circulation amplification factor
    :return:
    """
    # Getting the axis ranges
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=False)

    # Selecting which model data we want to plot based on the diffusion type
    swb, kpp, _ = utils_v.boolean_diff_type('all')
    # Selecting which model data we want to plot based on the diffusion scheme
    boundary_list = ['Ceiling', 'Ceiling_Markov', 'Ceiling_Markov', 'Ceiling_Markov', 'Ceiling_Markov',
                     'Ceiling_Markov', 'Ceiling_Markov']
    alpha_list = [0, 0.0, 0.1, 0.3, 0.5, 0.7, 0.95]
    line_style = ['-', '--', '--', '--', '--', '--', '--']

    # Creating the axis, one for KPP and one for Kukulka
    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size, shape=(1, 2), plot_num=2,
                             all_x_labels=True)

    # Setting the wind speed, and adding the field observations to the axis
    mean_wind = np.mean(utils.beaufort_limits()[4])
    for axis in ax:
        data_line, data_label = utils_v.add_observations(axis, norm_depth=False, alpha=alpha,
                                                         wind_range=utils.beaufort_limits()[4])

    # First, plotting the KPP data:
    for count, boundary in enumerate(boundary_list):
        if kpp:
            profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                          output_step=output_step, diffusion_type='KPP',
                                                          boundary=boundary, alpha_list=alpha_list[count],
                                                          wave_roughness=wave_roughness, theta=theta)
            for counter in range(len(profile_dict['concentration_list'])):
                concentration = profile_dict['concentration_list'][counter]
                depth = profile_dict['depth_bins']
                ax[0].plot(concentration, depth,
                           label=utils_v.label_alpha_comparison(boundary=boundary, alpha=alpha_list[count]),
                           linestyle=line_style[count], color=utils_v.return_color(count))
                if add_variability:
                    std = profile_dict['std_list'][counter]
                    upper_limit, lower_limit = concentration + std, concentration - std
                    ax[0].fill_betweenx(depth, lower_limit, upper_limit, alpha=0.2, color=utils_v.return_color(count))

        if swb:
            profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                          output_step=output_step, diffusion_type='SWB',
                                                          boundary=boundary, alpha_list=alpha_list[count],
                                                          wave_roughness=wave_roughness, theta=theta)
            for counter in range(len(profile_dict['concentration_list'])):
                concentration = profile_dict['concentration_list'][counter]
                depth = profile_dict['depth_bins']
                ax[1].plot(concentration, depth,
                           label=utils_v.label_alpha_comparison(boundary=boundary, alpha=alpha_list[count]),
                           linestyle=line_style[count], color=utils_v.return_color(count))
                if add_variability:
                    std = profile_dict['std_list'][counter]
                    upper_limit, lower_limit = concentration + std, concentration - std
                    ax[1].fill_betweenx(depth, lower_limit, upper_limit, alpha=0.2, color=utils_v.return_color(count))

    lines, labels = ax[1].get_legend_handles_labels()

    # Adding the legends
    ax[0].legend(data_line, data_label, fontsize=legend_size, loc='lower right')
    ax[1].legend(lines[:-6], labels[:-6], fontsize=legend_size, loc='lower right')

    # Adding subplot titles
    ax[0].set_title(r'(a) KPP', fontsize=ax_label_size)
    ax[1].set_title(r'(b) SWB', fontsize=ax_label_size)

    # Saving the figure
    str_format = w_rise_list[0], settings.MLD, settings.dt_int.seconds
    plt.savefig(settings.figure_dir + 'markov_diffusion_check_w_rise={}_mld={}_dt={}'.format(*str_format) + '.png',
                bbox_inches='tight', dpi=600)