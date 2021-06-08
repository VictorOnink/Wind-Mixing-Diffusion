import settings
import matplotlib.pyplot as plt
import utils
import utils.utils_physics
import visualization.utils_visualization
from visualization import utils_visualization as utils_v
import numpy as np


def markov_alpha_dependence(w_rise_list, selection='w_10', close_up=None, y_label='Depth (m)', alpha=0.3,
                            x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(8, 8),
                            ax_label_size=16, legend_size=12, single_select=1,
                            output_step=-1, diffusion_type='SWB'):
    """
    A figure that shows how M-1 profiles with various alpha values compare with the M-0 profiles for Beaufort 4 conditions
    :param w_rise_list: list of rise velocities
    :param selection: selection variable for loading the rise velocity
    :param close_up: close up of the y axis as (max, min)
    :param y_label: label of the y axis
    :param alpha: opaqueness of the field data markers
    :param x_label: label of x axis
    :param fig_size: size of hte figure
    :param ax_label_size: fontsize of the axes labels
    :param legend_size: fontsize of the legend
    :param single_select: selection index related to 'selection'
    :param output_step: which time index we plot from the output, default is the final output time
    :param diffusion_type: which diffusion type to plot, either KPP or SWB
    :return:
    """
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=False)

    # Selecting which model data we want to plot based on the diffusion type
    swb, kpp, artificial = utils_v.boolean_diff_type(diffusion_type)
    # Selecting which model data we want to plot based on the diffusion scheme
    boundary_list = ['Reflect', 'Reflect_Markov', 'Reflect_Markov', 'Reflect_Markov', 'Reflect_Markov',
                     'Reflect_Markov', 'Reflect_Markov']
    alpha_list = [0, 0.0, 0.1, 0.3, 0.5, 0.7, 0.95]

    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size)
    line_style = ['-', '--', '--', '--', '--', '--', '--']

    # Loading the field data for Beaufort 4 wind conditions
    mean_wind = np.mean(utils.utils_physics.beaufort_limits()[4])
    _, _ = utils_v.add_observations(ax, norm_depth=False, alpha=alpha, wind_range=utils.utils_physics.beaufort_limits()[4])

    # Looping through all the M-0 and M-1 cases
    for count, boundary in enumerate(boundary_list):
        # Plotting the distribution according to the SWB diffusion profile
        if swb:
            profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                          output_step=output_step, diffusion_type='SWB',
                                                          boundary=boundary, alpha_list=alpha_list[count])
            for counter in range(len(profile_dict['concentration_list'])):
                ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                        label=utils_v.label_alpha_comparison(boundary=boundary, alpha=alpha_list[count]),
                        linestyle=line_style[count], color=visualization.utils_visualization.return_color(count))
        # Plotting the distribution according to the KPP diffusion profile
        if kpp:
            profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                          output_step=output_step, diffusion_type='KPP',
                                                          boundary=boundary, alpha_list=alpha_list[count])
            for counter in range(len(profile_dict['concentration_list'])):
                ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                        label=utils_v.label_alpha_comparison(boundary=boundary, alpha=alpha_list[count]),
                        linestyle=line_style[count], color=visualization.utils_visualization.return_color(count))
        # Plotting the distribution according to the artificial diffusion profile
        if artificial:
            profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                          output_step=output_step, diffusion_type='artificial',
                                                          boundary=boundary, alpha_list=alpha_list[count])
            for counter in range(len(profile_dict['concentration_list'])):
                ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                        label=utils_v.label_alpha_comparison(boundary=boundary, alpha=alpha_list[count]),
                        linestyle=line_style[count], color=visualization.utils_visualization.return_color(count))

    lines, labels = ax.get_legend_handles_labels()
    # Adding the legend
    ax.legend(lines[:-5], labels[:-5], fontsize=legend_size, loc='lower right')
    ax.set_title(r'u$_{10}$ = 5.4-7.9 m s$^{-1}$ - $\alpha$', fontsize=ax_label_size)

    # Saving the figure
    str_format = diffusion_type, w_rise_list[0], settings.MLD
    plt.savefig(settings.figure_dir + '{}_alpha_check_w_rise={}_mld={}'.format(*str_format) + '.png',
                bbox_inches='tight', dpi=600)