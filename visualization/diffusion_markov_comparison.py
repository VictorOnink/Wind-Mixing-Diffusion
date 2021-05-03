import settings
import matplotlib.pyplot as plt
import utils
import utils.utils_physics
import visualization.utils_visualization
from visualization import utils_visualization as utils_v
import numpy as np


def diffusion_markov_comparison(w_rise_list, selection='w_10', close_up=None, y_label='Depth (m)', alpha=0.3,
                                x_label=r'Normalised Concentrations', fig_size=(10, 6),
                                ax_label_size=16, legend_size=10, single_select=1,
                                output_step=-1):
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=False)

    # Selecting which model data we want to plot based on the diffusion type
    kukulka, kpp, artificial = utils_v.boolean_diff_type('all')
    # Selecting which model data we want to plot based on the diffusion scheme
    boundary_list = ['Ceiling', 'Ceiling_Markov', 'Ceiling_Markov', 'Ceiling_Markov', 'Ceiling_Markov',
                     'Ceiling_Markov', 'Ceiling_Markov']
    alpha_list = [0, 0.0, 0.1, 0.3, 0.5, 0.7, 0.95]
    line_style = ['-', '--', '--', '--', '--', '--', '--']

    # Creating the axis, one for KPP and one for Kukulka
    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size, shape=(1, 2), plot_num=2,
                             all_x_labels=True)

    mean_wind = np.mean(utils.utils_physics.beaufort_limits()[4])
    for axis in ax:
        data_line, data_label = utils_v.add_observations(axis, norm_depth=False, alpha=alpha,
                                                         wind_range=utils.utils_physics.beaufort_limits()[4])

    # First, plotting the KPP data:
    for count, boundary in enumerate(boundary_list):
        if kpp:
            profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                          output_step=output_step, diffusion_type='KPP',
                                                          boundary=boundary, alpha_list=alpha_list[count])
            for counter in range(len(profile_dict['concentration_list'])):
                ax[0].plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                           label=utils_v.label_alpha_comparison(boundary=boundary, alpha=alpha_list[count]),
                           linestyle=line_style[count], color=visualization.utils_visualization.return_color(count))
        if kukulka:
            profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                          output_step=output_step, diffusion_type='Kukulka',
                                                          boundary=boundary, alpha_list=alpha_list[count])
            for counter in range(len(profile_dict['concentration_list'])):
                ax[1].plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                           label=utils_v.label_alpha_comparison(boundary=boundary, alpha=alpha_list[count]),
                           linestyle=line_style[count], color=visualization.utils_visualization.return_color(count))

    lines, labels = ax[1].get_legend_handles_labels()
    # Adding the legends
    ax[0].legend(data_line, data_label, fontsize=legend_size, loc='lower right')
    ax[1].legend(lines[:-5], labels[:-5], fontsize=legend_size, loc='lower right')

    # Adding subplot titles
    ax[0].set_title(r'(a) KPP', fontsize=ax_label_size)
    ax[1].set_title(r'(b) SWB', fontsize=ax_label_size)

    # Saving the figure
    str_format = w_rise_list[0], settings.MLD, settings.dt_int.seconds
    plt.savefig(settings.figure_dir + 'markov_diffusion_check_w_rise={}_mld={}_dt={}'.format(*str_format) + '.png',
                bbox_inches='tight', dpi=600)