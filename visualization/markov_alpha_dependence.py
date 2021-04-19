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
                            output_step=-1, diffusion_type='Kukulka'):
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=False)

    # Selecting which model data we want to plot based on the diffusion type
    kukulka, kpp, artificial = utils_v.boolean_diff_type(diffusion_type)
    # Selecting which model data we want to plot based on the diffusion scheme
    boundary_list = ['Reflect', 'Reflect_Markov', 'Reflect_Markov', 'Reflect_Markov', 'Reflect_Markov',
                     'Reflect_Markov', 'Reflect_Markov']
    alpha_list = [0, 0.0, 0.1, 0.3, 0.5, 0.7, 0.95]

    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size)
    line_style = ['-', '--', '--', '--', '--', '--', '--']

    mean_wind = np.mean(utils.utils_physics.beaufort_limits()[4])
    _, _ = utils_v.add_observations(ax, norm_depth=False, alpha=alpha, wind_range=utils.utils_physics.beaufort_limits()[4])

    for count, boundary in enumerate(boundary_list):
        # Plotting the distribution according to the KPP parametrization
        if kukulka:
            profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                          output_step=output_step, diffusion_type='Kukulka',
                                                          boundary=boundary, alpha_list=alpha_list[count])
            for counter in range(len(profile_dict['concentration_list'])):
                ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                        label=utils_v.label_alpha_comparison(boundary=boundary, alpha=alpha_list[count]),
                        linestyle=line_style[count], color=visualization.utils_visualization.return_color(count))
        if kpp:
            profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                          output_step=output_step, diffusion_type='KPP',
                                                          boundary=boundary, alpha_list=alpha_list[count])
            for counter in range(len(profile_dict['concentration_list'])):
                ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                        label=utils_v.label_alpha_comparison(boundary=boundary, alpha=alpha_list[count]),
                        linestyle=line_style[count], color=visualization.utils_visualization.return_color(count))
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