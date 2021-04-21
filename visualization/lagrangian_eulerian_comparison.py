import settings
import matplotlib.pyplot as plt
import utils
import visualization.utils_visualization
from visualization import utils_visualization as utils_v
import numpy as np


def lagrangian_eulerian_comparison(w_rise_list, alpha_list, output_step=-1, single_select=0, norm_depth=False,
                                   y_label='Depth (m)', close_up=None, x_label=r'Normalised Concentrations ($C/C_{max}$)',
                                   selection='w_rise', fig_size=(16, 20), ax_label_size=16, legend_size=10,
                                   boundary='Reflect'):
    if norm_depth:
        y_label = 'Depth/MLD'
        correction = settings.MLD
    else:
        correction = 1.0
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=norm_depth)

    # Get the base figure axis
    shape = (5, 2)
    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size, shape=shape, plot_num=10,
                             all_x_labels=True)

    for row in range(shape[0]):
        wind_range = utils.beaufort_limits()[row + 1]
        mean_wind = np.mean(wind_range)

        # Adding plot titles
        title_dict = {0: 'a', 1: 'b', 2: 'c', 3: 'd', 4: 'e', 5: 'f', 6: 'g', 7: 'h', 8: 'i', 9: 'j'}
        ax[2 * row + 0].set_title('({}) SWB, '.format(title_dict[2 * row + 0]) + r'$u_{10}$ = ' +
                                  '{:.2f}'.format(mean_wind) + r' m s$^{-1}$', fontsize=ax_label_size)
        ax[2 * row + 1].set_title('({}) KPP, '.format(title_dict[2 * row + 1]) + r'$u_{10}$ = ' +
                                  '{:.2f}'.format(mean_wind) + r' m s$^{-1}$', fontsize=ax_label_size)

        # Plotting the distribution according to the Kukulka parametrization, which goes in Axis 0
        profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                      output_step=output_step, diffusion_type='Kukulka',
                                                      boundary=boundary, alpha_list=alpha_list)
        for counter in range(len(profile_dict['concentration_list'])):
            _, w_r = profile_dict['parameter_kukulka'][counter]
            ax[2 * row + 0].plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'] / correction,
                                 label=line_labels(w_rise=w_r, boundary_type=boundary),
                                 linestyle='-', color=visualization.utils_visualization.return_color(counter))

        for counter, w_r in enumerate(np.abs(w_rise_list)):
            eul_dict = utils.load_obj(
                utils.get_eulerian_output_name(w_10=mean_wind, w_rise=w_r, diffusion_type='Kukulka'))
            ax[2 * row + 0].plot(eul_dict['C'], eul_dict['Z'] / correction, linestyle='--',
                                 color=visualization.utils_visualization.return_color(counter),
                                 label=line_labels(w_rise=w_r, boundary_type='eulerian'))

        # Plotting the distribution according to the KPP parametrization, which goes in Axis 1
        profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                      output_step=output_step, diffusion_type='KPP',
                                                      boundary=boundary, alpha_list=alpha_list)
        for counter in range(len(profile_dict['concentration_list'])):
            _, w_r = profile_dict['parameter_kukulka'][counter]
            ax[2 * row + 1].plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'] / correction,
                                 label=line_labels(w_rise=w_r, boundary_type=boundary),
                                 linestyle='-', color=visualization.utils_visualization.return_color(counter))

        for counter, w_r in enumerate(np.abs(w_rise_list)):
            eul_dict = utils.load_obj(utils.get_eulerian_output_name(w_10=mean_wind, w_rise=w_r, diffusion_type='KPP'))
            ax[2 * row + 1].plot(eul_dict['C'], eul_dict['Z'] / correction, linestyle='--',
                                 color=visualization.utils_visualization.return_color(counter),
                                 label=line_labels(w_rise=w_r, boundary_type='eulerian'))

    lines, labels = ax[1].get_legend_handles_labels()
    ax[0].legend(lines, labels, fontsize=legend_size, loc='lower right')

    plt.savefig(save_figure_name(alpha=alpha_list[0], boundary=boundary), bbox_inches='tight')


def line_labels(w_rise, boundary_type):
    w_rise = np.abs(w_rise)
    boundary_dict = {'Reflect': 'M-0', 'Reflect_Markov': 'M-1', 'eulerian': 'Eulerian'}
    boundary = boundary_dict[boundary_type]
    filename = '{}, '.format(boundary) + r'$w_r=$' + '{}'.format(w_rise) + r' m s$^{-1}$'
    return filename


def save_figure_name(boundary, alpha):
    boundary_dict = {'Reflect': 'M0', 'Reflect_Markov': 'M1'}
    dt = settings.dt_int.seconds
    if boundary_dict[boundary] is 'M0':
        filename = settings.figure_dir + 'Eulerian_comparison/euler_comp_{}_dt={}.png'.format(boundary_dict[boundary], dt)
    else:
        filename = settings.figure_dir + 'Eulerian_comparison/euler_comp_{}_alpha={}_dt={}.png'.format(boundary_dict[boundary], alpha, dt)
    return filename
