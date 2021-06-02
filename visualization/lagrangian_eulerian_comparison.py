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
    """
    Here we plot comparisons of the Eulerian and Lagrangian parcels concentration profiles for different wind conditions
    and rise velocities
    :param w_rise_list: the list of rise velocities
    :param alpha_list: the list of alpha values for the M-1 simulations
    :param output_step: set at -1, so we plot the final concentration profile in the parcels simulation
    :param single_select: selection variable in loading the parcels concentration profiles
    :param norm_depth: if True, the depths are all normalized by the MLD
    :param y_label: y axis label
    :param close_up: setting the limits of the y axis as (max min)
    :param x_label: x axis label
    :param selection: is w_rise since in each subplot we are plotting a range of w_rise values for one wind condition
    :param fig_size: figure size
    :param ax_label_size: fontsize of axis labels
    :param legend_size: fontsize of legend
    :param boundary: which boundary condition we are plotting for, and whether M-0 or M-1
    :return:
    """
    if norm_depth:
        y_label = 'Depth/MLD'
        correction = settings.MLD
    else:
        correction = 1.0
    # Getting the y and x axis limits
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=norm_depth)

    # Get the base figure axis, with each row corresponding to a different surface wind speed, and first column being
    # for SWB diffusion, and the second for KPP diffusion
    shape = (5, 2)
    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size, shape=shape, plot_num=10,
                             all_x_labels=True)

    # Looping through the rows, which correspond to increasing with speeds
    for row in range(shape[0]):
        # Selecting the wind speed
        wind_range = utils.beaufort_limits()[row + 1]
        mean_wind = np.mean(wind_range)

        # Adding plot titles
        title_dict = {0: 'a', 1: 'b', 2: 'c', 3: 'd', 4: 'e', 5: 'f', 6: 'g', 7: 'h', 8: 'i', 9: 'j'}
        ax[2 * row + 0].set_title('({}) SWB, '.format(title_dict[2 * row + 0]) + r'$u_{10}$ = ' +
                                  '{:.2f}'.format(mean_wind) + r' m s$^{-1}$', fontsize=ax_label_size)
        ax[2 * row + 1].set_title('({}) KPP, '.format(title_dict[2 * row + 1]) + r'$u_{10}$ = ' +
                                  '{:.2f}'.format(mean_wind) + r' m s$^{-1}$', fontsize=ax_label_size)

        # Plotting the distribution according to the SWB parametrization, which goes in Axis 0
        # First, the parcels concentration profile
        profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                      output_step=output_step, diffusion_type='SWB',
                                                      boundary=boundary, alpha_list=alpha_list)
        for counter in range(len(profile_dict['concentration_list'])):
            _, w_r = profile_dict['parameter_kukulka'][counter]
            ax[2 * row + 0].plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'] / correction,
                                 label=line_labels(w_rise=w_r, boundary_type=boundary),
                                 linestyle='-', color=visualization.utils_visualization.return_color(counter))
        # Next, the Eulerian concentration profile
        for counter, w_r in enumerate(np.abs(w_rise_list)):
            eul_dict = utils.load_obj(
                utils.get_eulerian_output_name(w_10=mean_wind, w_rise=w_r, diffusion_type='SWB'))
            ax[2 * row + 0].plot(eul_dict['C'], eul_dict['Z'] / correction, linestyle='--',
                                 color=visualization.utils_visualization.return_color(counter),
                                 label=line_labels(w_rise=w_r, boundary_type='eulerian'))

        # Plotting the distribution according to the KPP parametrization, which goes in Axis 1, with first the parcels
        # concentration profile
        profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                      output_step=output_step, diffusion_type='KPP',
                                                      boundary=boundary, alpha_list=alpha_list)
        for counter in range(len(profile_dict['concentration_list'])):
            _, w_r = profile_dict['parameter_kukulka'][counter]
            ax[2 * row + 1].plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'] / correction,
                                 label=line_labels(w_rise=w_r, boundary_type=boundary),
                                 linestyle='-', color=visualization.utils_visualization.return_color(counter))
        # Next, the Eulerian concentration profile
        for counter, w_r in enumerate(np.abs(w_rise_list)):
            eul_dict = utils.load_obj(utils.get_eulerian_output_name(w_10=mean_wind, w_rise=w_r, diffusion_type='KPP'))
            ax[2 * row + 1].plot(eul_dict['C'], eul_dict['Z'] / correction, linestyle='--',
                                 color=visualization.utils_visualization.return_color(counter),
                                 label=line_labels(w_rise=w_r, boundary_type='eulerian'))
    # adding in a lagend to the first subplot
    lines, labels = ax[1].get_legend_handles_labels()
    ax[0].legend(lines, labels, fontsize=legend_size, loc='lower right')

    # Saving the figure
    plt.savefig(save_figure_name(alpha=alpha_list[0], boundary=boundary), bbox_inches='tight')


def line_labels(w_rise, boundary_type):
    """
    Labelling the lines in the subplots based on the boundary type (which temporarily includes the Eulerian option) and
    the rise velocity
    """
    w_rise = np.abs(w_rise)
    boundary_dict = {'Reflect': 'M-0', 'Reflect_Markov': 'M-1', 'eulerian': 'Eulerian'}
    boundary = boundary_dict[boundary_type]
    filename = '{}, '.format(boundary) + r'$w_r=$' + '{}'.format(w_rise) + r' m s$^{-1}$'
    return filename


def save_figure_name(boundary, alpha):
    """
    Setting the figure file name depending on whether we were comparing the Eulerian model results with M-0 or M-1
    simulations, and for M-1 simulation specifying the alpha term
    :param boundary:
    :param alpha:
    :return:
    """
    boundary_dict = {'Reflect': 'M0', 'Reflect_Markov': 'M1'}
    dt = settings.dt_int.seconds
    if boundary_dict[boundary] is 'M0':
        filename = settings.figure_dir + 'Eulerian_comparison/euler_comp_{}_dt={}.png'.format(boundary_dict[boundary], dt)
    else:
        filename = settings.figure_dir + 'Eulerian_comparison/euler_comp_{}_alpha={}_dt={}.png'.format(boundary_dict[boundary], alpha, dt)
    return filename
