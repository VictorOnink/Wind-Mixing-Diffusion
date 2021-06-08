import settings
import matplotlib.pyplot as plt
import utils
import utils.utils_physics
from visualization import utils_visualization as utils_v
import numpy as np


def multiple_boundary_condition_comparison(selection='w_rise', close_up=None, output_step=-1,
                                  y_label='Depth (m)', x_label=r'Normalised Concentrations', fig_size=(10, 8),
                                  ax_label_size=16, legend_size=12, single_select=0, beaufort=4):
    """
    Figure comparing the vertical concentration profiles with Ceiling and Reflecting boundary conditions
    :param selection: selection criteria for loading parcels concentration profiles
    :param close_up: setting the y axis limit, with (max, min)
    :param output_step: determining which time index is shown, default is the concentration profile at the end of the
                        simulation
    :param y_label: label of the y axis
    :param x_label: label of the x axis
    :param fig_size: size of the figure
    :param ax_label_size: fontsize of the axis labels
    :param legend_size: fontsize of the legend
    :param single_select: selection index related to 'selection'
    :param beaufort: sea state/wind state which is being shown
    :return:
    """
    # Setting the axis limits
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=False)

    # Selecting the boundary conditions and rise velocities that we want to include in the figure
    boundary_list = ['Reflect', 'Ceiling']
    linestyle = ['-', '--']
    w_rise_list = [-0.03, -0.003, -0.0003]

    # Creating the axis
    shape = (1, 2)
    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size, shape=shape, plot_num=2,
                             all_x_labels=True)

    # Setting the wind speed
    mean_wind = np.mean(utils.utils_physics.beaufort_limits()[beaufort])

    for count, boundary in enumerate(boundary_list):
        # Plotting the distribution according to the SWB parametrization
        profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                      output_step=output_step, diffusion_type='SWB',
                                                      boundary=boundary, alpha_list=[0])
        for counter in range(len(profile_dict['concentration_list'])):
            ax[1].plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                       linestyle=linestyle[count], color=utils_v.return_color(counter))

        # Plotting the distribution according to the KPP parametrization
        profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                      output_step=output_step, diffusion_type='KPP',
                                                      boundary=boundary, alpha_list=[0])
        for counter in range(len(profile_dict['concentration_list'])):
            ax[0].plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                       linestyle=linestyle[count], color=utils_v.return_color(counter))

    # Adding the legend first for the color scheme (reflecting the rise velocity), then the linestyle (reflecting the
    # boundary condition)
    for index, w_r in enumerate(w_rise_list):
        ax[0].plot([], [], linestyle='-', color=utils_v.return_color(index),
                   label=r'$w_{rise}$ = ' + '{}'.format(abs(w_r)) + r' m s$^{-1}$')
    ax[0].legend(fontsize=legend_size, loc='lower right')
    for index, bound in enumerate(['Reflecting BC', 'Ceiling BC']):
        ax[1].plot([], [], linestyle=linestyle[index], color='k', label=bound)
    ax[1].legend(fontsize=legend_size, loc='lower right')

    # Adding subplot titles
    ax[0].set_title(r'(a) KPP', fontsize=ax_label_size)
    ax[1].set_title(r'(b) SWB', fontsize=ax_label_size)

    # Saving the figure
    plt.savefig(saving_filename_boundary(settings.figure_dir + '/Boundary Conditions/', close_up, beaufort),
                bbox_inches='tight', dpi=600)


def saving_filename_boundary(save_location, close_up, beafort):
    """ Setting the filename of the figure """
    if close_up is None:
        return save_location + 'Boundary_comparison_Bft={}.png'.format(beafort)
    else:
        ymax, ymin = close_up
        return save_location + 'Boundary_comparison_Bft={}_max={}_min={}.png'.format(beafort, ymax, ymin)