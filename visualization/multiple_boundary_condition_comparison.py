import settings
import matplotlib.pyplot as plt
import utils
import utils.utils_physics
from visualization import utils_visualization as utils_v
import numpy as np


def multiple_boundary_condition_comparison(selection='w_rise', close_up=None, output_step=-1,
                                  y_label='Depth (m)', x_label=r'Normalised Concentrations', fig_size=(10, 8),
                                  ax_label_size=16, legend_size=12, single_select=0, beaufort=4,
                                  interval=1, alpha=0.3):
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=False)

    # Selecting which model data we want to plot based on the diffusion scheme
    boundary_list = ['Reflect', 'Ceiling']
    linestyle = ['-', '--']
    w_rise_list = [-0.03, -0.003, -0.0003]

    # Creating the axis
    shape = (1, 2)
    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size, shape=shape, plot_num=2,
                             all_x_labels=True)

    mean_wind = np.mean(utils.utils_physics.beaufort_limits()[beaufort])
    # _, _ = utils_v.add_observations(ax, norm_depth=False, alpha=alpha, wind_range=utils.utils_physics.beaufort_limits()[4])

    for count, boundary in enumerate(boundary_list):
        # Plotting the distribution according to the Kukulka parametrization
        profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                      output_step=output_step, diffusion_type='Kukulka',
                                                      boundary=boundary, alpha_list=[0])
        for counter in range(len(profile_dict['concentration_list'])):
            w_rise = w_rise_list[counter]
            ax[1].plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                       linestyle=linestyle[count], color=utils_v.return_color(counter))

        # Plotting the distribution according to the KPP parametrization
        profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                      output_step=output_step, diffusion_type='KPP',
                                                      boundary=boundary, alpha_list=[0])
        for counter in range(len(profile_dict['concentration_list'])):
            w_rise = w_rise_list[counter]
            ax[0].plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                       linestyle=linestyle[count], color=utils_v.return_color(counter))

    # Adding the legend first for the color scheme, then for the boundary conditions
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
    plt.savefig(saving_filename_boundary(settings.figure_dir + '/Boundary Conditions/', selection, close_up, beaufort),
                bbox_inches='tight', dpi=600)


def saving_filename_boundary(save_location, selection, close_up, beafort):
    if close_up is None:
        return save_location + 'Boundary_comparison_Bft={}.png'.format(beafort)
    else:
        ymax, ymin = close_up
        return save_location + 'Boundary_comparison_Bft={}_max={}_min={}.png'.format(beafort, ymax, ymin)


def label_boundary(w_10, w_rise, diffusion_type, boundary):
    w_rise=np.abs(w_rise)
    boundary_dict = {'Mixed': 'M-0 Mixed', 'Reflect': 'M-0 Reflect', 'Reduce_dt': 'M-0 Reduce dt',
                     'Mixed_Markov': 'M-1 Mixed', 'Reflect_Markov':'M-1 Reflect',
                     'Reduce_dt_Markov':'M-1 Reduce dt', 'Ceiling': 'M-0 Ceiling'}
    if diffusion_type == 'Kukulka':
        return r'SWB, w$_{r}$ '+'= {}'.format(w_rise) + ' m s$^{-1}$, ' + boundary_dict[boundary]
    elif diffusion_type == 'KPP':
        return r'KPP, w$_{r}$ '+'= {}'.format(w_rise) + 'm s$^{-1}$, ' + boundary_dict[boundary]