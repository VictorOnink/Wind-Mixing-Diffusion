import settings
import matplotlib.pyplot as plt
import utils
import utils.utils_physics
import visualization.utils_visualization
from visualization import utils_visualization as utils_v
import numpy as np


def boundary_condition_comparison(w_rise_list, alpha_list, selection='w_rise', close_up=None, output_step=-1,
                                  y_label='Depth (m)', x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(8, 8),
                                  ax_label_size=16, legend_size=12, single_select=0, beaufort=4,
                                  diffusion_type='KPP', interval=1, alpha=0.3):
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=False)

    # Selecting which model data we want to plot based on the diffusion type
    kukulka, kpp, _ = utils_v.boolean_diff_type(diffusion_type)
    # Selecting which model data we want to plot based on the diffusion scheme
    boundary_list = ['Reflect', 'Ceiling', 'Mixed']
    linestyle = ['-', '--', ':']

    # Creating the axis
    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size)

    mean_wind = np.mean(utils.utils_physics.beaufort_limits()[beaufort])
    # _, _ = utils_v.add_observations(ax, norm_depth=False, alpha=alpha, wind_range=utils.utils_physics.beaufort_limits()[4])

    for count, boundary in enumerate(boundary_list):
        # Plotting the distribution according to the Kukulka parametrization
        if kukulka:
            profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                          output_step=output_step, diffusion_type='Kukulka',
                                                          boundary=boundary, alpha_list=alpha_list)
            for counter in range(len(profile_dict['concentration_list'])):
                w_rise = w_rise_list[counter]
                ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                        label=label_boundary(w_10=mean_wind, w_rise=w_rise, diffusion_type='Kukulka', boundary=boundary),
                        linestyle=linestyle[count], color=utils_v.return_color(counter))

        # Plotting the distribution according to the KPP parametrization
        if kpp:
            profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                          output_step=output_step, diffusion_type='KPP',
                                                          boundary=boundary, alpha_list=alpha_list)
            for counter in range(len(profile_dict['concentration_list'])):
                w_rise = w_rise_list[counter]
                ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                        label=label_boundary(w_10=mean_wind, w_rise=w_rise, diffusion_type='KPP', boundary=boundary),
                        linestyle=linestyle[count], color=utils_v.return_color(counter))

    lines, labels = ax.get_legend_handles_labels()
    # Adding the legend
    ax.legend(lines, labels, fontsize=legend_size, loc='lower right')
    ax.set_title(r'w$_{10}$ = 5.4-7.9 m s$^{-1}$ - Boundary Conditions', fontsize=ax_label_size)

    # Saving the figure
    plt.savefig(saving_filename_boundary(settings.figure_dir + '/Boundary Conditions/', selection, close_up,
                                         diffusion_type), bbox_inches='tight', dpi=600)


def saving_filename_boundary(save_location, selection, close_up, diffusion_type):
    if close_up is None:
        return save_location + diffusion_type + '_boundary_full_variable={}.png'.format(selection)
    else:
        ymax, ymin = close_up
        return save_location + diffusion_type + '_boundary_max={}_min={}_variable={}.png'.format(ymax, ymin, selection)


def label_boundary(w_10, w_rise, diffusion_type, boundary):
    w_rise=np.abs(w_rise)
    boundary_dict = {'Mixed': 'M-0 Mixed', 'Reflect': 'M-0 Reflect', 'Reduce_dt': 'M-0 Reduce dt',
                     'Mixed_Markov': 'M-1 Mixed', 'Reflect_Markov':'M-1 Reflect',
                     'Reduce_dt_Markov':'M-1 Reduce dt', 'Ceiling': 'M-0 Ceiling'}
    if diffusion_type == 'Kukulka':
        return r'SWB, w$_{rise}$ '+'= {}'.format(w_rise) + ' m s$^{-1}$, ' + boundary_dict[boundary]
    elif diffusion_type == 'KPP':
        return r'KPP, w$_{rise}$ '+'= {}'.format(w_rise) + 'm s$^{-1}$, ' + boundary_dict[boundary]