import settings
import matplotlib.pyplot as plt
from visualization import utils_visualization as utils_v
import numpy as np
import utils


def timestep_comparison(selection='w_10', close_up=None,
                        y_label='Depth (m)', x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(12, 8),
                        ax_label_size=14, legend_size=10, single_select=0, mld=settings.MLD,
                        diffusion_type='KPP', interval=1, boundary='Ceiling', wave_roughness=False):
    """
    The timestep here refers to the concentration profile over time in a simulation
    :param selection: selection criteria for loading the parcels concentration profiles
    :param close_up: setting the axis limits of the Y axis as (max, min)
    :param y_label: label of the y axis
    :param x_label: label of the x axis
    :param fig_size: size of the figure
    :param ax_label_size: fontsize of the axis labels
    :param legend_size: fontsize of hte legend
    :param single_select: selection index related to 'selection'
    :param mld: mixed layer depth
    :param diffusion_type: type of diffusion, KPP or SWB
    :param interval: interval of the profiles we plot in time (1 = every output timestep, 2 = every second output
                     timestep, etc.)
    :param boundary: which boundary condition, and M-0 or M-1
    :param wave_roughness: if True, have surface roughness be wave height dependent
    :return:
    """
    # Load the concentration profiles


    # Creating the axis
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=False)
    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size, shape=(1, 2), plot_num=2,
                             all_x_labels=True, legend_axis=True, width_ratios=[1, 1, 0.3])

    # Loading the field data for Beaufort 4 wind conditions
    # _, _ = utils_v.add_observations(ax, norm_depth=False, wind_range=utils.beaufort_limits()[4])

    # Plotting the distribution according to the KPP parametrization
    profile_dict = utils_v.get_concentration_list([6.65], [-0.003], selection, single_select, 'KPP',
                                                  all_timesteps=True, boundary=boundary, mld=mld, alpha_list=[0],
                                                  wave_roughness=wave_roughness)
    subdivision = len(profile_dict['concentration_list']) // interval
    for counter in range(0, len(profile_dict['concentration_list']), interval):
        ax[0].plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                  label=label_time_step(counter, interval),
                  linestyle='-', color=utils_v.discrete_color_from_cmap(index=counter, subdivisions=subdivision))
    ax[0].set_title('(a) KPP', fontsize=ax_label_size)

    # Plotting the distribution according to the SWB parametrization
    profile_dict = utils_v.get_concentration_list([6.65], [-0.003], selection, single_select, 'SWB',
                                                  all_timesteps=True, boundary=boundary, mld=mld, alpha_list=[0],
                                                  wave_roughness=wave_roughness)
    for counter in range(0, len(profile_dict['concentration_list']), interval):
        ax[1].plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                  label=label_time_step(counter, interval),
                  linestyle='-', color=utils_v.discrete_color_from_cmap(index=counter, subdivisions=subdivision))
    ax[1].set_title('(b) SWB', fontsize=ax_label_size)

    lines, labels = ax[0].get_legend_handles_labels()

    # Adding the legend
    ax[-1].legend(lines, labels, fontsize=legend_size, loc='upper right')

    # Saving the figure
    plt.savefig(saving_filename_time_step(settings.figure_dir, close_up, diffusion_type),
                bbox_inches='tight', dpi=600)


def label_time_step(steps, interval):
    """ the labels of hte various lines"""
    t = steps * interval * settings.dt_out.seconds // 3600
    return 't = {} hours'.format(t)


def saving_filename_time_step(save_location, close_up, diffusion_type):
    """" The filename of hte figure """
    if close_up is None:
        return save_location + 'time_step_comparison_full.png'
    else:
        ymax, ymin = close_up
        return save_location + 'time_step_comparison_max={}_min={}.png'.format(ymax, ymin)