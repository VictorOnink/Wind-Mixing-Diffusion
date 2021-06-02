import settings
import matplotlib.pyplot as plt
import utils
import visualization.utils_visualization
from visualization import utils_visualization as utils_v


def timestep_comparison(w_10_list, w_rise_list, alpha_list, selection='w_10', close_up=None,
                        y_label='Depth (m)', x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(8, 8),
                        ax_label_size=16, legend_size=12, single_select=0, mld=settings.MLD,
                        diffusion_type='SWB', interval=1, boundary='Mixed', diffusion_curve=True):
    """
    The timestep here refers to the concentration profile over time in a simulation
    :param w_10_list: list of wind speeds
    :param w_rise_list: list of rise velocities
    :param alpha_list: list of alpha values
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
    :param diffusion_curve: if True, add the diffusion profile to the figure
    :return:
    """
    # Load the concentration profiles
    profile_dict = utils_v.get_concentration_list(w_10_list, w_rise_list, selection, single_select, diffusion_type,
                                                  all_timesteps=True, boundary=boundary, mld=mld, alpha_list=alpha_list)

    # Creating the axis
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=False)
    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size)

    # Plotting the distribution according to the SWB parametrization
    if diffusion_type is 'SWB':
        for counter in range(0, len(profile_dict['concentration_list']), interval):
            ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                    label=label_time_step(counter, interval),
                    linestyle='-', color=visualization.utils_visualization.return_color(counter))
    # Plotting the distribution according to the KPP parametrization
    if diffusion_type is 'KPP':
        for counter in range(0, len(profile_dict['concentration_list']), interval):
            ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                    label=label_time_step(counter, interval),
                    linestyle='-', color=visualization.utils_visualization.return_color(counter))

    lines, labels = ax.get_legend_handles_labels()

    # Plotting the diffusion curve
    if diffusion_curve:
        w_10, w_rise = profile_dict['parameter_concentrations'][0]
        ax2 = utils_v.diffusion_curve_axis(ax, ax_label_size, w_10, w_rise, profile_dict, diffusion_type, 'black')
        lines2, labels2 = ax2.get_legend_handles_labels()
        lines += lines2
        labels += labels2

    # Adding the legend
    ax.legend(lines, labels, fontsize=legend_size, loc='lower right')

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
        return save_location + diffusion_type + '_time_step_full.png'
    else:
        ymax, ymin = close_up
        return save_location + diffusion_type + '_time_step_max={}_min={}.png'.format(ymax, ymin)