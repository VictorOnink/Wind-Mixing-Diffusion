import settings
import matplotlib.pyplot as plt
import utils
import visualization.utils_visualization
from visualization import utils_visualization as utils_v


def timestep_comparison(w_10_list, w_rise_list, alpha_list, selection='w_10', close_up=None,
                        y_label='Depth (m)', x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(8, 8),
                        ax_label_size=16, legend_size=12, single_select=0, mld=settings.MLD,
                        diffusion_type='Kukulka', interval=1, boundary='Mixed', diffusion_curve=True):
    # Load the relevant data for the figure
    profile_dict = utils_v.get_concentration_list(w_10_list, w_rise_list, selection, single_select, diffusion_type,
                                                  all_timesteps=True, boundary=boundary, mld=mld, alpha_list=alpha_list)
    # Preparing for the actual plotting
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=False)

    # Creating the axis
    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size)

    # Plotting the distribution according to the Kukulka parametrization
    if diffusion_type is 'Kukulka':
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
    t = steps * interval * settings.dt_out.seconds // 3600
    return 't = {} hours'.format(t)


def saving_filename_time_step(save_location, close_up, diffusion_type):
    if close_up is None:
        return save_location + diffusion_type + '_time_step_full.png'
    else:
        ymax, ymin = close_up
        return save_location + diffusion_type + '_time_step_max={}_min={}.png'.format(ymax, ymin)