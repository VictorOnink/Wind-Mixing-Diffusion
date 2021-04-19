import matplotlib.pyplot as plt

import settings
import utils
import utils.utils_physics
from visualization import utils_visualization as utils_v


def plot_field_data_overview(norm_depth=False, wind_sort=False, y_label='Depth (m)', close_up=None,
                             x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(16, 8), ax_label_size=16,
                             legend_size=12, alpha=0.3):
    if norm_depth:
        y_label = 'Depth/MLD'
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=norm_depth)

    # Plotting all data points, with no sorting based on wind conditions
    if not wind_sort:
        # Get the base figure axis
        ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size)

        # Plotting the data points
        legend_line, legend_label = utils_v.add_observations(ax, norm_depth=norm_depth, alpha=alpha)

        # Adding the legend
        ax.legend(legend_line, legend_label, fontsize=legend_size, loc='lower right')

    # Plotting data points, split over multiple plots according to Beaufort wind scale
    if wind_sort:
        # Titles for the subplots
        sub_titles = [r'(a) u$_{10}$=0.2-1.5 m s$^{-1}$', r'(b) u$_{10}$=1.5-3.3 m s$^{-1}$',
                      r'(c) u$_{10}$=3.3-5.4 m s$^{-1}$', r'(d) u$_{10}$=5.4-7.9 m s$^{-1}$',
                      r'(e) u$_{10}$=7.9-10.7 m s$^{-1}$',
                      r' ']
        # Get the base figure axis
        plot_num = 6
        ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size, shape=(2, 3), plot_num=plot_num)
        beaufort = utils.utils_physics.beaufort_limits()
        for scale in range(plot_num):
            line, label = utils_v.add_observations(ax[scale], norm_depth=norm_depth, wind_range=beaufort[scale + 1],
                                                   alpha=alpha)
            ax[scale].set_title(sub_titles[scale], fontsize=ax_label_size)

        # Adding the legend
        ax[-1].legend(line, label, fontsize=legend_size, loc='lower right')
        ax[-1].axis('off')

    plt.savefig(field_data_figure_names(close_up, wind_sort, norm_depth), bbox_inches='tight')


def field_data_figure_names(close_up=None, wind_sort=False, norm_depth=False, output_type='.png'):
    figure_name = settings.figure_dir + '/Field Data/field_data_'
    if close_up is not None:
        max, min = close_up
        figure_name += 'max_{}_min_{}'.format(max, min)
    if wind_sort:
        figure_name += '_wind_sort'
    if norm_depth:
        figure_name += '_normalized_depth'
    return figure_name + output_type