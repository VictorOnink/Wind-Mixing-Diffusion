import settings
import matplotlib.pyplot as plt
import utils
import utils.utils_physics
import visualization.utils_visualization
from visualization import utils_visualization as utils_v
import numpy as np


def eulerian_field_data_comparison(w_10_list, w_rise_list, alpha_list, selection='w_10', output_step=-1,
                                   single_select=1,
                                   norm_depth=False, wind_sort=True, y_label='Depth (m)', close_up=None,
                                   x_label=r'Normalised Plastic Concentrations ($C/C_0$)', fig_size=(16, 8),
                                   ax_label_size=16,
                                   legend_size=10, diffusion_type='Kukulka', boundary='Reflect_Markov', alpha=0.3,
                                   beaufort=4):
    if norm_depth:
        y_label = 'Depth/MLD'
        correction = settings.MLD
    else:
        correction = 1.0
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=norm_depth)

    # Selecting which model data we want to plot based on the diffusion type
    kukulka, kpp, artificial = utils_v.boolean_diff_type(diffusion_type)

    # Plotting data points for just one set of wind conditions
    if not wind_sort:
        wind_range = utils.utils_physics.beaufort_limits()[beaufort]

        # Get the base figure axis
        ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size)

        # Plotting the field data points
        _, _ = utils_v.add_observations(ax, norm_depth=norm_depth, alpha=alpha, wind_range=wind_range)

        mean_wind = np.mean(utils.utils_physics.beaufort_limits()[beaufort])

        if kukulka:
            for counter, w_r in enumerate(np.abs(w_rise_list)):
                eul_dict = utils.load_obj(
                    utils.get_eulerian_output_name(w_10=mean_wind, w_rise=w_r, diffusion_type='Kukulka'))
                ax.plot(eul_dict['C'], eul_dict['Z'] / correction, linestyle='--',
                        color=visualization.utils_visualization.return_color(counter),
                        label=line_labels(w_rise=w_r, boundary_type='eulerian', diffusion_type='Kukulka'))

        # Plotting the distribution according to the KPP parametrization
        if kpp:
            for counter, w_r in enumerate(np.abs(w_rise_list)):
                eul_dict = utils.load_obj(
                    utils.get_eulerian_output_name(w_10=mean_wind, w_rise=w_r, diffusion_type='KPP'))
                ax.plot(eul_dict['C'], eul_dict['Z'] / correction, linestyle='--',
                        color=visualization.utils_visualization.return_color(counter),
                        label=line_labels(w_rise=w_r, boundary_type='eulerian', diffusion_type='KPP'))
        lines, labels = ax.get_legend_handles_labels()

        # Adding the legend
        ax.legend(fontsize=legend_size, loc='lower right')

    # Plotting data points, split over multiple plots according to Beaufort wind scale
    if wind_sort:
        # Titles for the subplots
        sub_titles = [r'(a) u$_{10}$=0.2-1.5 m s$^{-1}$', r'(b) u$_{10}$=1.5-3.3 m s$^{-1}$',
                      r'(c) u$_{10}$=3.3-5.4 m s$^{-1}$', r'(d) u$_{10}$=5.4-7.9 m s$^{-1}$',
                      r'(e) u$_{10}$=7.9-10.7 m s$^{-1}$', r'(f) w$_r$=0.0003 m s$^{-1}$']
        # Get the base figure axis
        plot_num = 6
        ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size, shape=(2, 3),
                                 plot_num=plot_num,
                                 legend_axis=True)

        beaufort = utils.utils_physics.beaufort_limits()
        # The figures for the different wind conditions, where we only want to have this for the higher rise velocities
        w_rise_list = utils_v.rise_velocity_selector(size_class='large', w_rise=w_rise_list)
        for scale in range(plot_num - 1):
            _, _ = utils_v.add_observations(ax[scale], norm_depth=norm_depth, wind_range=beaufort[scale + 1],
                                            alpha=alpha)
            mean_wind = np.mean(beaufort[scale + 1])
            ax[scale].set_title(sub_titles[scale], fontsize=ax_label_size)
            if kukulka:
                for counter, w_r in enumerate(np.abs(w_rise_list)):
                    eul_dict = utils.load_obj(
                        utils.get_eulerian_output_name(w_10=mean_wind, w_rise=w_r, diffusion_type='Kukulka'))
                    ax[scale].plot(eul_dict['C'], eul_dict['Z'] / correction, linestyle='-',
                            color=visualization.utils_visualization.return_color(counter),
                            label=line_labels(w_rise=w_r, boundary_type='eulerian', diffusion_type='Kukulka'))
            if kpp:
                for counter, w_r in enumerate(np.abs(w_rise_list)):
                    eul_dict = utils.load_obj(
                        utils.get_eulerian_output_name(w_10=mean_wind, w_rise=w_r, diffusion_type='KPP'))
                    ax[scale].plot(eul_dict['C'], eul_dict['Z'] / correction, linestyle='--',
                            color=visualization.utils_visualization.return_color(counter),
                            label=line_labels(w_rise=w_r, boundary_type='eulerian', diffusion_type='KPP'))
            lines, labels = ax[scale].get_legend_handles_labels()
        # Now, a final plot showing just the w_r=0.00003 case
        scale = plot_num - 1
        ax[scale].set_title(sub_titles[scale], fontsize=ax_label_size)
        if kukulka:
            for counter, w_10 in enumerate(np.abs(w_10_list)):
                eul_dict = utils.load_obj(
                    utils.get_eulerian_output_name(w_10=w_10, w_rise=0.0003, diffusion_type='Kukulka'))
                ax[scale].plot(eul_dict['C'], eul_dict['Z'] / correction, linestyle='-',
                               color=utils_v.discrete_color_from_cmap(counter, len(w_10_list)),
                               label=line_labels(w_rise=w_r, boundary_type='eulerian', diffusion_type='Kukulka'))
        if kpp:
            for counter, w_10 in enumerate(np.abs(w_10_list)):
                eul_dict = utils.load_obj(
                    utils.get_eulerian_output_name(w_10=w_10, w_rise=0.0003, diffusion_type='KPP'))
                ax[scale].plot(eul_dict['C'], eul_dict['Z'] / correction, linestyle='--',
                               color=utils_v.discrete_color_from_cmap(counter, len(w_10_list)),
                               label=line_labels(w_rise=w_r, boundary_type='eulerian', diffusion_type='KPP'))
        lines_w_10, labels_w_10 = ax[scale].get_legend_handles_labels()
        # Adding the legend for the split wind plots
        ax[-2].legend(lines, labels, fontsize=legend_size, loc='upper left')
        ax[-2].axis('off')
        # Adding the legend for the w_r=0.0003 case
        ax[-1].legend(lines_w_10, labels_w_10, fontsize=legend_size, loc='upper left')
        ax[-1].axis('off')

    # Saving the figure
    save_name = model_euler_field_data_comparison_name(diffusion_type, boundary, alpha_list, close_up=close_up,
                                                 wind_sort=wind_sort, norm_depth=norm_depth, beaufort=beaufort)
    plt.savefig(save_name, bbox_inches='tight')


def line_labels(w_rise, boundary_type, diffusion_type):
    w_rise = np.abs(w_rise)
    boundary_dict = {'Reflect': 'M-0', 'Reflect_Markov': 'M-1', 'eulerian': 'Eulerian'}
    boundary = boundary_dict[boundary_type]
    diffusion_dict = {'KPP': 'KPP', 'Kukulka': 'SWB'}
    diffusion = diffusion_dict[diffusion_type]
    filename = '{}, {}, '.format(diffusion, boundary) + r'$w_r=$' + '{}'.format(w_rise) + r' m s$^{-1}$'
    return filename


def model_euler_field_data_comparison_name(diffusion_type, boundary, alpha_list,close_up=None, wind_sort=False,
                                           norm_depth=False, output_type='.png', beaufort=1):
    diff_dict = {'Kukulka': 'Kukulka', 'KPP': 'KPP', 'all': 'Kukulka_KPP'}
    figure_name = settings.figure_dir + 'Eulerian_comparison/euler_field_data_{}'.format(diff_dict[diffusion_type])
    if close_up is not None:
        max, min = close_up
        figure_name += '_max_{}_min_{}'.format(max, min)
    if wind_sort:
        figure_name += '_wind_sort'
    elif not wind_sort:
        figure_name += '_Bft{}_'.format(beaufort)
    if norm_depth:
        figure_name += '_normalized_depth'
    return figure_name + output_type