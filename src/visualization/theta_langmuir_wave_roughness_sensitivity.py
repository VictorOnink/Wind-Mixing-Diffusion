import matplotlib.pyplot as plt
import utils
import settings
from visualization import utils_visualization as utils_v
import numpy as np


def theta_langmuir_wave_roughness_sensitivity(w_10_list, w_rise, theta_list, output_step=-1,
                                              single_select=1, add_variability=True, y_label='Depth (m)', close_up=None,
                                              x_label=r'Normalised Concentrations', fig_size=(16, 8), ax_label_size=16,
                                              legend_size=12, boundary='Ceiling', alpha=0.15, with_observations=True):
    # Setting the axes ranges
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=False)

    # Titles for the subplots
    sub_titles = [r'(a) u$_{10}$=0.2-1.5 m s$^{-1}$', r'(b) u$_{10}$=1.5-3.3 m s$^{-1}$',
                  r'(c) u$_{10}$=3.3-5.4 m s$^{-1}$', r'(d) u$_{10}$=5.4-7.9 m s$^{-1}$',
                  r'(e) u$_{10}$=7.9-10.7 m s$^{-1}$']

    # Get the base figure axis
    plot_num = 6
    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size, shape=(2, 3), plot_num=plot_num,
                             legend_axis=False)

    # set subtitles
    for scale in range(sub_titles.__len__()):
        ax[scale].set_title(sub_titles[scale], fontsize=ax_label_size)

    # Load wind conditions
    beaufort = utils.beaufort_limits()

    line_style = {False: '-', True: 'dotted'}

    for scale in range(w_10_list.__len__()):
        if with_observations:
            _, _ = utils_v.add_observations(ax[scale], norm_depth=False, wind_range=beaufort[scale + 1],
                                            alpha=alpha, mean_concentrations=True)
        for wave_roughness in [True, False]:
            for index_theta, theta in enumerate(theta_list):
                plot_color = utils_v.discrete_color_from_cmap(index=index_theta, subdivisions=theta_list.__len__())
                profile_dict = utils_v.get_concentration_list(w_10_list, [w_rise], 'all', single_select,
                                                              output_step=output_step, diffusion_type='KPP',
                                                              boundary=boundary, alpha_list=[0], theta=theta,
                                                              wave_roughness=wave_roughness)
                _, w_rise = profile_dict['parameter_concentrations'][scale]
                concentration = profile_dict['concentration_list'][scale]
                depth = profile_dict['depth_bins']
                ax[scale].plot(concentration, depth, label=label_theta(theta, wave_roughness),
                               linestyle=line_style[wave_roughness], color=plot_color)
                if add_variability:
                    std = profile_dict['std_list'][scale]
                    upper_limit, lower_limit = concentration + std, concentration - std
                    ax[scale].fill_betweenx(depth, lower_limit, upper_limit, alpha=0.2, color=plot_color)

        lines, labels = ax[scale].get_legend_handles_labels()

    # Labels for field data
    if with_observations:
        field_lines, field_labels = lines[-6:], labels[-6:]
    # Labels for wind conditions
    wind_lines = [plt.plot([], [], c=utils_v.discrete_color_from_cmap(index=index, subdivisions=theta_list.__len__()),
                           label=r'$\theta=$' + '{}'.format(theta), linestyle='-')[0]
                  for index, theta in enumerate(theta_list)]
    labels, styles = [r'$z_0=$ Zhao & Li (2019)', r'$z_0=0.1\times H_s$'], ['-', 'dotted']
    roughness_lines = [plt.plot([], [], c='k', label=label, linestyle=line)[0] for label, line in zip(labels, styles)]

    # Adding the legend
    if with_observations:
        handles = field_lines + wind_lines + roughness_lines
    else:
        handles = wind_lines + roughness_lines
    ax[-1].legend(handles=handles, fontsize=legend_size, loc='upper left')
    ax[-1].axis('off')

    # Saving the figure
    save_name = settings.figure_dir + 'KPP_theta_influence_w_r={}.png'.format(w_rise)
    plt.savefig(save_name, bbox_inches='tight')


def label_theta(theta, wave_roughness):
    label = r'$\theta=$' + '{}'.format(theta)
    if wave_roughness:
        label += r', $z_0=0.1\times H_S$'
    return label


def label_model_field_comparison(w_rise=None, w_10=None, diffusion_type=None, boundary=None):
    boundary_dict = {'Ceiling': 'M-0', 'Ceiling_Markov': 'M-1'}
    if diffusion_type is 'SWB':
        diff = 'SWB'
    elif diffusion_type is 'KPP':
        diff = 'KPP'
    if w_10 is None:
        w_rise = np.abs(w_rise)
        return diff + ', {}'.format(boundary_dict[boundary]) + ', w$_{rise}$ ' + '= {} m s'.format(w_rise) + r'$^{-1}$'
    else:
        return diff + ', {}'.format(boundary_dict[boundary]) + ', u$_{10}$ ' + '= {:.2f} m s'.format(w_10) + r'$^{-1}$'


def model_field_data_comparison_name(diffusion_type, boundary, alpha_list,close_up=None, wind_sort=False,
                                     norm_depth=False, output_type='.png', beaufort=1, theta=1):
    diff_dict = {'SWB': 'SWB', 'KPP': 'KPP', 'all': 'SWB_KPP'}
    figure_name = settings.figure_dir + 'model_field_data_{}_{}'.format(diff_dict[diffusion_type], boundary)
    if 'Markov' in boundary:
        figure_name += '_alpha_{}_'.format(alpha_list[0])
    if close_up is not None:
        max, min = close_up
        figure_name += '_max_{}_min_{}'.format(max, min)
    if wind_sort:
        figure_name += '_wind_sort'
    elif not wind_sort:
        figure_name += '_Bft{}'.format(beaufort)
    if norm_depth:
        figure_name += '_normalized_depth'
    figure_name += '_theta={}'.format(theta)
    return figure_name + '_dt={}'.format(settings.dt_int.seconds) + output_type


