import matplotlib.pyplot as plt
import utils, settings
from visualization import utils_visualization as utils_v
import numpy as np


def plot_model_field_data_comparison(w_10_list, w_rise_list, alpha_list, selection='w_10', output_step=-1,
                                     single_select=1,
                                     norm_depth=False, wind_sort=True, y_label='Depth (m)', close_up=None,
                                     x_label=r'Normalised Concentrations', fig_size=(16, 8), ax_label_size=16,
                                     legend_size=12, diffusion_type='SWB', boundary='Ceiling', alpha=0.3,
                                     beaufort=4):
    if norm_depth:
        y_label = 'Depth/MLD'
        correction = settings.MLD
    else:
        correction = 1.0
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=norm_depth)

    # Selecting which model data we want to plot based on the diffusion type
    swb, kpp, artificial = utils_v.boolean_diff_type(diffusion_type)
    # Selecting which model data we want to plot based on the diffusion scheme
    if boundary is 'Ceiling' or 'Ceiling_Markov':
        boundary_list = [boundary]
    if boundary is 'all':
        boundary_list = ['Ceiling', 'Ceiling_Markov']

    # Plotting data points for just one set of wind conditions
    if not wind_sort:
        wind_range = utils.beaufort_limits()[beaufort]

        # Get the base figure axis
        ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size)

        # Plotting the field data points
        _, _ = utils_v.add_observations(ax, norm_depth=norm_depth, alpha=alpha, wind_range=wind_range)

        mean_wind = np.mean(utils.beaufort_limits()[beaufort])
        for boundary in boundary_list:
            # Plotting the distribution according to the SWB parametrization
            if swb:
                profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                              output_step=output_step, diffusion_type='SWB',
                                                              boundary=boundary, alpha_list=alpha_list)
                for counter in range(len(profile_dict['concentration_list'])):
                    ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'] / correction,
                            label=utils_v.label_SWB(selection=selection,
                                                    parameters=profile_dict['parameter_concentrations'][counter]),
                            linestyle='-', color=utils_v.return_color(counter))

            # Plotting the distribution according to the KPP parametrization
            if kpp:
                profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                              output_step=output_step, diffusion_type='KPP',
                                                              boundary=boundary, alpha_list=alpha_list)
                for counter in range(len(profile_dict['concentration_list'])):
                    ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'] / correction,
                            label=utils_v.label_KPP(selection=selection,
                                                    parameters=profile_dict['parameter_concentrations'][counter]),
                            linestyle='--', color=utils_v.return_color(counter))
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
        ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size, shape=(2, 3), plot_num=plot_num,
                                 legend_axis=True)

        beaufort = utils.beaufort_limits()
        # The figures for the different wind conditions, where we only want to have this for the higher rise velocities
        w_rise_list = utils_v.rise_velocity_selector(size_class='large', w_rise_list=w_rise_list)
        for scale in range(plot_num - 1):
            _, _ = utils_v.add_observations(ax[scale], norm_depth=norm_depth, wind_range=beaufort[scale + 1],
                                            alpha=alpha)
            ax[scale].set_title(sub_titles[scale], fontsize=ax_label_size)
            for boundary in boundary_list:
                if swb:
                    profile_dict = utils_v.get_concentration_list(w_10_list, w_rise_list, 'all', single_select,
                                                                  output_step=output_step, diffusion_type='SWB',
                                                                  boundary=boundary, alpha_list=alpha_list)
                    for model in range(len(w_rise_list)):
                        _, w_rise = profile_dict['parameter_concentrations'][scale + model * len(w_10_list)]
                        concentration = profile_dict['concentration_list'][scale + model * len(w_10_list)]
                        ax[scale].plot(concentration, profile_dict['depth_bins'] / correction,
                                       label=label_model_field_comparison(w_rise=w_rise, diffusion_type='SWB', boundary=boundary),
                                       linestyle=utils_v.determine_linestyle(boundary, boundary_list, kpp, swb,
                                                                             'SWB'),
                                       color=utils_v.return_color(model))
                if kpp:
                    profile_dict = utils_v.get_concentration_list(w_10_list, w_rise_list, 'all', single_select,
                                                                  output_step=output_step, diffusion_type='KPP',
                                                                  boundary=boundary, alpha_list=alpha_list)
                    for model in range(len(w_rise_list)):
                        _, w_rise = profile_dict['parameter_concentrations'][scale + model * len(w_10_list)]
                        concentration = profile_dict['concentration_list'][scale + model * len(w_10_list)]
                        ax[scale].plot(concentration, profile_dict['depth_bins'] / correction,
                                       label=label_model_field_comparison(w_rise=w_rise, diffusion_type='KPP', boundary=boundary),
                                       linestyle=utils_v.determine_linestyle(boundary, boundary_list, kpp, swb,
                                                                             'KPP'),
                                       color=utils_v.return_color(model))
            lines, labels = ax[scale].get_legend_handles_labels()
        # Now, a final plot showing just the w_r=0.00003 case
        scale = plot_num - 1
        ax[scale].set_title(sub_titles[scale], fontsize=ax_label_size)
        for boundary in boundary_list:
            if swb:
                profile_dict = utils_v.get_concentration_list(w_10_list, [-0.0003], 'w_10', single_select=0,
                                                              output_step=output_step, diffusion_type='SWB',
                                                              boundary=boundary, alpha_list=alpha_list)
                for model in range(len(w_10_list)):
                    w_10, w_rise = profile_dict['parameter_concentrations'][model]
                    concentration = profile_dict['concentration_list'][model]
                    ax[scale].plot(concentration, profile_dict['depth_bins'] / correction,
                                   label=label_model_field_comparison(w_10=w_10, diffusion_type='SWB', boundary=boundary),
                                   linestyle=utils_v.determine_linestyle(boundary, boundary_list, kpp, swb,
                                                                         'SWB'),
                                   color=utils_v.discrete_color_from_cmap(model, len(w_10_list)))
            if kpp:
                profile_dict = utils_v.get_concentration_list(w_10_list, [-0.0003], 'w_10', single_select=0,
                                                              output_step=output_step, diffusion_type='KPP',
                                                              boundary=boundary, alpha_list=alpha_list)
                for model in range(len(w_10_list)):
                    w_10, w_rise = profile_dict['parameter_concentrations'][model]
                    concentration = profile_dict['concentration_list'][model]
                    ax[scale].plot(concentration, profile_dict['depth_bins'] / correction,
                                   label=label_model_field_comparison(w_10=w_10, diffusion_type='KPP', boundary=boundary),
                                   linestyle=utils_v.determine_linestyle(boundary, boundary_list, kpp, swb,
                                                                         'KPP'),
                                   color=utils_v.discrete_color_from_cmap(model, len(w_10_list)))
            lines_w_10, labels_w_10 = ax[scale].get_legend_handles_labels()
        # Adding the legend for the split wind plots
        ax[-2].legend(lines, labels, fontsize=legend_size, loc='upper left')
        ax[-2].axis('off')
        # Adding the legend for the w_r=0.0003 case
        ax[-1].legend(lines_w_10, labels_w_10, fontsize=legend_size, loc='upper left')
        ax[-1].axis('off')

    # Saving the figure
    save_name = model_field_data_comparison_name(diffusion_type, boundary, alpha_list, close_up=close_up,
                                                 wind_sort=wind_sort, norm_depth=norm_depth, beaufort=beaufort)
    plt.savefig(save_name, bbox_inches='tight')


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
                                     norm_depth=False, output_type='.png', beaufort=1):
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
        figure_name += '_Bft{}_'.format(beaufort)
    if norm_depth:
        figure_name += '_normalized_depth'
    return figure_name + '_dt={}'.format(settings.dt_int.seconds) + output_type