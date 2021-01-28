import settings
import matplotlib.pyplot as plt
import utils
import utils_visualization as utils_v
import numpy as np


def basic_profile_figure(w_10_list, w_rise_list, selection='w_10', close_up=None,
                         y_label='Depth (m)', x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(8, 8),
                         ax_label_size=16, legend_size=12, single_select=1,
                         output_step=-1, diffusion_type='Kukulka', boundary='Mixed', diffusion_curve=True):
    # Preparing for the actual plotting
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=False)
    # Creating the axis
    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size)

    # Plotting the distribution according to the Kukulka parametrization
    if diffusion_type is 'Kukulka':
        profile_dict = utils_v.get_concentration_list(w_10_list, w_rise_list, selection, single_select,
                                                      output_step=output_step, diffusion_type=diffusion_type,
                                                      boundary=boundary)
        for counter in range(len(profile_dict['concentration_list'])):
            ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                    label=utils_v.label_kukulka(selection=selection,
                                                parameters=profile_dict['parameter_kukulka'][counter]),
                    linestyle='-', color=utils.return_color(counter))

    # Plotting the distribution according to the KPP parametrization
    if diffusion_type is 'KPP':
        profile_dict = utils_v.get_concentration_list(w_10_list, w_rise_list, selection, single_select,
                                                      output_step=output_step, diffusion_type=diffusion_type,
                                                      boundary=boundary)
        for counter in range(len(profile_dict['concentration_list'])):
            ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                    label=utils_v.label_KPP(parameters=profile_dict['parameter_kukulka'][counter], selection=selection),
                    linestyle='-', color=utils.return_color(counter))
    lines, labels = ax.get_legend_handles_labels()

    # Plotting the diffusion curve
    kukulka, kpp = utils_v.boolean_diff_type(diffusion_type)
    if diffusion_curve:
        if kukulka:
            profile_dict = utils_v.get_concentration_list(w_10_list, w_rise_list, selection, single_select,
                                                          output_step=output_step, diffusion_type='Kukulka',
                                                          boundary=boundary)
            for counter in range(len(profile_dict['parameter_concentrations'])):
                w_10, w_rise = profile_dict['parameter_concentrations'][counter]
                ax2 = utils_v.diffusion_curve_axis(ax, ax_label_size, w_10, w_rise, profile_dict, 'Kukulka',
                                                   utils.return_color(counter))
                lines2, labels2 = ax2.get_legend_handles_labels()
                lines += lines2
                labels += labels2
        if kpp:
            profile_dict = utils_v.get_concentration_list(w_10_list, w_rise_list, selection, single_select,
                                                          output_step=output_step, diffusion_type='KPP',
                                                          boundary=boundary)
            if selection is 'w_rise':
                range_lim = 1
            if selection is 'w_10':
                range_lim = len(profile_dict['parameter_concentrations'])
            for counter in range(range_lim):
                w_10, w_rise = profile_dict['parameter_concentrations'][counter]
                ax2 = utils_v.diffusion_curve_axis(ax, ax_label_size, w_10, w_rise, profile_dict, "KPP",
                                                   'black')
                lines2, labels2 = ax2.get_legend_handles_labels()
                lines += lines2
                labels += labels2

    # Adding the legend
    ax.legend(lines, labels, fontsize=legend_size, loc='lower right')

    # Saving the figure
    plt.savefig(utils_v.saving_filename_basic_profile(settings.figure_dir, selection, close_up, diffusion_type),
                bbox_inches='tight', dpi=600)


def just_diffusion_profile(w_10_list, y_label='Depth (m)', x_label=r'$K_z$ ($10^{-2}$ m$^2$ s$^{-1}$)',
                           fig_size=(8, 8), ax_label_size=16, legend_size=12):
    ymax, ymin = 0, -1 * (settings.MLD + 10)
    ax_range = (2.5, 0, ymax, ymin)
    depth = np.linspace(ymax, np.abs(ymin), 1000)

    # Creating the axis
    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size)

    # Plotting the diffusion profile according to the Kukulka approach
    for count, w_10 in enumerate(w_10_list):
        profile = utils.get_vertical_diffusion_profile(w_10, depth, 'Kukulka')
        ax.plot(profile * 100, -1 * depth, color=utils.return_color(count), linestyle='-',
                label=utils_v.label_diffusivity_profile(w_10, 'Kukulka'))

    # Plotting the diffusion profile according the KPP approach
    for count, w_10 in enumerate(w_10_list):
        profile = utils.get_vertical_diffusion_profile(w_10, depth, 'KPP')
        ax.plot(profile * 100, -1 * depth, color=utils.return_color(count), linestyle='--',
                label=utils_v.label_diffusivity_profile(w_10, 'KPP'))

    # Adding the legend
    ax.legend(fontsize=legend_size, loc='lower right')

    plt.savefig(settings.figure_dir + 'diffusion_profile_MLD=20.png', bbox_inches='tight')


def timestep_comparison(w_10_list, w_rise_list, selection='w_10', close_up=None,
                        y_label='Depth (m)', x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(8, 8),
                        ax_label_size=16, legend_size=12, single_select=0, mld=settings.MLD,
                        diffusion_type='Kukulka', interval=1, boundary='Mixed', diffusion_curve=True):

    # Load the relevant data for the figure
    profile_dict = utils_v.get_concentration_list(w_10_list, w_rise_list, selection, single_select, diffusion_type,
                                                  all_timesteps=True, boundary=boundary, mld=mld)
    # Preparing for the actual plotting
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=False)

    # Creating the axis
    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size)

    # Plotting the distribution according to the Kukulka parametrization
    if diffusion_type is 'Kukulka':
        for counter in range(0, len(profile_dict['concentration_list']), interval):
            ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                    label=utils_v.label_time_step(counter, interval),
                    linestyle='-', color=utils.return_color(counter))
    # Plotting the distribution according to the KPP parametrization
    if diffusion_type is 'KPP':
        for counter in range(0, len(profile_dict['concentration_list']), interval):
            ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                    label=utils_v.label_time_step(counter, interval),
                    linestyle='-', color=utils.return_color(counter))

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
    plt.savefig(utils_v.saving_filename_time_step(settings.figure_dir, selection, close_up, diffusion_type),
                bbox_inches='tight', dpi=600)


def boundary_condition_comparison(w_rise_list, selection='w_10', close_up=None, output_step=-1,
                                  y_label='Depth (m)', x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(8, 8),
                                  ax_label_size=16, legend_size=12, single_select=0,
                                  diffusion_type='KPP', interval=1, alpha=0.3):
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=False)

    # Selecting which model data we want to plot based on the diffusion type
    kukulka, kpp = utils_v.boolean_diff_type(diffusion_type)
    # Selecting which model data we want to plot based on the diffusion scheme
    boundary_list = ['Mixed', 'Reflect', 'Reduce_dt', 'Mixed_Markov', 'Reflect_Markov', 'Reduce_dt_Markov']

    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size)
    line_style = ['-', '-', '-', '--', '--', '--']

    mean_wind = np.mean(utils.beaufort_limits()[4])
    _, _ = utils_v.add_observations(ax, norm_depth=False, alpha=alpha, wind_range=utils.beaufort_limits()[4])

    for count, boundary in enumerate(boundary_list):
        # Plotting the distribution according to the Kukulka parametrization
        if kukulka:
            profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                          output_step=output_step, diffusion_type='Kukulka',
                                                          boundary=boundary)
            for counter in range(len(profile_dict['concentration_list'])):
                ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                        label=utils_v.label_boundary(w_10=mean_wind, diffusion_type='KPP', boundary=boundary),
                        linestyle=line_style[count], color=utils.return_color(count % 3))

        # Plotting the distribution according to the KPP parametrization
        if kpp:
            profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                          output_step=output_step, diffusion_type='KPP',
                                                          boundary=boundary)
            for counter in range(len(profile_dict['concentration_list'])):
                ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                        label=utils_v.label_boundary(w_10=mean_wind, diffusion_type='KPP', boundary=boundary),
                        linestyle=line_style[count], color=utils.return_color(count % 3))

    lines, labels = ax.get_legend_handles_labels()
    # Adding the legend
    ax.legend(lines, labels, fontsize=legend_size, loc='lower right')
    ax.set_title(r'w$_{10}$ = 5.4-7.9 m s$^{-1}$ - Boundary Conditions', fontsize=ax_label_size)

    # Saving the figure
    plt.savefig(utils_v.saving_filename_boundary(settings.figure_dir + '/Boundary Conditions/', selection, close_up,
                                                 diffusion_type),
                bbox_inches='tight', dpi=600)


def plot_field_data_overview(norm_depth=False, wind_sort=False, y_label='Depth (m)', close_up=None,
                             x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(12, 8), ax_label_size=16,
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
        sub_titles = [r'(a) w$_{10}=0.2-1.5 m s$^{-1}$', r'(b) w$_{10}=1.5-3.3 m s$^{-1}$',
                      r'(c) w$_{10}=3.3-5.4 m s$^{-1}$', r'(d) w$_{10}=5.4-7.9 m s$^{-1}$',
                      r'(e) w$_{10}=7.9-10.7 m s$^{-1}$',
                      ' ']
        # Get the base figure axis
        plot_num = 6
        ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size, shape=(2, 3), plot_num=plot_num)
        beaufort = utils.beaufort_limits()
        for scale in range(plot_num):
            line, label = utils_v.add_observations(ax[scale], norm_depth=norm_depth, wind_range=beaufort[scale + 1],
                                                   alpha=alpha)
            ax[scale].set_title(sub_titles[scale], fontsize=ax_label_size)

        # Adding the legend
        ax[-1].legend(line, label, fontsize=legend_size, loc='lower right')
        ax[-1].axis('off')

    plt.savefig(utils_v.field_data_figure_names(close_up, wind_sort, norm_depth), bbox_inches='tight')


def plot_model_field_data_comparison(w_10_list, w_rise_list, selection='w_10', output_step=-1, single_select=1,
                                     norm_depth=False, wind_sort=True, y_label='Depth (m)', close_up=None,
                                     x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(16, 8), ax_label_size=16,
                                     legend_size=12, diffusion_type='Kukulka', boundary='Reflect_Markov', alpha=0.3,
                                     beaufort=4):
    if norm_depth:
        y_label = 'Depth/MLD'
        correction = settings.MLD
    else:
        correction = 1.0
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=norm_depth)

    # Selecting which model data we want to plot based on the diffusion type
    kukulka, kpp = utils_v.boolean_diff_type(diffusion_type)
    # Selecting which model data we want to plot based on the diffusion scheme
    if boundary is 'Reflect' or 'Reflect_Markov':
        boundary_list = [boundary]
    if boundary is 'all':
        boundary_list = ['Reflect', 'Reflect_Markov']

    # Plotting data points for just one set of wind conditions
    if not wind_sort:
        wind_range = utils.beaufort_limits()[beaufort]

        # Get the base figure axis
        ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size)

        # Plotting the field data points
        _, _ = utils_v.add_observations(ax, norm_depth=norm_depth, alpha=alpha, wind_range=wind_range)

        mean_wind = np.mean(utils.beaufort_limits()[beaufort])
        for boundary in boundary_list:
            # Plotting the distribution according to the Kukulka parametrization
            if kukulka:
                profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                              output_step=output_step, diffusion_type='Kukulka',
                                                              boundary=boundary)
                for counter in range(len(profile_dict['concentration_list'])):
                    ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'] / correction,
                            label=utils_v.label_kukulka(selection=selection,
                                                        parameters=profile_dict['parameter_kukulka'][counter]),
                            linestyle='-', color=utils.return_color(counter))

            # Plotting the distribution according to the KPP parametrization
            if kpp:
                profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                              output_step=output_step, diffusion_type='KPP',
                                                              boundary=boundary)
                for counter in range(len(profile_dict['concentration_list'])):
                    ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'] / correction,
                            label=utils_v.label_KPP(selection=selection,
                                                    parameters=profile_dict['parameter_kukulka'][counter]),
                            linestyle='--', color=utils.return_color(counter))
        lines, labels = ax.get_legend_handles_labels()

        # Adding the legend
        ax.legend(fontsize=legend_size, loc='lower right')

    # Plotting data points, split over multiple plots according to Beaufort wind scale
    if wind_sort:
        # Titles for the subplots
        sub_titles = [r'(a) w$_{10}=0.2-1.5 m s$^{-1}$', r'(b) w$_{10}=1.5-3.3 m s$^{-1}$',
                      r'(c) w$_{10}=3.3-5.4 m s$^{-1}$', r'(d) w$_{10}=5.4-7.9 m s$^{-1}$',
                      r'(e) w$_{10}=7.9-10.7 m s$^{-1}$',
                      ' ']
        # Get the base figure axis
        plot_num = 6
        ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size, shape=(2, 3), plot_num=plot_num)

        beaufort = utils.beaufort_limits()
        for scale in range(plot_num - 1):
            _, _ = utils_v.add_observations(ax[scale], norm_depth=norm_depth, wind_range=beaufort[scale + 1],
                                            alpha=alpha)
            ax[scale].set_title(sub_titles[scale], fontsize=ax_label_size)
            for boundary in boundary_list:
                if kukulka:
                    profile_dict = utils_v.get_concentration_list(w_10_list, w_rise_list, 'all', single_select,
                                                                  output_step=output_step, diffusion_type='Kukulka',
                                                                  boundary=boundary)
                    for model in range(len(w_rise_list)):
                        _, w_rise = profile_dict['parameter_kukulka'][scale + model * len(w_10_list)]
                        concentration = profile_dict['concentration_list'][scale + model * len(w_10_list)]
                        ax[scale].plot(concentration, profile_dict['depth_bins'] / correction,
                                       label=utils_v.label_model_field_comparison(w_rise, 'Kukulka', boundary),
                                       linestyle=utils_v.determine_linestyle(boundary, boundary_list, kpp, kukulka,
                                                                             'Kukulka'),
                                       color=utils.return_color(model))
                if kpp:
                    profile_dict = utils_v.get_concentration_list(w_10_list, w_rise_list, 'all', single_select,
                                                                  output_step=output_step, diffusion_type='KPP',
                                                                  boundary=boundary)
                    for model in range(len(w_rise_list)):
                        _, w_rise = profile_dict['parameter_kukulka'][scale + model * len(w_10_list)]
                        concentration = profile_dict['concentration_list'][scale + model * len(w_10_list)]
                        ax[scale].plot(concentration, profile_dict['depth_bins'] / correction,
                                       label=utils_v.label_model_field_comparison(w_rise, 'KPP', boundary),
                                       linestyle=utils_v.determine_linestyle(boundary, boundary_list, kpp, kukulka,
                                                                             'KPP'),
                                       color=utils.return_color(model))
            lines, labels = ax[scale].get_legend_handles_labels()

        # Adding the legend to the last plot, and hiding the grid
        ax[-1].legend(lines, labels, fontsize=legend_size, loc='lower right')
        ax[-1].axis('off')

    # Saving the figure
    save_name = utils_v.model_field_data_comparison_name(diffusion_type, boundary, close_up=close_up,
                                                         wind_sort=wind_sort, norm_depth=norm_depth, beaufort=beaufort)
    plt.savefig(save_name, bbox_inches='tight')


def mld_depth_influence(w_rise_list, MLD_list, selection='w_rise', output_step=-1, single_select=0,
                        y_label='Depth/MLD', close_up=None, beaufort=5,
                        x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(16, 8), ax_label_size=16,
                        legend_size=12, diffusion_type='KPP', boundary='Reflect', alpha=0.3):
    correction = settings.MLD
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=True)

    # Selecting which model data we want to plot based on the diffusion type
    kukulka, kpp = utils_v.boolean_diff_type(diffusion_type)

    # Get the base figure axis
    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size)

    # Plotting the field data points
    wind_range = utils.beaufort_limits()[beaufort]
    _, _ = utils_v.add_observations(ax, norm_depth=True, alpha=alpha, wind_range=wind_range)

    line_style = ['-', '--', '-.']
    for count_mld, mld in enumerate(MLD_list):
        # Plotting the distribution according to the Kukulka parametrization
        if kukulka:
            profile_dict = utils_v.get_concentration_list([np.mean(wind_range)], w_rise_list, selection,
                                                          single_select,
                                                          output_step=output_step, diffusion_type='Kukulka',
                                                          boundary=boundary, mld=mld)
            for counter in range(len(profile_dict['concentration_list'])):
                ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'] / correction,
                        label=utils_v.label_MLD_Comparison(parameters=profile_dict['parameter_kukulka'][counter],
                                                           mld=mld, diffusion_type='Kukulka'),
                        linestyle=line_style[counter], color=utils.return_color(count_mld))

        # Plotting the distribution according to the KPP parametrization
        if kpp:
            profile_dict = utils_v.get_concentration_list([np.mean(wind_range)], w_rise_list, selection,
                                                          single_select,
                                                          output_step=output_step, diffusion_type='KPP',
                                                          boundary=boundary, mld=mld)
            for counter in range(len(profile_dict['concentration_list'])):
                ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'] / correction,
                        label=utils_v.label_MLD_Comparison(parameters=profile_dict['parameter_kukulka'][counter],
                                                           mld=mld, diffusion_type='KPP'),
                        linestyle=line_style[counter], color=utils.return_color(count_mld))
    # Adding the legend
    ax.legend(fontsize=legend_size, loc='lower right')

    ax.set_title('w$_{10}' + '={}-{}'.format(*wind_range) + ' m s$^{-1}$', fontsize=ax_label_size)

    plt.savefig(utils_v.mld_comparison_name(diffusion_type, boundary, beaufort=beaufort, close_up=close_up),
                bbox_inches='tight')
