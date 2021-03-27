import settings
import matplotlib.pyplot as plt
import utils
import utils_visualization as utils_v
import numpy as np
from netCDF4 import Dataset
import analysis


def basic_profile_figure(w_10_list, w_rise_list, alpha_list, selection='w_10', close_up=None,
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
                                                      alpha_list=alpha_list, output_step=output_step,
                                                      diffusion_type=diffusion_type, boundary=boundary)
        for counter in range(len(profile_dict['concentration_list'])):
            ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                    label=utils_v.label_kukulka(selection=selection,
                                                parameters=profile_dict['parameter_kukulka'][counter]),
                    linestyle='-', color=utils.return_color(counter))

    # Plotting the distribution according to the KPP parametrization
    if diffusion_type is 'KPP':
        profile_dict = utils_v.get_concentration_list(w_10_list, w_rise_list, selection, single_select,
                                                      output_step=output_step, diffusion_type=diffusion_type,
                                                      boundary=boundary, alpha_list=alpha_list)
        for counter in range(len(profile_dict['concentration_list'])):
            ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                    label=utils_v.label_KPP(parameters=profile_dict['parameter_kukulka'][counter], selection=selection),
                    linestyle='-', color=utils.return_color(counter))
    if diffusion_type is 'artificial':
        profile_dict = utils_v.get_concentration_list(w_10_list, w_rise_list, selection, single_select,
                                                      output_step=output_step, diffusion_type=diffusion_type,
                                                      boundary=boundary, alpha_list=alpha_list)
        for counter in range(len(profile_dict['concentration_list'])):
            ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                    label=utils_v.label_KPP(parameters=profile_dict['parameter_kukulka'][counter], selection=selection),
                    linestyle='-', color=utils.return_color(counter))
    lines, labels = ax.get_legend_handles_labels()

    # Plotting the diffusion curve
    kukulka, kpp, artificial = utils_v.boolean_diff_type(diffusion_type)
    if diffusion_curve:
        if kukulka:
            profile_dict = utils_v.get_concentration_list(w_10_list, w_rise_list, selection, single_select,
                                                          output_step=output_step, diffusion_type='Kukulka',
                                                          boundary=boundary)
            for counter in range(len(profile_dict['parameter_concentrations'])):
                w_10, w_rise = profile_dict['parameter_concentrations'][counter]
                ax2 = utils_v.diffusion_curve_axis(ax, ax_label_size, w_10, profile_dict, 'Kukulka',
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
            elif selection is 'w_10':
                range_lim = len(profile_dict['parameter_concentrations'])

            for counter in range(range_lim):
                w_10, w_rise = profile_dict['parameter_concentrations'][counter]
                ax2 = utils_v.diffusion_curve_axis(ax, ax_label_size, w_10, profile_dict, "KPP",
                                                   'black', gradient=True)
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
    ax_range = (1.55, 0, ymax, ymin)
    depth = np.linspace(ymax, np.abs(ymin), 1000)

    # Creating the axis
    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size)

    # Plotting the diffusion profile according to the Kukulka approach
    for count, w_10 in enumerate(w_10_list):
        profile = utils.get_vertical_diffusion_profile(w_10, depth, 'Kukulka')
        ax.plot(profile * 100, -1 * depth, color=utils_v.discrete_color_from_cmap(count, len(w_10_list)), linestyle='-',
                label=utils_v.label_diffusivity_profile(w_10, 'Kukulka'))

    # Plotting the diffusion profile according the KPP approach
    for count, w_10 in enumerate(w_10_list):
        profile = utils.get_vertical_diffusion_profile(w_10, depth, 'KPP')
        ax.plot(profile * 100, -1 * depth, color=utils_v.discrete_color_from_cmap(count, len(w_10_list)),
                linestyle='--', label=utils_v.label_diffusivity_profile(w_10, 'KPP'))

    # Adding the legend
    ax.legend(fontsize=legend_size, loc='lower right')

    plt.savefig(settings.figure_dir + 'diffusion_profile_MLD=20.png', bbox_inches='tight')


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
    plt.savefig(utils_v.saving_filename_time_step(settings.figure_dir, close_up, diffusion_type),
                bbox_inches='tight', dpi=600)


def boundary_condition_comparison(w_rise_list, alpha_list, selection='w_10', close_up=None, output_step=-1,
                                  y_label='Depth (m)', x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(8, 8),
                                  ax_label_size=16, legend_size=12, single_select=0,
                                  diffusion_type='KPP', interval=1, alpha=0.3):
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=False)

    # Selecting which model data we want to plot based on the diffusion type
    kukulka, kpp, artificial = utils_v.boolean_diff_type(diffusion_type)
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
                                                          boundary=boundary, alpha_list=alpha_list)
            for counter in range(len(profile_dict['concentration_list'])):
                ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                        label=utils_v.label_boundary(w_10=mean_wind, diffusion_type='Kukulka', boundary=boundary),
                        linestyle=line_style[count], color=utils.return_color(count % 3))

        # Plotting the distribution according to the KPP parametrization
        if kpp:
            profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                          output_step=output_step, diffusion_type='KPP',
                                                          boundary=boundary, alpha_list=alpha_list)
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
                                                 diffusion_type), bbox_inches='tight', dpi=600)


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
        beaufort = utils.beaufort_limits()
        for scale in range(plot_num):
            line, label = utils_v.add_observations(ax[scale], norm_depth=norm_depth, wind_range=beaufort[scale + 1],
                                                   alpha=alpha)
            ax[scale].set_title(sub_titles[scale], fontsize=ax_label_size)

        # Adding the legend
        ax[-1].legend(line, label, fontsize=legend_size, loc='lower right')
        ax[-1].axis('off')

    plt.savefig(utils_v.field_data_figure_names(close_up, wind_sort, norm_depth), bbox_inches='tight')


def plot_model_field_data_comparison(w_10_list, w_rise_list, alpha_list, selection='w_10', output_step=-1,
                                     single_select=1,
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
    kukulka, kpp, artificial = utils_v.boolean_diff_type(diffusion_type)
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
                                                              boundary=boundary, alpha_list=alpha_list)
                for counter in range(len(profile_dict['concentration_list'])):
                    ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'] / correction,
                            label=utils_v.label_kukulka(selection=selection,
                                                        parameters=profile_dict['parameter_kukulka'][counter]),
                            linestyle='-', color=utils.return_color(counter))

            # Plotting the distribution according to the KPP parametrization
            if kpp:
                profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                              output_step=output_step, diffusion_type='KPP',
                                                              boundary=boundary, alpha_list=alpha_list)
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
        sub_titles = [r'(a) u$_{10}$=0.2-1.5 m s$^{-1}$', r'(b) u$_{10}$=1.5-3.3 m s$^{-1}$',
                      r'(c) u$_{10}$=3.3-5.4 m s$^{-1}$', r'(d) u$_{10}$=5.4-7.9 m s$^{-1}$',
                      r'(e) u$_{10}$=7.9-10.7 m s$^{-1}$',
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
                                                                  boundary=boundary, alpha_list=alpha_list)
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
                                                                  boundary=boundary, alpha_list=alpha_list)
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
    save_name = utils_v.model_field_data_comparison_name(diffusion_type, boundary, alpha_list, close_up=close_up,
                                                         wind_sort=wind_sort, norm_depth=norm_depth, beaufort=beaufort)
    plt.savefig(save_name, bbox_inches='tight')


def mld_depth_influence(w_rise_list, MLD_list, alpha_list, selection='w_rise', output_step=-1, single_select=0,
                        y_label='Depth/MLD', close_up=None, beaufort=5,
                        x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(16, 8), ax_label_size=16,
                        legend_size=12, diffusion_type='KPP', boundary='Reflect', alpha=0.3):
    correction = settings.MLD
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=True)

    # Selecting which model data we want to plot based on the diffusion type
    kukulka, kpp, artificial = utils_v.boolean_diff_type(diffusion_type)

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
                                                          single_select, alpha_list=alpha_list,
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
                                                          single_select, alpha_list=alpha_list,
                                                          output_step=output_step, diffusion_type='KPP',
                                                          boundary=boundary, mld=mld)
            for counter in range(len(profile_dict['concentration_list'])):
                ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'] / correction,
                        label=utils_v.label_MLD_Comparison(parameters=profile_dict['parameter_kukulka'][counter],
                                                           mld=mld, diffusion_type='KPP'),
                        linestyle=line_style[counter], color=utils.return_color(count_mld))
    # Adding the legend
    ax.legend(fontsize=legend_size, loc='lower right')

    ax.set_title(r'u$_{10}$' + '={}-{}'.format(*wind_range) + r' m s$^{-1}$', fontsize=ax_label_size)

    plt.savefig(utils_v.mld_comparison_name(diffusion_type, boundary, beaufort=beaufort, close_up=close_up),
                bbox_inches='tight')


def Markov_alpha_dependence(w_rise_list, selection='w_10', close_up=None, y_label='Depth (m)', alpha=0.3,
                            x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(8, 8),
                            ax_label_size=16, legend_size=12, single_select=1,
                            output_step=-1, diffusion_type='Kukulka'):
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=False)

    # Selecting which model data we want to plot based on the diffusion type
    kukulka, kpp, artificial = utils_v.boolean_diff_type(diffusion_type)
    # Selecting which model data we want to plot based on the diffusion scheme
    boundary_list = ['Reflect', 'Reflect_Markov', 'Reflect_Markov', 'Reflect_Markov', 'Reflect_Markov',
                     'Reflect_Markov', 'Reflect_Markov']
    alpha_list = [0, 0.0, 0.1, 0.3, 0.5, 0.7, 0.95]

    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size)
    line_style = ['-', '--', '--', '--', '--', '--', '--']

    mean_wind = np.mean(utils.beaufort_limits()[4])
    _, _ = utils_v.add_observations(ax, norm_depth=False, alpha=alpha, wind_range=utils.beaufort_limits()[4])

    for count, boundary in enumerate(boundary_list):
        # Plotting the distribution according to the KPP parametrization
        if kukulka:
            profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                          output_step=output_step, diffusion_type='Kukulka',
                                                          boundary=boundary, alpha_list=alpha_list[count])
            for counter in range(len(profile_dict['concentration_list'])):
                ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                        label=utils_v.label_alpha_comparison(boundary=boundary, alpha=alpha_list[count]),
                        linestyle=line_style[count], color=utils.return_color(count))
        if kpp:
            profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                          output_step=output_step, diffusion_type='KPP',
                                                          boundary=boundary, alpha_list=alpha_list[count])
            for counter in range(len(profile_dict['concentration_list'])):
                ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                        label=utils_v.label_alpha_comparison(boundary=boundary, alpha=alpha_list[count]),
                        linestyle=line_style[count], color=utils.return_color(count))
        if artificial:
            profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                          output_step=output_step, diffusion_type='artificial',
                                                          boundary=boundary, alpha_list=alpha_list[count])
            for counter in range(len(profile_dict['concentration_list'])):
                ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                        label=utils_v.label_alpha_comparison(boundary=boundary, alpha=alpha_list[count]),
                        linestyle=line_style[count], color=utils.return_color(count))

    lines, labels = ax.get_legend_handles_labels()
    # Adding the legend
    ax.legend(lines[:-5], labels[:-5], fontsize=legend_size, loc='lower right')
    ax.set_title(r'u$_{10}$ = 5.4-7.9 m s$^{-1}$ - $\alpha$', fontsize=ax_label_size)

    # Saving the figure
    str_format = diffusion_type, w_rise_list[0], settings.MLD
    plt.savefig(settings.figure_dir + '{}_alpha_check_w_rise={}_mld={}'.format(*str_format) + '.png',
                bbox_inches='tight', dpi=600)


def diffusion_markov_comparison(w_rise_list, selection='w_10', close_up=None, y_label='Depth (m)', alpha=0.3,
                                x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(10, 6),
                                ax_label_size=16, legend_size=12, single_select=1,
                                output_step=-1):
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=False)

    # Selecting which model data we want to plot based on the diffusion type
    kukulka, kpp, artificial = utils_v.boolean_diff_type('all')
    # Selecting which model data we want to plot based on the diffusion scheme
    boundary_list = ['Reflect', 'Reflect_Markov', 'Reflect_Markov', 'Reflect_Markov', 'Reflect_Markov',
                     'Reflect_Markov', 'Reflect_Markov']
    alpha_list = [0, 0.0, 0.1, 0.3, 0.5, 0.7, 0.95]
    line_style = ['-', '--', '--', '--', '--', '--', '--']

    # Creating the axis, one for KPP and one for Kukulka
    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size, shape=(1, 2), plot_num=2,
                             all_x_labels=True)

    mean_wind = np.mean(utils.beaufort_limits()[4])
    for axis in ax:
        data_line, data_label = utils_v.add_observations(axis, norm_depth=False, alpha=alpha,
                                                         wind_range=utils.beaufort_limits()[4])

    # First, plotting the KPP data:
    for count, boundary in enumerate(boundary_list):
        if kpp:
            profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                          output_step=output_step, diffusion_type='KPP',
                                                          boundary=boundary, alpha_list=alpha_list[count])
            for counter in range(len(profile_dict['concentration_list'])):
                ax[0].plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                           label=utils_v.label_alpha_comparison(boundary=boundary, alpha=alpha_list[count]),
                           linestyle=line_style[count], color=utils.return_color(count))
        if kukulka:
            profile_dict = utils_v.get_concentration_list([mean_wind], w_rise_list, selection, single_select,
                                                          output_step=output_step, diffusion_type='Kukulka',
                                                          boundary=boundary, alpha_list=alpha_list[count])
            for counter in range(len(profile_dict['concentration_list'])):
                ax[1].plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                          label=utils_v.label_alpha_comparison(boundary=boundary, alpha=alpha_list[count]),
                          linestyle=line_style[count], color=utils.return_color(count))

    lines, labels = ax[1].get_legend_handles_labels()
    # Adding the legends
    ax[0].legend(data_line, data_label, fontsize=legend_size, loc='lower right')
    ax[1].legend(lines[:-5], labels[:-5], fontsize=legend_size, loc='lower right')

    # Adding subplot titles
    ax[0].set_title(r'KPP', fontsize=ax_label_size)
    ax[1].set_title(r'PZK', fontsize=ax_label_size)

    # Saving the figure
    str_format = w_rise_list[0], settings.MLD
    plt.savefig(settings.figure_dir + 'markov_diffusion_check_w_rise={}_mld={}'.format(*str_format) + '.png',
                bbox_inches='tight', dpi=600)


def sanity_check(w_10, w_rise, diffusion_type, boundary, alpha_list):
    dataset = Dataset(utils.get_parcels_output_name(w_10, w_rise, diffusion_type, boundary=boundary, mld=settings.MLD,
                                                    alpha=alpha_list[0]))
    start, end = 0, -1
    depth = dataset.variables['z'][0, start:end]
    ratio = dataset.variables['ratio'][0, start:end]
    sig2 = dataset.variables['sig2_store'][0, start:end]
    w_T = dataset.variables['w_T_store'][0, start:end]
    w_T2 = np.square(w_T)
    dWz = dataset.variables['dWz_store'][0, start:end]
    dHis = dataset.variables['d_History_store'][0, start:end]
    dGrad = dataset.variables['d_Gradient_store'][0, start:end]
    dDiff = dataset.variables['d_diff_store'][0, start:end]

    dim = 8
    fig = plt.figure(figsize=(dim, dim))

    var = [depth, sig2, w_T, w_T2, ratio, dHis, dGrad, dDiff]
    title = ['depth', 'sig2', 'w_T', 'w_T2', 'ratio', 'dHis', 'dGrad', 'dDiff']
    shape = (len(var), 1)
    count = 1

    for row in range(shape[0]):
        for column in range(shape[1]):
            ax_sub = fig.add_subplot(shape[0], shape[1], count)
            ax_sub.plot(var[count - 1])
            ax_sub.set_ylabel(title[count - 1])
            count += 1
    # ax_sub.set_ylim([0, 600])

    plt.tight_layout()
    plt.savefig(settings.figure_dir + '/sanity')


def markov_0_RMSE_comparison(y_label=r'$u_{10}$ (m s$^{-1}$)', x_label=r'RMSE', fig_size=(8, 8),
                             ax_label_size=16, legend_size=12):
    w_10 = [0.85, 2.4, 4.35, 6.65, 9.3]
    w_r = [-0.03, -0.003, -0.0003]
    boundary = 'Reflect'
    diffusion = ['KPP', 'Kukulka']
    diffusion_offset = {'KPP': 0.1, 'Kukulka': -0.1}
    # Selecting the marker type according which form of diffusion it is:
    marker_type = {'KPP': 'o', 'Kukulka': 'X', 'PZK': 'X'}
    # Selecting the marker color according to which rise velocity
    marker_color = {-0.03: 'tab:blue', -0.003: 'tab:red', -0.0003: 'tab:green'}

    # Looping through the simulations, and retrieving the RMSE values for them
    point_list = []
    for index_w10, wind in enumerate(w_10):
        for rise in w_r:
            for diffusion_type in diffusion:
                RMSE = analysis.determine_RMSE(wind, rise, diffusion_type, boundary, alpha=0.0, output=True)
                plot_tuple = RMSE, index_w10 + 1 + diffusion_offset[diffusion_type], marker_type[diffusion_type], marker_color[rise]
                point_list.append(plot_tuple)

    # Creating the axis
    fig = plt.figure(figsize=fig_size)
    ax = fig.add_subplot(111)
    # Y axis = Depth axis
    ax.set_ylabel(y_label, fontsize=ax_label_size)
    ax.set_ylim((0, 6))
    ax.tick_params(axis='both', labelsize=ax_label_size)
    # X axis = Concentration axis
    ax.set_xlabel(x_label, fontsize=ax_label_size)
    ax.set_xlim((0, 1))

    # Now, plotting the points
    for point in point_list:
        RMSE, index_w10, marker, color = point
        ax.plot(RMSE, index_w10, color=color, marker=marker, alpha=0.7, markersize=10)

    # Now, altering the Y axis to list the wind speeds instead of the simple labels 1 - 5
    ax.set_yticks(range(7))
    ax.set_yticklabels(['', 0.85, 2.40, 4.35, 6.65, 9.30, ''])

    # Next, adding a legend to explain the color scheme and the marker type
    label_marker = ['KPP', 'PZK']
    marker = [plt.plot([], [], c='k', markersize=10, marker=marker_type[label], label=label, linestyle='')[0] for label in
              label_marker]

    def label_color(rise):
        return r'$w_r$ = {}'.format(rise) + ' m s$^{-1}$'

    color = [plt.plot([], [], c=marker_color[rise], markersize=10, marker='o', label=label_color(rise), linestyle='')[0] for rise in
              w_r]

    ax.legend(handles=marker + color, fontsize=legend_size, loc='lower right')

    plt.savefig(settings.figure_dir + 'model_evaluation_markov_0.png', bbox_inches='tight', dpi=600)

