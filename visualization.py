import settings
import matplotlib.pyplot as plt
import utils
import utils_visualization as utils_v


def basic_profile_figure(w_10_list, w_rise_list, selection='w_10', close_up=None,
                         y_label='Depth (m)', x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(8, 8),
                         ax_label_size=16, legend_size=12, single_select=1,
                         output_step=-1, diffusion_type='Kukulka', boundary='Mixed', diffusion_curve=True):
    # Load the relevant data for the figure
    profile_dict = utils_v.get_concentration_list(w_10_list, w_rise_list, selection, single_select,
                                                  output_step=output_step, diffusion_type=diffusion_type,
                                                  boundary=boundary)

    # Preparing for the actual plotting
    range_dict = utils_v.get_axes_range(profile_dict['depth_bins'], profile_dict['concentration_list'])
    xmax, xmin = range_dict['max_count'], range_dict['min_count']
    if close_up == None:
        ymax, ymin = range_dict['max_depth'], range_dict['min_depth']
    else:
        # Allowing for easy close up for a specific part of the depth profile
        ymax, ymin = close_up
    ax_range = (xmax, xmin, ymax, ymin)

    # Creating the axis
    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size)

    # Plotting the distribution according to the Kukulka parametrization
    if diffusion_type is 'Kukulka':
        for counter in range(len(profile_dict['concentration_list'])):
            ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                    label=utils_v.label_kukulka(selection=selection,
                                                parameters=profile_dict['parameter_kukulka'][counter]),
                    linestyle='-', color=utils.return_color(counter))

    # Plotting the distribution according to the KPP parametrization
    if diffusion_type is 'KPP':
        for counter in range(len(profile_dict['concentration_list'])):
            ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                    label=utils_v.label_KPP(parameters=profile_dict['parameter_kukulka'][counter]),
                    linestyle='-', color=utils.return_color(counter))
    lines, labels = ax.get_legend_handles_labels()

    # Plotting the diffusion curve
    if diffusion_curve:
        if diffusion_type is 'Kukulka':
            for counter in range(len(profile_dict['parameter_concentrations'])):
                w_10, w_rise = profile_dict['parameter_concentrations'][counter]
                ax2 = utils_v.diffusion_curve_axis(ax, ax_label_size, w_10, w_rise, profile_dict, diffusion_type,
                                                   utils.return_color(counter))
                lines2, labels2 = ax2.get_legend_handles_labels()
                lines += lines2
                labels += labels2
        elif diffusion_type is 'KPP':
            if selection is 'w_rise':
                range_lim = 1
            if selection is 'w_10':
                range_lim = len(profile_dict['parameter_concentrations'])
            for counter in range(range_lim):
                w_10, w_rise = profile_dict['parameter_concentrations'][counter]
                ax2 = utils_v.diffusion_curve_axis(ax, ax_label_size, w_10, w_rise, profile_dict, diffusion_type,
                                                   'black')
                lines2, labels2 = ax2.get_legend_handles_labels()
                lines += lines2
                labels += labels2

    # Adding the legend
    ax.legend(lines, labels, fontsize=legend_size, loc='lower right')

    # Saving the figure
    plt.savefig(utils_v.saving_filename_basic_profile(settings.figure_dir, selection, close_up, diffusion_type),
                bbox_inches='tight', dpi=600)


def timestep_comparison(w_10_list, w_rise_list, selection='w_10', close_up=None,
                        y_label='Depth (m)', x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(8, 8),
                        ax_label_size=16, legend_size=12, single_select=0,
                        diffusion_type='Kukulka', interval=1, boundary='Mixed', diffusion_curve=True):
    # Load the relevant data for the figure
    profile_dict = utils_v.get_concentration_list(w_10_list, w_rise_list, selection, single_select, diffusion_type,
                                                  all_timesteps=True, boundary=boundary)
    # Preparing for the actual plotting
    range_dict = utils_v.get_axes_range(profile_dict['depth_bins'], profile_dict['concentration_list'])
    xmax, xmin = range_dict['max_count'], range_dict['min_count']
    if close_up is None:
        ymax, ymin = range_dict['max_depth'], range_dict['min_depth']
    else:
        # Allowing for easy close up for a specific part of the depth profile
        ymax, ymin = close_up
    ax_range = (xmax, xmin, ymax, ymin)

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


def boundary_condition_comparison(w_10_list, w_rise_list, selection='w_10', close_up=None,
                                  y_label='Depth (m)', x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(8, 8),
                                  ax_label_size=16, legend_size=12, kukulka=True, model=True, single_select=0,
                                  diffusion_type='KPP', interval=1, boundary='Mixed', diffusion_curve=True):
    # Load the relevant data for the figure
    profile_dict_mix = utils_v.get_concentration_list(w_10_list, w_rise_list, selection, single_select,
                                                      diffusion_type, boundary='Mixed')
    profile_dict_zero = utils_v.get_concentration_list(w_10_list, w_rise_list, selection, single_select,
                                                       diffusion_type, boundary='Zero_Ceiling')
    profile_dict_reflect = utils_v.get_concentration_list(w_10_list, w_rise_list, selection, single_select,
                                                          diffusion_type, boundary='Reflect')
    profile_dict_reduce = utils_v.get_concentration_list(w_10_list, w_rise_list, selection, single_select,
                                                         diffusion_type, boundary='Reduce_dt')
    profile_dict_markov1 = utils_v.get_concentration_list(w_10_list, w_rise_list, selection, single_select,
                                                          diffusion_type, boundary='Reflect_Markov')

    # Preparing for the actual plotting
    range_dict = utils_v.get_axes_range(profile_dict_mix['depth_bins'], profile_dict_mix['concentration_list'])
    xmax, xmin = range_dict['max_count'], range_dict['min_count']
    if close_up == None:
        ymax, ymin = range_dict['max_depth'], range_dict['min_depth']
    else:
        # Allowing for easy close up for a specific part of the depth profile
        ymax, ymin = close_up
    ax_range = (xmax, xmin, ymax, ymin)

    # Creating the axis
    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size)

    # First the mixed boundary layer
    for counter in range(len(profile_dict_mix['concentration_list'])):
        _, w_10, _ = profile_dict_mix['parameter_concentrations'][counter]
        ax.plot(profile_dict_mix['concentration_list'][counter], profile_dict_mix['depth_bins'],
                label=utils_v.label_boundary(w_10, diffusion_type, 'Random Mixed Layer'),
                color=utils.return_color(0))

    # Next the zero ceiling boundary condition
    for counter in range(len(profile_dict_zero['concentration_list'])):
        _, w_10, _ = profile_dict_mix['parameter_concentrations'][counter]
        ax.plot(profile_dict_zero['concentration_list'][counter], profile_dict_zero['depth_bins'],
                label=utils_v.label_boundary(w_10, diffusion_type, 'Zero Ceiling'),
                color=utils.return_color(1))

    # Then, the reflecting boundary condition
    for counter in range(len(profile_dict_reflect['concentration_list'])):
        _, w_10, _ = profile_dict_reflect['parameter_concentrations'][counter]
        ax.plot(profile_dict_reflect['concentration_list'][counter], profile_dict_reflect['depth_bins'],
                label=utils_v.label_boundary(w_10, diffusion_type, 'Reflect'),
                color=utils.return_color(2))
    lines, labels = ax.get_legend_handles_labels()

    # Then, the reducing dt boundary condition
    for counter in range(len(profile_dict_reduce['concentration_list'])):
        _, w_10, _ = profile_dict_reduce['parameter_concentrations'][counter]
        ax.plot(profile_dict_reduce['concentration_list'][counter], profile_dict_reduce['depth_bins'],
                label=utils_v.label_boundary(w_10, diffusion_type, 'Reduce dt'),
                color=utils.return_color(3))

    # Then, the Markov 1 reflecting boundary condition
    for counter in range(len(profile_dict_markov1['concentration_list'])):
        _, w_10, _ = profile_dict_markov1['parameter_concentrations'][counter]
        ax.plot(profile_dict_markov1['concentration_list'][counter], profile_dict_markov1['depth_bins'],
                label=utils_v.label_boundary(w_10, diffusion_type, 'Markov-1, Reflect'),
                color=utils.return_color(4))
    lines, labels = ax.get_legend_handles_labels()

    # Plotting the diffusion curve
    if diffusion_curve:
        w_10, w_rise = profile_dict_mix['parameter_concentrations'][0]
        ax2 = utils_v.diffusion_curve_axis(ax, ax_label_size, w_10, w_rise, profile_dict_mix, diffusion_type, 'black')
        lines2, labels2 = ax2.get_legend_handles_labels()
        lines += lines2
        labels += labels2

    # Adding the legend
    ax.legend(lines, labels, fontsize=legend_size, loc='lower right')
    # Saving the figure
    plt.savefig(utils_v.saving_filename_boundary(settings.figure_dir, selection, close_up, diffusion_type),
                bbox_inches='tight', dpi=600)


def plot_field_data_overview(norm_depth=False, wind_sort=False, y_label='Depth (m)', close_up=None,
                             x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(12, 8), ax_label_size=16,
                             legend_size=12):
    if norm_depth:
        if close_up == None:
            ymax, ymin = 0, -1.5
        else:
            ymax, ymin = close_up
        y_label = 'Depth/MLD'
    else:
        if close_up == None:
            ymax, ymin = 0, -1 * settings.MLD
        else:
            # Allowing for easy close up for a specific part of the depth profile
            ymax, ymin = close_up
    ax_range = (1, 0, ymax, ymin)

    # Plotting all data points, with no sorting based on wind conditions
    if not wind_sort:
        # Get the base figure axis
        ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size)

        # Plotting the data points
        legend_line, legend_label = utils_v.add_observations(ax, norm_depth=norm_depth)

        # Adding the legend
        ax.legend(legend_line, legend_label, fontsize=legend_size, loc='lower right')

    # Plotting data points, split over multiple plots according to Beaufort wind scale
    if wind_sort:
        # Titles for the subplots
        sub_titles = ['(a) Beaufort 1', '(b) Beaufort 2', '(c) Beaufort 3', '(d) Beaufort 4', '(e) Beaufort 5',
                      '(f) Beaufort 6']
        # Get the base figure axis
        plot_num = 6
        ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size, shape=(2, 3), plot_num=plot_num)
        beaufort = utils.beaufort_limits()
        for scale in range(plot_num):
            line, label = utils_v.add_observations(ax[scale], norm_depth=norm_depth, wind_range=beaufort[scale + 1])
            ax[scale].set_title(sub_titles[scale], fontsize=ax_label_size)

        # Adding the legend
        ax[-1].legend(line, label, fontsize=legend_size, loc='lower right')

    plt.savefig(utils_v.field_data_figure_names(close_up, wind_sort, norm_depth), bbox_inches='tight')


def plot_model_field_data_comparison(w_10_list, w_rise_list, selection='w_10', output_step=-1, single_select=1,
                                     norm_depth=False, wind_sort=True, y_label='Depth (m)', close_up=None,
                                     x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(16, 8), ax_label_size=16,
                                     legend_size=12, diffusion_type='Kukulka', boundary='Reflect_Markov',alpha=0.3):
    if norm_depth:
        if close_up == None:
            ymax, ymin = 0, -1.5
        else:
            ymax, ymin = close_up
        y_label = 'Depth/MLD'
    else:
        if close_up == None:
            ymax, ymin = 0, -1 * settings.max_depth
        else:
            # Allowing for easy close up for a specific part of the depth profile
            ymax, ymin = close_up
    ax_range = (1, 0, ymax, ymin)

    # Plotting all data points, with no sorting based on wind conditions
    if not wind_sort:
        profile_dict = utils_v.get_concentration_list(w_10_list, w_rise_list, selection, single_select,
                                                      output_step=output_step, diffusion_type=diffusion_type,
                                                      boundary=boundary)
        # Get the base figure axis
        ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size)

        # Plotting the field data points
        legend_line, legend_label = utils_v.add_observations(ax, norm_depth=norm_depth, alpha=alpha)

        # Plotting the distribution according to the Kukulka parametrization
        if diffusion_type is 'Kukulka':
            for counter in range(len(profile_dict['concentration_list'])):
                ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                        label=utils_v.label_kukulka(selection=selection,
                                                    parameters=profile_dict['parameter_kukulka'][counter]),
                        linestyle='-', color=utils.return_color(counter))

        # Plotting the distribution according to the KPP parametrization
        if diffusion_type is 'KPP':
            for counter in range(len(profile_dict['concentration_list'])):
                ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                        label=utils_v.label_KPP(parameters=profile_dict['parameter_kukulka'][counter]),
                        linestyle='-', color=utils.return_color(counter))
        lines, labels = ax.get_legend_handles_labels()

        # Adding the legend
        ax.legend(fontsize=legend_size, loc='lower right')

    # Plotting data points, split over multiple plots according to Beaufort wind scale
    if wind_sort:
        profile_dict = utils_v.get_concentration_list(w_10_list, w_rise_list, 'all', single_select,
                                                      output_step=output_step, diffusion_type=diffusion_type,
                                                      boundary=boundary)
        # Titles for the subplots
        sub_titles = ['(a) Beaufort 1', '(b) Beaufort 2', '(c) Beaufort 3', '(d) Beaufort 4', '(e) Beaufort 5',
                      '(f) Beaufort 6']
        # Get the base figure axis
        plot_num = 6
        ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size, shape=(2, 3), plot_num=plot_num)
        beaufort = utils.beaufort_limits()

        for scale in range(plot_num - 1):
            line, label = utils_v.add_observations(ax[scale], norm_depth=norm_depth, wind_range=beaufort[scale + 1],
                                                   alpha=alpha)
            ax[scale].set_title(sub_titles[scale], fontsize=ax_label_size)
            for model in range(len(w_rise_list)):
                _, w_rise = profile_dict['parameter_kukulka'][scale + model * len(w_10_list)]
                concentration = profile_dict['concentration_list'][scale + model * len(w_10_list)]
                ax[scale].plot(concentration, profile_dict['depth_bins'],
                               label=utils_v.label_model_field_comparison(w_rise, diffusion_type, boundary),
                               linestyle='-', color=utils.return_color(model))
            lines, labels = ax[scale].get_legend_handles_labels()

        # Adding the legend to the last plot, and hiding the grid
        ax[-1].legend(lines, labels, fontsize=legend_size, loc='lower right')
        ax[-1].axis('off')

    # Saving the figure
    save_name = utils_v.model_field_data_comparison_name(diffusion_type, boundary, close_up=close_up,
                                                         wind_sort=wind_sort, norm_depth=norm_depth)
    plt.savefig(save_name, bbox_inches='tight')

