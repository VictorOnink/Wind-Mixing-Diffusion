import settings
import matplotlib.pyplot as plt
import utils
import visualization.utils_visualization
from visualization import utils_visualization as utils_v


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
                    linestyle='-', color=visualization.utils_visualization.return_color(counter))

    # Plotting the distribution according to the KPP parametrization
    if diffusion_type is 'KPP':
        profile_dict = utils_v.get_concentration_list(w_10_list, w_rise_list, selection, single_select,
                                                      output_step=output_step, diffusion_type=diffusion_type,
                                                      boundary=boundary, alpha_list=alpha_list)
        for counter in range(len(profile_dict['concentration_list'])):
            ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                    label=utils_v.label_KPP(parameters=profile_dict['parameter_kukulka'][counter], selection=selection),
                    linestyle='-', color=visualization.utils_visualization.return_color(counter))
    if diffusion_type is 'artificial':
        profile_dict = utils_v.get_concentration_list(w_10_list, w_rise_list, selection, single_select,
                                                      output_step=output_step, diffusion_type=diffusion_type,
                                                      boundary=boundary, alpha_list=alpha_list)
        for counter in range(len(profile_dict['concentration_list'])):
            ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                    label=utils_v.label_KPP(parameters=profile_dict['parameter_kukulka'][counter], selection=selection),
                    linestyle='-', color=visualization.utils_visualization.return_color(counter))
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
                                                   visualization.utils_visualization.return_color(counter))
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
    plt.savefig(saving_filename_basic_profile(settings.figure_dir, selection, close_up, diffusion_type),
                bbox_inches='tight', dpi=600)


def saving_filename_basic_profile(save_location, selection, close_up, diffusion_type):
    if close_up is None:
        return save_location + diffusion_type + '_Depth_profile_full_variable={}.png'.format(selection)
    else:
        ymax, ymin = close_up
        return save_location + diffusion_type + '_Depth_profile_max={}_min={}_variable={}.png'.format(ymax, ymin,
                                                                                                      selection)