import settings
import numpy as np
import matplotlib.pyplot as plt
import utils
import math


def label_boundary(w_10, diffusion_type, boundary):
    if diffusion_type == 'Kukulka':
        return 'Kukulka et al. (2012),  w_10 = {}'.format(w_10) + ' m s$^{-1}$, ' + boundary
    elif diffusion_type == 'KPP':
        return r'KPP, w_10 = {}'.format(w_10) + 'm s$^{-1}$, MLD = ' + '{} m'.format(settings.MLD) + ', ' + boundary


def diffusion_curve_axis(ax, ax_label_size, w_10, w_rise, profile_dict, diffusion_type, color):
    ax2 = ax.twiny()
    # Get the proper diffusion curve
    depth = profile_dict['depth_bins']
    profile = utils.get_vertical_diffusion_profile(w_10, depth * -1, diffusion_type)
    # X axis = Concentration axis
    ax2.set_xlabel(r'$K_z$ (m$^2$ s$^{-1}$)', fontsize=ax_label_size)
    ax2.set_xlim((0, 0.05))
    # The actual plotting
    ax2.plot(profile, depth, color=color, linestyle='dotted', label=label_diffusivity_profile(w_10, diffusion_type))
    return ax2


def label_time_step(steps, interval):
    t = steps * interval * settings.dt_out.seconds // 3600
    return 't = {} hours'.format(t)


def label_diffusivity_profile(w_10, diffusion_type):
    if diffusion_type == 'Kukulka':
        return 'Kukulka et al. (2012) & Poulain $K_z$,  w_10 = {}'.format(w_10) + ' m s$^{-1}$'
    elif diffusion_type == 'KPP':
        return r'KPP $K_z$, w_10 = {}'.format(w_10) + 'm s$^{-1}$, MLD = ' + '{} m'.format(settings.MLD)


def base_figure(fig_size, ax_range, y_label, x_label, ax_label_size, shape=(1, 1), plot_num=1):
    xmax, xmin, ymax, ymin = ax_range
    fig = plt.figure(figsize=fig_size)
    if shape == (1, 1):
        ax = fig.add_subplot(111)
        # Y axis = Depth axis
        ax.set_ylabel(y_label, fontsize=ax_label_size)
        ax.set_ylim((ymin, ymax))
        ax.tick_params(axis='both', labelsize=ax_label_size)
        # X axis = Concentration axis
        ax.set_xlabel(x_label, fontsize=ax_label_size)
        ax.set_xlim((xmin, xmax))
        return ax
    else:
        ax, fig_index = [], 1
        for row in range(shape[0]):
            for column in range(shape[1]):
                ax_sub = fig.add_subplot(shape[0], shape[1], fig_index)
                # Only add y labels if we are in the first column
                ax_sub.set_ylim((ymin, ymax))
                ax_sub.tick_params(axis='both', labelsize=ax_label_size)
                if column == 0:
                    ax_sub.set_ylabel(y_label, fontsize=ax_label_size)
                else:
                    ax_sub.tick_params(labelleft=False)
                # Only add x labels if we are in the bottom row:
                ax_sub.set_xlim((xmin, xmax))
                if row == (shape[0] - 1) and column % 2 is 1:
                    ax_sub.set_xlabel(x_label, fontsize=ax_label_size)
                else:
                    ax_sub.tick_params(labelbottom=False)
                # Add the axis to the list
                ax.append(ax_sub)
                # Continuing on to the next plot
                fig_index += 1
                if fig_index > plot_num:
                    break
        return tuple(ax)


def saving_filename_basic_profile(save_location, selection, close_up, diffusion_type):
    if close_up is None:
        return save_location + diffusion_type + '_Depth_profile_full_variable={}.png'.format(selection)
    else:
        ymax, ymin = close_up
        return save_location + diffusion_type + '_Depth_profile_max={}_min={}_variable={}.png'.format(ymax, ymin,
                                                                                                      selection)


def saving_filename_time_step(save_location, selection, close_up, diffusion_type):
    if close_up is None:
        return save_location + diffusion_type + '_time_step_full_variable={}.png'.format(selection)
    else:
        ymax, ymin = close_up
        return save_location + diffusion_type + '_time_step_max={}_min={}_variable={}.png'.format(ymax, ymin, selection)


def saving_filename_boundary(save_location, selection, close_up, diffusion_type):
    if close_up is None:
        return save_location + diffusion_type + '_boundary_full_variable={}.png'.format(selection)
    else:
        ymax, ymin = close_up
        return save_location + diffusion_type + '_boundary_max={}_min={}_variable={}.png'.format(ymax, ymin, selection)


def label_profile(selection, parameters):
    w_10, w_rise = parameters
    return r'$w_{10} = $ ' + '{}'.format(w_10) + r' m s$^{-1}$, $w_r = $ ' + '{}'.format(w_rise) + ' m s$^{-1}$'


def label_kukulka(selection, parameters):
    w_10, w_rise = parameters
    if selection is 'w_rise':
        return 'Kukulka et al. (2012) & Poulain, w_rise = {} m s'.format(w_rise) + r'$^{-1}$'
    elif selection is 'w_10':
        return 'Kukulka et al. (2012) & Poulain, w_10 = {} m s'.format(w_10) + r'$^{-1}$'


def label_KPP(parameters):
    w_10, w_rise = parameters
    return r'KPP, w_10 = {}'.format(w_10) + 'm s$^{-1}$, MLD = ' + '{} m'.format(settings.MLD)


def get_axes_range(depth_bins, concentrations):
    output_dict = {'max_depth': np.round(depth_bins.max(), -1), 'min_depth': np.round(depth_bins.min(), -1),
                   'max_count': 0, 'min_count': 0}
    for conc in concentrations:
        if conc.max() > output_dict['max_count']:
            output_dict['max_count'] = conc.max()
    return output_dict


def get_concentration_list(w_10_list, w_rise_list, selection, single_select, diffusion_type, output_step=-1,
                           all_timesteps=False, boundary='Mixed'):
    output_dic = {'concentration_list': [], 'parameter_concentrations': [],
                  'parameter_kukulka': []}
    if selection == 'w_10':
        w_rise_list = [w_rise_list[single_select]]
    elif selection == 'w_rise':
        w_10_list = [w_10_list[single_select]]

    for w_rise in w_rise_list:
        for w_10 in w_10_list:
            # Loading the dictionary containing the concentrations
            input_dir = utils.load_obj(
                utils.get_concentration_output_name(w_10, w_rise, diffusion_type, boundary))
            # Selecting the timeslice of interest
            if output_step == -1:
                concentration = [input_dir['last_time_slice']]
            else:
                concentration = [output_step]
            if all_timesteps:
                concentration = range(input_dir['last_time_slice'] + 1)
            for c in concentration:
                output_dic['concentration_list'].append(input_dir[c] / input_dir[c].max())
            # Get the depth bins of the concentrations
            output_dic['depth_bins'] = input_dir['bin_edges'][:-1] * -1
            # Saving the physical parameters of the concentration distribution we have just loaded
            output_dic['parameter_concentrations'].append((w_10, w_rise))
            # Getting the Kukulka profile and the associated parameters
            output_dic['parameter_kukulka'].append((w_10, w_rise))
    return output_dic


def add_observations(ax, sources=None, wind_range=None, norm_depth=False, alpha=0.5):
    # If no individual sources as specified, just assume I want to include all of them
    if sources is None:
        sources = ['Kooi', 'Pieper', 'Zettler', 'Kukulka']
    data_markers = ['v', 'P', '*', 'o']

    data_labels = {'Kooi': 'Kooi et al. (2016)', 'Pieper': 'Pieper et al. (2020)', 'Zettler': 'PE448 South Atlantic',
                   'Kukulka': 'Kukulka et al. (2012)'}

    # Adding the actual plotting
    for count, source in enumerate(sources):
        # Load the correct data dictionary
        data_dict = utils.load_obj(utils.get_data_output_name(source))
        source_conc = data_dict['concentration']
        # Normalized depths vs un-normalized
        if norm_depth:
            source_depth = data_dict['depth_norm']
        else:
            source_depth = data_dict['depth']
        if wind_range is not None:
            w_min, w_max = wind_range
            wind = data_dict['wind_speed']
            wind_select = (w_min < wind) & (wind <= w_max)
            if sum(wind_select) is 0:
                break
            else:
                source_depth = source_depth[wind_select]
                source_conc = source_conc[wind_select]

        # Making the actual plot
        ax.scatter(source_conc, -1 * source_depth, marker=data_markers[count], c=utils.return_color(count),
                   alpha=alpha, label=data_labels[source])

    lines, labels = ax.get_legend_handles_labels()
    return lines, labels


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
