import settings
import numpy as np
import matplotlib.pyplot as plt
import utils
from matplotlib.gridspec import GridSpec


def label_boundary(w_10, diffusion_type, boundary):
    boundary_dict = {'Mixed': 'M-0 Mixed', 'Reflect': 'M-0 Reflect', 'Reduce_dt': 'M-0 Reduce dt',
                     'Mixed_Markov': 'M-1 Mixed', 'Reflect_Markov':'M-1 Reflect',
                     'Reduce_dt_Markov':'M-1 Reduce dt'}
    if diffusion_type == 'Kukulka':
        return r'SWB, u$_{10}$ '+'= {:.2f}'.format(w_10) + ' m s$^{-1}$, ' + boundary_dict[boundary]
    elif diffusion_type == 'KPP':
        return r'KPP, u$_{10}$ '+'= {:.2f}'.format(w_10) + 'm s$^{-1}$, MLD = ' + '{} m'.format(settings.MLD) + ', ' \
               + boundary_dict[boundary]


def diffusion_curve_axis(ax, ax_label_size, w_10, profile_dict, diffusion_type, color, linestyle='dotted',
                         gradient=False):
    ax2 = ax.twiny()
    # Get the proper diffusion curve
    depth = profile_dict['depth_bins']
    profile = utils.get_vertical_diffusion_profile(w_10, depth * -1, diffusion_type)
    # X axis = Concentration axis
    ax2.set_xlabel(r'$K_z$ (m$^2$ s$^{-1}$)', fontsize=ax_label_size)
    ax2.set_xlim((0, 0.05))
    # The actual plotting
    ax2.plot(profile, depth, color=color, linestyle=linestyle, label=label_diffusivity_profile(w_10, diffusion_type))
    if gradient:
        ax2.set_xlim((-0.005, 0.05))
        profile = utils.get_vertical_diffusion_gradient_profile(w_10, depth * -1, diffusion_type)
        ax2.plot(profile, depth, color='red', linestyle='-.',
                 label=label_diffusivity_profile(w_10, diffusion_type))
    return ax2


def label_diffusivity_profile(w_10, diffusion_type):
    if diffusion_type == 'Kukulka':
        return r'SWB, u$_{10}$' + ' = {:.2f}'.format(w_10) + ' m s$^{-1}$'
    elif diffusion_type == 'KPP':
        return r'KPP, u$_{10}$ '+'= {:.2f}'.format(w_10) + 'm s$^{-1}$, MLD = ' + '{} m'.format(settings.MLD)


def base_figure(fig_size, ax_range, y_label, x_label, ax_label_size, shape=(1, 1), plot_num=1, all_x_labels=False,
                legend_axis=False):
    xmax, xmin, ymax, ymin = ax_range
    fig = plt.figure(figsize=fig_size)
    if shape == (1, 1):
        if legend_axis:
            grid = GridSpec(nrows=shape[0], ncols=shape[1] + 1, figure=fig)
        else:
            grid = GridSpec(nrows=shape[0], ncols=shape[1], figure=fig)
        ax_sub = fig.add_subplot(grid[0, 0])
        # Y axis = Depth axis
        ax_sub.set_ylabel(y_label, fontsize=ax_label_size)
        ax_sub.set_ylim((ymin, ymax))
        ax_sub.tick_params(axis='both', labelsize=ax_label_size)
        # X axis = Concentration axis
        ax_sub.set_xlabel(x_label, fontsize=ax_label_size)
        ax_sub.set_xlim((xmin, xmax))
        ax_sub.set_xscale('log')
        if not legend_axis:
            return ax_sub
        else:
            ax = [ax_sub]
            ax_legend = fig.add_subplot(grid[0, -1])
            ax_legend.set_axis_off()
            ax.append(ax_legend)
            return ax
    else:
        if legend_axis:
            grid = GridSpec(nrows=shape[0], ncols=shape[1] + 1, figure=fig)
        else:
            grid = GridSpec(nrows=shape[0], ncols=shape[1], figure=fig)
        ax, fig_index = [], 1
        for row in range(shape[0]):
            for column in range(shape[1]):
                # ax_sub = fig.add_subplot(shape[0], shape[1], fig_index)
                ax_sub = fig.add_subplot(grid[row, column])
                # Only add y labels if we are in the first column
                ax_sub.set_ylim((ymin, ymax))
                ax_sub.tick_params(axis='both', labelsize=ax_label_size)
                if column == 0:
                    ax_sub.set_ylabel(y_label, fontsize=ax_label_size)
                else:
                    ax_sub.tick_params(labelleft=False)
                # Only add x labels if we are in the bottom row:
                ax_sub.set_xlim((xmin, xmax))
                ax_sub.set_xscale('log')
                if row == (shape[0] - 1):
                    if not all_x_labels and column % 2 is 1:
                        ax_sub.set_xlabel(x_label, fontsize=ax_label_size)
                    elif all_x_labels:
                        ax_sub.set_xlabel(x_label, fontsize=ax_label_size)
                    else:
                        ax_sub.tick_params(labelbottom=False)
                else:
                    ax_sub.tick_params(labelbottom=False)
                # Add the axis to the list
                ax.append(ax_sub)
                # Continuing on to the next plot
                fig_index += 1
                if fig_index > plot_num:
                    break
        if legend_axis:
            for rows in range(shape[0]):
                ax_legend = fig.add_subplot(grid[rows, -1])
                ax_legend.set_axis_off()
                ax.append(ax_legend)
        return tuple(ax)


def label_kukulka(selection, parameters):
    w_10, w_rise = parameters
    w_rise = np.abs(w_rise)
    if selection is 'w_rise':
        return r'SWB, w$_{rise}$ '+'= {} m s'.format(w_rise) + r'$^{-1}$'
    elif selection is 'w_10':
        return r'SWB, u$_{10}$ '+'= {} m s'.format(w_10) + r'$^{-1}$'


def label_KPP(selection, parameters, mld=settings.MLD):
    w_10, w_rise = parameters
    w_rise = np.abs(w_rise)
    if selection is 'w_rise':
        return r'KPP, w$_{rise}$ '+'= {}'.format(w_rise) + 'm s$^{-1}$, MLD = ' + '{:.1f} m'.format(mld)
    elif selection is 'w_10':
        return r'KPP, u$_{10}$ '+'= {}'.format(w_10) + 'm s$^{-1}$, MLD = ' + '{:.1f} m'.format(mld)


def label_alpha_comparison(boundary, alpha):
    if boundary is 'Reflect':
        return 'M-0'
    else:
        return r'M-1, $\alpha$ = ' + '{}'.format(alpha)


def get_axes_range(close_up, norm_depth, delta_x=0.01, delta_y=0.5):
    # The x-axis shows the relative concentration
    xmax, xmin = 1 + delta_x, 0 - delta_x
    if norm_depth:
        delta_y = 0.1
    # If we have a close_up plot
    if close_up is not None:
        ymax, ymin = close_up
        ymax += delta_y
    else:
        # If we have normalized depths
        if norm_depth:
            ymax, ymin = 0 + delta_y, -1
        else:
            ymax, ymin = 0 + delta_y, -1 * settings.max_depth - delta_y
    return xmax, xmin, ymax, ymin


def determine_linestyle(boundary, boundary_list, kpp, kukulka, diffusion_type):
    if len(boundary_list) is 1:
        if kukulka and kpp:
            if diffusion_type is 'Kukulka':
                return '-'
            elif diffusion_type is 'KPP':
                return '--'
        elif kukulka:
            return '-'
        elif kpp:
            return '--'
    else:
        line_style = {'Reflect': '-', 'Reflect_Markov': '--'}
        return line_style[boundary]


def boolean_diff_type(diffusion_type):
    if diffusion_type is 'Kukulka':
        kukulka, kpp, artificial = True, False, False
    elif diffusion_type is 'KPP':
        kukulka, kpp, artificial = False, True, False
    elif diffusion_type is 'all':
        kukulka, kpp, artificial = True, True, False
    elif diffusion_type is 'artificial':
        kukulka, kpp, artificial = False, False, True
    return kukulka, kpp, artificial


def get_concentration_list(w_10_list, w_rise_list, selection, single_select, diffusion_type, alpha_list, output_step=-1,
                           all_timesteps=False, boundary='Mixed', mld=settings.MLD, dt=settings.dt_int.seconds):
    output_dic = {'concentration_list': [], 'parameter_concentrations': [],
                  'parameter_kukulka': []}
    if selection == 'w_10':
        w_rise_list = [w_rise_list[single_select]]
    elif selection == 'w_rise':
        w_10_list = [w_10_list[single_select]]
    # selection == 'all' will just return all simulations for a particular diffusion_type and boundary
    if type(alpha_list) is list:
        alpha = alpha_list[0]
    else:
        alpha = alpha_list

    for w_rise in w_rise_list:
        for w_10 in w_10_list:
            # Loading the dictionary containing the concentrations
            input_dir = utils.load_obj(utils.get_concentration_output_name(w_10, w_rise, diffusion_type, boundary,
                                                                           mld=mld, alpha=alpha, dt=dt))
            # Selecting the timeslice of interest
            if output_step == -1:
                concentration = [input_dir['last_time_slice']]
            else:
                concentration = [output_step]
            if all_timesteps:
                concentration = range(input_dir['last_time_slice'] + 1)
            for c in concentration:
                output_dic['concentration_list'].append(input_dir[c] / input_dir[c].sum())
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
        sources = ['Kooi', 'Pieper', 'Zettler', 'Kukulka', 'Egger']
    data_markers = ['v', 'P', '*', 'o', 's']

    data_labels = {'Kooi': 'Kooi et al. (2016)', 'Pieper': 'Pieper et al. (2020)',
                   'Zettler': 'Amaral-Zettler (Unpublished data)',
                   'Kukulka': 'Kukulka et al. (2012)', 'Egger': 'Egger et al. (2020)'}

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
        ax.scatter(source_conc, -1 * source_depth, marker=data_markers[count], c='black',
                   alpha=alpha, label=data_labels[source])

    lines, labels = ax.get_legend_handles_labels()
    return lines, labels


def discrete_color_from_cmap(index, subdivisions, cmap='viridis_r'):
    cmap_steps = plt.cm.get_cmap(cmap, subdivisions)
    return cmap_steps(index)


def rise_velocity_selector(size_class, w_rise):
    if size_class is 'large':
        if -0.0003 in w_rise:
            w_rise.remove(-0.0003)
    return w_rise


def return_color(index):
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray',
              'tab:olive', 'tab:cyan']
    num = len(colors)
    return colors[index % num]