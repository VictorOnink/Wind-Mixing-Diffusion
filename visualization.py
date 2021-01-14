import settings
import numpy as np
import matplotlib.pyplot as plt
import utils
import math
from datetime import timedelta


def quick_plot(x, y):
    plt.plot(x[:-1], y, '.', c='r')
    plt.xlim((x.min(), x.max()))
    plt.ylim(y.min(), y.max())
    plt.show()


def basic_profile_figure(k_z_list, w_10_list, w_rise_list, selection='k_z', close_up=None,
                         y_label='Depth (m)', x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(8, 8),
                         ax_label_size=16, legend_size=12, kukulka=True, model=True, single_select=1,
                         output_step=-1, diffusion_type='Rouse', boundary='Mixed', diffusion_curve=True):
    # Load the relevant data for the figure
    profile_dict = get_concentration_list(k_z_list, w_10_list, w_rise_list, selection, single_select,
                                          output_step=output_step, diffusion_type=diffusion_type, boundary=boundary)
    # Preparing for the actual plotting
    range_dict = get_axes_range(profile_dict['depth_bins'], profile_dict['concentration_list'])
    xmax, xmin = range_dict['max_count'], range_dict['min_count']
    if close_up == None:
        ymax, ymin = range_dict['max_depth'], range_dict['min_depth']
    else:
        # Allowing for easy close up for a specific part of the depth profile
        ymax, ymin = close_up
    ax_range = (xmax, xmin, ymax, ymin)
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray',
              'tab:olive', 'tab:cyan']

    # Creating the axis
    ax = base_figure(fig_size, ax_range, y_label, x_label, ax_label_size)

    # Plotting the distribution from the Rouse approach
    if model:
        for counter in range(len(profile_dict['concentration_list'])):
            ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                    label=label_profile(selection, parameters=profile_dict['parameter_concentrations'][counter]),
                    color=colors[counter])
    # Plotting the distribution according to the Kukulka parametrization
    if kukulka:
        for counter in range(len(profile_dict['kukulka_list'])):
            ax.plot(profile_dict['kukulka_list'][counter], profile_dict['depth_bins'],
                    label=label_kukulka(parameters=profile_dict['parameter_kukulka'][counter]),
                    linestyle='--', color=colors[counter])
    lines, labels = ax.get_legend_handles_labels()

    # Plotting the diffusion curve
    if diffusion_curve:
        for counter in range(len(profile_dict['concentration_list'])):
            k_z, w_10, w_rise = profile_dict['parameter_concentrations'][counter]
            ax2 = diffusion_curve_axis(ax, ax_label_size, k_z, w_10, w_rise, profile_dict, diffusion_type, colors[counter])
            lines2, labels2 = ax2.get_legend_handles_labels()
            lines += lines2
            labels += labels2

    # Adding the legend
    ax.legend(lines, labels, fontsize=legend_size, loc='lower right')

    # Saving the figure
    plt.savefig(saving_filename_basic_profile(settings.figure_dir, selection, close_up, diffusion_type),
                bbox_inches='tight', dpi=600)


def timestep_comparison(k_z_list, w_10_list, w_rise_list, selection='k_z', close_up=None,
                        y_label='Depth (m)', x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(8, 8),
                        ax_label_size=16, legend_size=12, kukulka=True, model=True, single_select=0,
                        diffusion_type='Rouse', interval=1, boundary='Mixed', diffusion_curve=True):
    # Load the relevant data for the figure
    profile_dict = get_concentration_list(k_z_list, w_10_list, w_rise_list, selection, single_select, diffusion_type,
                                          all_timesteps=True, boundary=boundary)
    # Preparing for the actual plotting
    range_dict = get_axes_range(profile_dict['depth_bins'], profile_dict['concentration_list'])
    xmax, xmin = range_dict['max_count'], range_dict['min_count']
    if close_up == None:
        ymax, ymin = range_dict['max_depth'], range_dict['min_depth']
    else:
        # Allowing for easy close up for a specific part of the depth profile
        ymax, ymin = close_up
    ax_range = (xmax, xmin, ymax, ymin)
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray',
              'tab:olive', 'tab:cyan']
    # Creating the axis
    ax = base_figure(fig_size, ax_range, y_label, x_label, ax_label_size)

    # Plotting the modeled distribution
    if model:
        steps = 0
        for counter in range(0, len(profile_dict['concentration_list']), interval):
            ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                    label=label_time_step(steps, interval),
                    color=colors[steps % len(colors)])
            steps += 1
    # Plotting the distribution according to the Kukulka parametrization
    if kukulka:
        for counter in range(0, len(profile_dict['kukulka_list']), interval):
            ax.plot(profile_dict['kukulka_list'][counter], profile_dict['depth_bins'],
                    label=label_kukulka(parameters=profile_dict['parameter_kukulka'][counter]),
                    linestyle='--', color=colors[counter])
    lines, labels = ax.get_legend_handles_labels()

    # Plotting the diffusion curve
    if diffusion_curve:
        k_z, w_10, w_rise = profile_dict['parameter_concentrations'][0]
        ax2 = diffusion_curve_axis(ax, ax_label_size, k_z, w_10, w_rise, profile_dict, diffusion_type, 'black')
        lines2, labels2 = ax2.get_legend_handles_labels()
        lines += lines2
        labels += labels2

    # Adding the legend
    ax.legend(lines, labels, fontsize=legend_size, loc='lower right')
    # Saving the figure
    plt.savefig(saving_filename_time_step(settings.figure_dir, selection, close_up, diffusion_type),
                bbox_inches='tight', dpi=600)


def boundary_condition_comparison(k_z_list, w_10_list, w_rise_list, selection='k_z', close_up=None,
                        y_label='Depth (m)', x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(8, 8),
                        ax_label_size=16, legend_size=12, kukulka=True, model=True, single_select=0,
                        diffusion_type='KPP', interval=1, boundary='Mixed', diffusion_curve=True):
    # Load the relevant data for the figure
    profile_dict_mix = get_concentration_list(k_z_list, w_10_list, w_rise_list, selection, single_select,
                                              diffusion_type, boundary='Mixed')
    profile_dict_zero = get_concentration_list(k_z_list, w_10_list, w_rise_list, selection, single_select,
                                               diffusion_type, boundary='Zero_Ceiling')
    profile_dict_reflect = get_concentration_list(k_z_list, w_10_list, w_rise_list, selection, single_select,
                                                  diffusion_type, boundary='Reflect')
    profile_dict_reduce = get_concentration_list(k_z_list, w_10_list, w_rise_list, selection, single_select,
                                                  diffusion_type, boundary='Reduce_dt')
    profile_dict_markov1 = get_concentration_list(k_z_list, w_10_list, w_rise_list, selection, single_select,
                                                  diffusion_type, boundary='Reflect_Markov')

    # Preparing for the actual plotting
    range_dict = get_axes_range(profile_dict_mix['depth_bins'], profile_dict_mix['concentration_list'])
    xmax, xmin = range_dict['max_count'], range_dict['min_count']
    if close_up == None:
        ymax, ymin = range_dict['max_depth'], range_dict['min_depth']
    else:
        # Allowing for easy close up for a specific part of the depth profile
        ymax, ymin = close_up
    ax_range = (xmax, xmin, ymax, ymin)
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray',
              'tab:olive', 'tab:cyan']
    # Creating the axis
    ax = base_figure(fig_size, ax_range, y_label, x_label, ax_label_size)

    # First the mixed boundary layer
    for counter in range(len(profile_dict_mix['concentration_list'])):
        _, w_10, _ = profile_dict_mix['parameter_concentrations'][counter]
        ax.plot(profile_dict_mix['concentration_list'][counter], profile_dict_mix['depth_bins'],
                label=label_boundary(w_10, diffusion_type, 'Random Mixed Layer'),
                color=colors[0])

    # Next the zero ceiling boundary condition
    for counter in range(len(profile_dict_zero['concentration_list'])):
        _, w_10, _ = profile_dict_mix['parameter_concentrations'][counter]
        ax.plot(profile_dict_zero['concentration_list'][counter], profile_dict_zero['depth_bins'],
                label=label_boundary(w_10, diffusion_type, 'Zero Ceiling'),
                color=colors[1])

    # Then, the reflecting boundary condition
    for counter in range(len(profile_dict_reflect['concentration_list'])):
        _, w_10, _ = profile_dict_reflect['parameter_concentrations'][counter]
        ax.plot(profile_dict_reflect['concentration_list'][counter], profile_dict_reflect['depth_bins'],
                label=label_boundary(w_10, diffusion_type, 'Reflect'),
                color=colors[2])
    lines, labels = ax.get_legend_handles_labels()

    # Then, the reducing dt boundary condition
    for counter in range(len(profile_dict_reduce['concentration_list'])):
        _, w_10, _ = profile_dict_reduce['parameter_concentrations'][counter]
        ax.plot(profile_dict_reduce['concentration_list'][counter], profile_dict_reduce['depth_bins'],
                label=label_boundary(w_10, diffusion_type, 'Reduce dt'),
                color=colors[3])

    # Then, the Markov 1 reflecting boundary condition
    for counter in range(len(profile_dict_markov1['concentration_list'])):
        _, w_10, _ = profile_dict_markov1['parameter_concentrations'][counter]
        ax.plot(profile_dict_markov1['concentration_list'][counter], profile_dict_markov1['depth_bins'],
                label=label_boundary(w_10, diffusion_type, 'Markov-1, Reflect'),
                color=colors[4])
    lines, labels = ax.get_legend_handles_labels()


    # Plotting the diffusion curve
    if diffusion_curve:
        k_z, w_10, w_rise = profile_dict_mix['parameter_concentrations'][0]
        ax2 = diffusion_curve_axis(ax, ax_label_size, k_z, w_10, w_rise, profile_dict_mix, diffusion_type, 'black')
        lines2, labels2 = ax2.get_legend_handles_labels()
        lines += lines2
        labels += labels2

    # Adding the legend
    ax.legend(lines, labels, fontsize=legend_size, loc='lower right')
    # Saving the figure
    plt.savefig(saving_filename_boundary(settings.figure_dir, selection, close_up, diffusion_type),
                bbox_inches='tight', dpi=600)


def label_boundary(w_10, diffusion_type, boundary):
    if diffusion_type == 'Kukulka':
        return 'Kukulka et al. (2012),  w_10 = {}'.format(w_10) + ' m s$^{-1}$, ' + boundary
    elif diffusion_type == 'KPP':
        return r'KPP, w_10 = {}'.format(w_10) + 'm s$^{-1}$, MLD = '+'{} m'.format(settings.MLD) + ', ' + boundary


def diffusion_curve_axis(ax, ax_label_size, k_z, w_10, w_rise, profile_dict, diffusion_type, color):
    ax2 = ax.twiny()
    # Get the proper diffusion curve
    depth = profile_dict['depth_bins']
    profile = utils.get_vertical_diffusion_profile(w_10, k_z, depth * -1, diffusion_type)
    # X axis = Concentration axis
    ax2.set_xlabel(r'$K_z$ (m$^2$ s$^{-1}$)', fontsize=ax_label_size)
    ax2.set_xlim((0, 0.1))
    # The actual plotting
    ax2.plot(profile, depth, color=color, linestyle='dotted', label=label_diffusivity_profile(w_10, diffusion_type))
    return ax2


def label_time_step(steps, interval):
    t = steps * interval * settings.dt_out.seconds// 3600
    return 't = {} hours'.format(t)


def label_diffusivity_profile(w_10, diffusion_type):
    if diffusion_type == 'Kukulka':
        return 'Kukulka et al. (2012) $K_z$,  w_10 = {}'.format(w_10) + ' m s$^{-1}$'
    elif diffusion_type == 'KPP':
        return r'KPP $K_z$, w_10 = {}'.format(w_10) + 'm s$^{-1}$, MLD = '+'{} m'.format(settings.MLD)


def base_figure(fig_size, ax_range, y_label, x_label, ax_label_size):
    xmax, xmin, ymax, ymin = ax_range
    fig = plt.figure(figsize=fig_size)
    ax = fig.add_subplot(111)
    # Y axis = Depth axis
    ax.set_ylabel(y_label, fontsize=ax_label_size)
    ax.set_ylim((ymin, ymax))
    ax.tick_params(axis='both', labelsize=ax_label_size)
    # X axis = Concentration axis
    ax.set_xlabel(x_label, fontsize=ax_label_size)
    ax.set_xlim((xmin, xmax))
    return ax


def saving_filename_basic_profile(save_location, selection, close_up, diffusion_type):
    if close_up is None:
        return save_location + diffusion_type + '_Depth_profile_full_variable={}.png'.format(selection)
    else:
        ymax, ymin = close_up
        return save_location + diffusion_type + '_Depth_profile_max={}_min={}_variable={}.png'.format(ymax, ymin, selection)


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
    if selection == 'k_z':
        k_z, w_10, w_rise = parameters
        return r'$K_0 = $ ' + '{:.1E}'.format(k_z) + r' m$^2$ s$^{-1}$, $w_{10} = $ ' + '{}'.format(w_10) + \
               r' m s$^{-1}$, $w_r = $ ' + '{}'.format(w_rise) + ' m s$^{-1}$'
    elif selection == 'w_10':
        k_z, w_10, w_rise = parameters
        return r'$w_{10} = $ ' + '{}'.format(w_10) + r' m s$^{-1}$, $w_r = $ ' + '{}'.format(w_rise) + ' m s$^{-1}$'
    elif selection == 'w_rise':
        k_z, w_10, w_rise = parameters
        return r'$w_{10} = $ ' + '{}'.format(w_10) + r' m s$^{-1}$, $w_r = $ ' + '{}'.format(w_rise) + ' m s$^{-1}$'


def label_kukulka(parameters):
    w_10, w_rise = parameters
    lambda_c = np.round(1. / utils.determine_kukulka_e_length(w_10, w_rise), decimals=1)
    return 'Kukulka et al. (2012), $\lambda = $' + '{}'.format(lambda_c) + ' m'


def label_rouse(parameters):
    k_z, w_10, w_rise = parameters
    k = settings.vk  # von Karman constant
    u_shear = math.sqrt(utils.determine_tau(w_10, settings.rho_a) / settings.rho_w)
    Ro = math.fabs(w_rise) / (k * u_shear)  # Rouse value
    return r'Rouse (Boudreau & Hill, 2020), Ro = {:.1f}, $\epsilon_0$ = {:.1E}'.format(Ro, k_z) + r'm$^2$ s$^{-1}$'


def get_axes_range(depth_bins, concentrations):
    output_dict = {'max_depth': np.round(depth_bins.max(), -1), 'min_depth': np.round(depth_bins.min(), -1),
                   'max_count': 0, 'min_count': 0}
    for conc in concentrations:
        if conc.max() > output_dict['max_count']:
            output_dict['max_count'] = conc.max()
    return output_dict


def get_concentration_list(k_z_list, w_10_list, w_rise_list, selection, single_select, diffusion_type, output_step=-1,
                           all_timesteps=False, boundary='Mixed'):
    output_dic = {'concentration_list': [], 'kukulka_list': [], 'rouse_list': [], 'parameter_concentrations': [],
                  'parameter_kukulka': []}
    if selection == 'k_z':
        w_10_list = [w_10_list[single_select]]
        w_rise_list = [w_rise_list[single_select]]
    elif selection == 'w_10':
        k_z_list = [k_z_list[single_select]]
        w_rise_list = [w_rise_list[single_select]]
    elif selection == 'w_rise':
        k_z_list = [k_z_list[single_select]]
        w_10_list = [w_10_list[single_select]]

    for w_rise in w_rise_list:
        for w_10 in w_10_list:
            for k_z in k_z_list:
                # Loading the dictionary containing the concentrations
                input_dir = utils.load_obj(utils.get_concentration_output_name(k_z, w_10, w_rise, diffusion_type, boundary))
                # Selecting the timeslice of interest
                if output_step == -1:
                    concentration = [input_dir['last_time_slice']]
                else:
                    concentration = [output_step]
                if all_timesteps:
                    concentration = range(input_dir['last_time_slice'] + 1)
                for c in concentration:
                    output_dic['concentration_list'].append(input_dir[c] / input_dir[c].max())
                # Get the depth bins of the concentrations and getting the Rouse profile corresponding to the simulation
                output_dic['depth_bins'] = input_dir['bin_edges'][:-1] * -1
                output_dic['rouse_list'].append(get_Rouse_profile(k_z, w_10, w_rise, output_dic['depth_bins']))
                # Saving the physical parameters of the concentration distribution we have just loaded
                output_dic['parameter_concentrations'].append((k_z, w_10, w_rise))
            # Getting the Kukulka profile and the associated parameters
            output_dic['kukulka_list'].append(get_kukulka_profile(w_10, w_rise, output_dic['depth_bins']))
            output_dic['parameter_kukulka'].append((w_10, w_rise))
    return output_dic


def get_kukulka_profile(w_10, w_rise, depth_bins):
    folding_scale = utils.determine_kukulka_e_length(w_10, w_rise)
    kukulka_concentration = np.exp(folding_scale * depth_bins)
    return kukulka_concentration


def get_Rouse_profile(k_z, w_10, w_rise, depth_bins):
    k = settings.vk  # von Karman constant
    u_shear = math.sqrt(utils.determine_tau(w_10, settings.rho_a) / settings.rho_w)
    Ro = np.abs(w_rise) / (k * u_shear)  # Rouse value
    rouse_profile = np.power(1 + k * u_shear * np.abs(depth_bins) / k_z, -Ro)
    return rouse_profile
