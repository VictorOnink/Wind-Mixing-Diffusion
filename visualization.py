import settings
import numpy as np
import matplotlib.pyplot as plt
import utils
import math


def quick_plot(x, y):
    plt.plot(x[:-1], y, '.', c='r')
    plt.xlim((x.min(), x.max()))
    plt.ylim(y.min(), y.max())
    plt.show()


def basic_profile_figure(k_z_list, w_10_list, w_rise_list, selection='k_z', close_up=None,
                         y_label='Depth (m)', x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(8, 8),
                         ax_label_size=16, legend_size=12, rouse=True, kukulka=True, model=True, single_select=1,
                         output_step=-1, diffusion_type='Rouse'):
    # Load the relevant data for the figure
    profile_dict = get_concentration_list(k_z_list, w_10_list, w_rise_list, selection, single_select,
                                          output_step=output_step, diffusion_type=diffusion_type)
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
    # Plotting the distribution according to the theoretical Rouse profile
    if rouse:
        for counter in range(len(profile_dict['rouse_list'])):
            ax.plot(profile_dict['rouse_list'][counter], profile_dict['depth_bins'],
                    label=label_rouse(parameters=profile_dict['parameter_concentrations'][counter]),
                    linestyle='-.', color=colors[counter])

    # Adding the legend
    ax.legend(fontsize=legend_size, loc='lower right')

    # Saving the figure
    plt.savefig(saving_filename_basic_profile(settings.figure_dir, selection, close_up, diffusion_type),
                bbox_inches='tight', dpi=600)


def timestep_comparison(k_z_list, w_10_list, w_rise_list, selection='k_z', close_up=None,
                        y_label='Depth (m)', x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(8, 8),
                        ax_label_size=16, legend_size=12, rouse=True, kukulka=True, model=True, single_select=0,
                        time_range=1, diffusion_type='Rouse', interval = 1):
    # Load the relevant data for the figure
    profile_dict = get_concentration_list(k_z_list, w_10_list, w_rise_list, selection, single_select, diffusion_type,
                                          all_timesteps=True)
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
        for counter in range(0, len(profile_dict['concentration_list']), interval):
            ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                    label=label_profile(selection, parameters=profile_dict['parameter_concentrations'][0]),
                    color=colors[counter%len(colors)])
    # Plotting the distribution according to the Kukulka parametrization
    if kukulka:
        for counter in range(0, len(profile_dict['kukulka_list']), interval):
            ax.plot(profile_dict['kukulka_list'][counter], profile_dict['depth_bins'],
                    label=label_kukulka(parameters=profile_dict['parameter_kukulka'][counter]),
                    linestyle='--', color=colors[counter])
    # Plotting the distribution according to the theoretical Rouse profile
    if rouse:
        for counter in range(0, len(profile_dict['rouse_list']), interval):
            ax.plot(profile_dict['rouse_list'][counter], profile_dict['depth_bins'],
                    label=label_rouse(parameters=profile_dict['parameter_concentrations'][counter]),
                    linestyle='-.', color=colors[counter])

    # Adding the legend
    # ax.legend(fontsize=legend_size, loc='lower right')
    plt.show()


def saving_filename_basic_profile(save_location, selection, close_up, diffusion_type):
    if close_up == None:
        return save_location + diffusion_type + '_Depth_profile_full_variable={}.png'.format(selection)
    else:
        ymax, ymin = close_up
        return save_location + diffusion_type + '_Depth_profile_max={}_min={}_variable={}.png'.format(ymax, ymin, selection)


def label_profile(selection, parameters):
    if selection == 'k_z':
        k_z, w_10, w_rise = parameters
        return r'$K_0 = $ ' + '{:.1E}'.format(k_z) + r' m$^2$ s$^{-1}$, $w_{10} = $ ' + '{}'.format(w_10) + \
               r' m s$^{-1}$, $w_r = $ ' + '{}'.format(w_rise) + ' m s$^{-1}$'
    elif selection == 'w_10':
        k_z, w_10, w_rise = parameters
        return r'$K_0 = $ ' + '{:.1E}'.format(k_z) + r' m$^2$ s$^{-1}$, $w_{10} = $ ' + '{}'.format(w_10) + \
               r' m s$^{-1}$, $w_r = $ ' + '{}'.format(w_rise) + ' m s$^{-1}$'
    elif selection == 'w_rise':
        k_z, w_10, w_rise = parameters
        return r'$K_0 = $ ' + '{:.1E}'.format(k_z) + r' m$^2$ s$^{-1}$, $w_{10} = $ ' + '{}'.format(w_10) + \
               r' m s$^{-1}$, $w_r = $ ' + '{}'.format(w_rise) + ' m s$^{-1}$'


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


def get_axes_range(depth_bins, concentrations):
    output_dict = {'max_depth': np.round(depth_bins.max(), -1), 'min_depth': np.round(depth_bins.min(), -1),
                   'max_count': 0, 'min_count': 0}
    for conc in concentrations:
        if conc.max() > output_dict['max_count']:
            output_dict['max_count'] = conc.max()
    return output_dict


def get_concentration_list(k_z_list, w_10_list, w_rise_list, selection, single_select, diffusion_type, output_step=-1,
                           all_timesteps=False):
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
                input_dir = utils.load_obj(utils.get_concentration_output_name(k_z, w_10, w_rise, diffusion_type))
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
