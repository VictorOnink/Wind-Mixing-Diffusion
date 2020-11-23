import settings
import numpy as np
import matplotlib.pyplot as plt
import utils
from decimal import Decimal


def quick_plot(x, y):
    plt.plot(x[:-1], y, '.', c='r')
    plt.xlim((x.min(), x.max()))
    plt.ylim(y.min(), y.max())
    plt.show()


def basic_profile_figure(k_z_list, w_10_list, w_rise_list, selection='k_z', close_up=None,
                         y_label='Depth (m)', x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(8, 8),
                         ax_label_size=16, legend_size=12):
    # Load the relevant data for the figure
    conc_list, depth_bins, kukulka_profiles, par_conc, par_kukulka = get_concentration_list(k_z_list, w_10_list,
                                                                                            w_rise_list, selection)
    # Preparing for the actual plotting
    range_dict = get_axes_range(depth_bins, conc_list)
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
    for counter in range(len(conc_list)):
        ax.plot(conc_list[counter], depth_bins, label=label_profile(selection, parameters=par_conc[counter]),
                color=colors[counter])
    # Plotting the distribution according to the Kukulka parametrization
    for counter in range(len(kukulka_profiles)):
        ax.plot(kukulka_profiles[counter], depth_bins, label=label_kukulka(parameters=par_kukulka[counter]),
                linestyle='--',color=colors[counter])
    ax.legend(fontsize=legend_size, loc='lower right')
    plt.savefig(saving_filename(settings.figure_dir, selection, close_up), bbox_inches='tight', dpi=600)


def saving_filename(save_location, selection, close_up):
    if close_up == None:
        return save_location + 'Depth_profile_full_variable={}.png'.format(selection)
    else:
        ymax, ymin = close_up
        return save_location + 'Depth_profile_max={}_min={}_variable={}.png'.format(ymax, ymin, selection)


def label_profile(selection, parameters):
    if selection == 'k_z':
        k_z, w_10, w_rise = parameters
        return r'$K_0 = $ ' + '{:.1E}'.format(k_z) + r' m$^2$ s$^{-1}$, $w_{10} = $ ' + '{}'.format(w_10) + \
               r' m s$^{-1}$, $w_r = $ '+'{}'.format(w_rise) + ' m s$^{-1}$'
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
    lambda_c = np.round(1. / utils.determine_kukulka_e_length(w_10, w_rise),decimals=1)
    # return 'Kukulka et al. (2012), $w_{10} = $ ' + '{}'.format(w_10) + r' m s$^{-1}$, $w_r = $ ' + '{}'.format(w_rise) \
    #        + r' m s$^{-1}$, $\lambda = $' + '{}'.format(lambda_c) + ' m'
    return 'Kukulka et al. (2012), $\lambda = $' + '{}'.format(lambda_c) + ' m'


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


def get_concentration_list(k_z_list, w_10_list, w_rise_list, selection):
    concentration_list = []
    kukulka_list = []
    parameter_concentration = []
    parameter_kukulka = []
    if selection == 'k_z':
        w_10_list = [w_10_list[1]]
        w_rise_list = [w_rise_list[1]]
    elif selection == 'w_10':
        k_z_list = [k_z_list[1]]
        w_rise_list = [w_rise_list[1]]
    elif selection == 'w_rise':
        k_z_list = [k_z_list[1]]
        w_10_list = [w_10_list[1]]
    for w_rise in w_rise_list:
        for w_10 in w_10_list:
            for k_z in k_z_list:
                input_dir = utils.load_obj(utils.get_concentration_output_name(k_z, w_10, w_rise))
                depth_bins = input_dir['bin_edges'][:-1] * -1
                concentration_list.append(input_dir['concentration'] / input_dir['concentration'].max())
                parameter_concentration.append((k_z, w_10, w_rise))
            kukulka_list.append(_get_kukulka_profile(w_10, w_rise, depth_bins))
            parameter_kukulka.append((w_10, w_rise))

    return concentration_list, depth_bins, kukulka_list, parameter_concentration, parameter_kukulka


def _get_kukulka_profile(w_10, w_rise, depth_bins):
    folding_scale = utils.determine_kukulka_e_length(w_10, w_rise)
    kukulka_concentration = np.exp(folding_scale * depth_bins)
    return kukulka_concentration
