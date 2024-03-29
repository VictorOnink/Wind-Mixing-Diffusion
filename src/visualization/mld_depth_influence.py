import matplotlib.pyplot as plt
import utils, settings
from visualization import utils_visualization as utils_v
import numpy as np


def mld_depth_influence(w_rise_list, MLD_list, alpha_list, selection='w_rise', output_step=-1, single_select=0,
                        y_label='Depth/MLD', close_up=None, beaufort=5,
                        x_label=r'Normalised Plastic Counts ($n/n_0$)', fig_size=(16, 8), ax_label_size=16,
                        legend_size=12, diffusion_type='KPP', boundary='Reflect', alpha=0.3):
    """
    Figure looking at the influence of the MLD on the diffusion profile, where everything is plotted with normalized
    depths (relative to the MLD)
    :param w_rise_list: list of rise velocities
    :param MLD_list: list of MLD depths
    :param alpha_list: list of alpha values for M-1 simulation
    :param selection: selection criteria for loading the parcels concentration data
    :param output_step: which time index is being plotted, with the default being the last of the simulation
    :param single_select: selection criteria related to 'selection'
    :param y_label: label of the y axis
    :param close_up: setting the limits of the y axis as (max, min0
    :param beaufort: wind force selecting which profiles to plot
    :param x_label: label of the x axis
    :param fig_size: size of the figure
    :param ax_label_size: fontsize of the axes labels
    :param legend_size: fontsize of the legend
    :param diffusion_type: type of diffusion, either KPP or SWB
    :param boundary: which boundary condition, and whether M-0 or M-1
    :param alpha: transparancy of the field data markers
    :return: 
    """
    correction = settings.MLD
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=True)

    # Selecting which model data we want to plot based on the diffusion type
    swb, kpp, artificial = utils_v.boolean_diff_type(diffusion_type)

    # Get the base figure axis
    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size)

    # Plotting the field data points
    wind_range = utils.beaufort_limits()[beaufort]
    _, _ = utils_v.add_observations(ax, norm_depth=True, alpha=alpha, wind_range=wind_range)

    line_style = ['-', '--', '-.']
    for count_mld, mld in enumerate(MLD_list):
        # Plotting the distribution according to the SWB parametrization
        if swb:
            profile_dict = utils_v.get_concentration_list([np.mean(wind_range)], w_rise_list, selection,
                                                          single_select, alpha_list=alpha_list,
                                                          output_step=output_step, diffusion_type='SWB',
                                                          boundary=boundary, mld=mld)
            for counter in range(len(profile_dict['concentration_list'])):
                ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'] / correction,
                        label=label_MLD_Comparison(parameters=profile_dict['parameter_SWB'][counter],
                                                   mld=mld, diffusion_type='SWB'),
                        linestyle=line_style[counter], color=utils_v.return_color(count_mld))

        # Plotting the distribution according to the KPP parametrization
        if kpp:
            profile_dict = utils_v.get_concentration_list([np.mean(wind_range)], w_rise_list, selection,
                                                          single_select, alpha_list=alpha_list,
                                                          output_step=output_step, diffusion_type='KPP',
                                                          boundary=boundary, mld=mld)
            for counter in range(len(profile_dict['concentration_list'])):
                ax.plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'] / correction,
                        label=label_MLD_Comparison(parameters=profile_dict['parameter_SWB'][counter],
                                                   mld=mld, diffusion_type='KPP'),
                        linestyle=line_style[counter], color=utils_v.return_color(count_mld))
    # Adding the legend
    ax.legend(fontsize=legend_size, loc='lower right')

    ax.set_title(r'u$_{10}$' + '={}-{}'.format(*wind_range) + r' m s$^{-1}$', fontsize=ax_label_size)

    plt.savefig(mld_comparison_name(diffusion_type, boundary, beaufort=beaufort, close_up=close_up),
                bbox_inches='tight')


def mld_comparison_name(diffusion_type, boundary, beaufort, close_up=None, output_type='.png'):
    """ Creating the name of the figure"""
    diff_dict = {'SWB': 'SWB', 'KPP': 'KPP', 'all': 'SWB_KPP'}
    figure_name = settings.figure_dir + 'norm_comparison_{}_{}_Bft{}'.format(diff_dict[diffusion_type], boundary,
                                                                             beaufort)
    if close_up is not None:
        max, min = close_up
        figure_name += '_max_{}_min_{}'.format(max, min)
    return figure_name + output_type


def label_MLD_Comparison(parameters, diffusion_type, mld=settings.MLD):
    """ Creating the labels of the plot """
    w_10, w_rise = parameters
    w_rise = np.abs(w_rise)
    if diffusion_type == 'KPP':
        return r'KPP, u$_{10}$ '+'= {:.2f}'.format(w_10) + 'm s$^{-1}$,' + 'w$_{rise}$ '+'= {}'.format(w_rise) + \
               'm s$^{-1}$, MLD = ' + '{:.1f} m'.format(mld)
    elif diffusion_type == 'SWB':
        return r'SWB, u$_{10}$ '+'= {:.2f}'.format(w_10) + 'm s$^{-1}$,' + 'w$_{rise}$ '+'= {}'.format(w_rise) + \
               'm s$^{-1}$, MLD = ' + '{:.1f} m'.format(mld)