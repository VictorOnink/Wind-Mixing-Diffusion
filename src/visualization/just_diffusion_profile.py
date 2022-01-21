import matplotlib.pyplot as plt
import utils, settings
from visualization import utils_visualization as utils_v
import numpy as np

from visualization.utils_visualization import label_diffusivity_profile


def just_diffusion_profile(w_10_list, y_label='Depth (m)', x_label=r'$K_z$ (m$^2$ s$^{-1}$)',
                           fig_size=(10, 10), ax_label_size=16, legend_size=12, theta=1,
                           wave_roughness=False, with_theta=True, kpp=True, swb=True, with_gamma=True):
    """
    A simple figure just showing the diffusion profiles for KPP and SWB diffusion with different wind conditions
    :param w_10_list: list of the wind speeds
    :param y_label: the y axis label
    :param x_label: the x axis label
    :param fig_size: the size of the figure
    :param ax_label_size: the fontsize of the axis labels
    :param legend_size: the fontsize of the legend
    :param theta: Langmuir circulation amplification factor
    :param wave_roughness: if True, set the roughness scale to be the waveheight for KPP diffusion
    :param with_theta: if True, plot the KPP diffusion profile for theta=5
    :param kpp: if True, plot KPP profiles
    :param swb: if True, plot SWB profiles
    :param with_gamma: if True, plot influence of SWB
    :return:
    """
    # Setting the range of the x and y axis
    ymax, ymin = 0, -1 * (settings.MLD + 10)
    ax_range = (1e-1, 1e-5, ymax, ymin)
    depth = np.linspace(ymax, np.abs(ymin), 1000)

    # Creating the axis
    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size, log_xscale=True)

    # Plotting the diffusion profile according to the SWB approach
    if swb:
        for count, w_10 in enumerate(w_10_list):
            profile = utils.get_vertical_diffusion_profile(w_10, depth, 'SWB', theta=theta, wave_roughness=False)
            ax.plot(profile, -1 * depth, color=utils_v.discrete_color_from_cmap(count, len(w_10_list)), linestyle='-',
                    label=label_diffusivity_profile(w_10, 'SWB'), linewidth=3)
        if with_gamma:
            for count, w_10 in enumerate(w_10_list):
                profile = utils.get_vertical_diffusion_profile(w_10, depth, 'SWB', theta=theta, wave_roughness=False,
                                                               gamma=2.0)
                ax.plot(profile, -1 * depth, color=utils_v.discrete_color_from_cmap(count, len(w_10_list)),
                        linestyle='-', label=label_diffusivity_profile(w_10, 'SWB', gamma=2.0))

    # Plotting the diffusion profile according the KPP approach with theta = 1
    if kpp:
        for count, w_10 in enumerate(w_10_list):
            profile = utils.get_vertical_diffusion_profile(w_10, depth, 'KPP', theta=1, wave_roughness=False)
            ax.plot(profile, -1 * depth, color=utils_v.discrete_color_from_cmap(count, len(w_10_list)),
                    linestyle='--', label=label_diffusivity_profile(w_10, 'KPP'))

    # Plotting the diffusion profile according the KPP approach with theta = 5
    if with_theta:
        for count, w_10 in enumerate(w_10_list):
            profile = utils.get_vertical_diffusion_profile(w_10, depth, 'KPP', theta=5, wave_roughness=False)
            ax.plot(profile, -1 * depth, color=utils_v.discrete_color_from_cmap(count, len(w_10_list)),
                    linestyle='dotted', label=label_diffusivity_profile(w_10, 'KPP'))

    # Plotting the KPP diffusion profile with wave surface roughness
    if wave_roughness:
        for count, w_10 in enumerate(w_10_list):
            profile = utils.get_vertical_diffusion_profile(w_10, depth, 'KPP', theta=1, wave_roughness=True)
            ax.plot(profile, -1 * depth, color=utils_v.discrete_color_from_cmap(count, len(w_10_list)),
                    linestyle='-.', label=label_diffusivity_profile(w_10, 'KPP'))

    # Creating a legend
    color_labels = [plt.plot([], [], c=utils_v.discrete_color_from_cmap(index_size, subdivisions=w_10_list.__len__()),
                                label=wind_label(w_10), linestyle='-')[0] for index_size, w_10 in enumerate(w_10_list)]
    line_labels = []
    if swb:
        line_labels += [plt.plot([], [], c='k', label=r'SWB, $\gamma = 1.0$', linestyle='-', linewidth=3.0)[0]]
        if with_gamma:
            line_labels += [plt.plot([], [], c='k', label=r'SWB, $\gamma = 2.0$', linestyle='-')[0]]
    if kpp:
        line_labels += [plt.plot([], [], c='k', label=r'KPP, MLD = 20 m, $\theta=1$', linestyle='--')[0]]
        if with_theta:
            line_labels += [plt.plot([], [], c='k', label=r'KPP, MLD = 20 m, $\theta=5$', linestyle='dotted')[0]]
        if wave_roughness:
            line_labels += [plt.plot([], [], c='k', label=r'KPP, MLD = 20 m, $z_0=0.1\times H_s$', linestyle='-.')[0]]

    ax.legend(handles=color_labels + line_labels, fontsize=legend_size, loc='lower right', ncol=2)

    plt.savefig(file_name(theta, wave_roughness), bbox_inches='tight')


def wind_label(w_10):
    return 'u$_{10}$' + ' = {:.2f}'.format(w_10) + ' m s$^{-1}$'


def file_name(theta, wave_roughness, file_type='.png'):
    filename = settings.figure_dir + 'diffusion_profile_MLD=20_theta={}'.format(theta)
    if wave_roughness:
        filename += '_wave_roughness'
    return filename + file_type

