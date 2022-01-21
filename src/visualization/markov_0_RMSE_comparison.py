import matplotlib.pyplot as plt
import numpy as np
import analysis, settings
from visualization import utils_visualization as utils_v


def markov_0_RMSE_comparison(x_label=r'$u_{10}$ (m s$^{-1}$)', y_label=r'RMSE', fig_size=(16, 8),
                             ax_label_size=16, legend_size=12, boundary='Ceiling', wave_roughness=False):
    """
    A figure showing the RMSE of M-0 simulations relative to the field data for the various wind conditions
    :param x_label: x axis label
    :param y_label: y axis label
    :param fig_size: figure size
    :param ax_label_size: axis label fontsize
    :param legend_size: legend fontsize
    :param wave_roughness: if True, have surface roughness be wave height dependent
    :return:
    """
    # Setting the wind speeds and rise velocities that we want to plot, along with the theta values
    w_10 = [0.85, 2.4, 4.35, 6.65, 9.3]
    theta_list = [1.0, 2.0, 3.0, 4.0, 5.0]
    gamma_list = [0.5, 1.0, 1.5, 2.0]
    # Setting the type of diffusion, and a small offset for the x axis so that the markers are not plotted on top of
    # each other
    diffusion_offset = {'KPP': 0.1, 'SWB': -0.1}

    # Looping through the simulations, and retrieving the RMSE values for them
    point_list_high_kpp, point_list_low_kpp = [], []
    point_list_high_SWB, point_list_low_SWB = [], []
    for index_w10, wind in enumerate(w_10):
        for index_theta, theta in enumerate(theta_list):
            color = utils_v.discrete_color_from_cmap(index=index_theta, subdivisions=theta_list.__len__())
            RMSE = analysis.determine_RMSE(wind, -0.03, 'KPP', boundary, alpha=0.0, output=True,
                                           wave_roughness=wave_roughness, theta=theta)
            plot_tuple = RMSE, index_w10 + 1 + diffusion_offset['KPP'], 'o', color
            point_list_high_kpp.append(plot_tuple)

            RMSE = analysis.determine_RMSE(wind, -0.003, 'KPP', boundary, alpha=0.0, output=True,
                                           wave_roughness=wave_roughness, theta=theta)
            plot_tuple = RMSE, index_w10 + 1 + diffusion_offset['KPP'], 'o', color
            point_list_low_kpp.append(plot_tuple)
    for index_w10, wind in enumerate(w_10):
        for index_gamma, gamma in enumerate(gamma_list):
            color = utils_v.discrete_color_from_cmap(index=index_gamma, subdivisions=gamma_list.__len__())
            RMSE = analysis.determine_RMSE(wind, -0.03, 'SWB', boundary, alpha=0.0, output=True,
                                           wave_roughness=wave_roughness, gamma=gamma)
            plot_tuple = RMSE, index_w10 + 1 + diffusion_offset['SWB'], 'X', color
            point_list_high_SWB.append(plot_tuple)

            RMSE = analysis.determine_RMSE(wind, -0.003, 'SWB', boundary, alpha=0.0, output=True,
                                           wave_roughness=wave_roughness, gamma=gamma)
            plot_tuple = RMSE, index_w10 + 1 + diffusion_offset['SWB'], 'X', color
            point_list_low_SWB.append(plot_tuple)

    # Creating the axis
    fig = plt.figure(figsize=fig_size)
    ax1 = fig.add_subplot(121)
    # X axis = wind speed axis
    ax1.set_xlabel(x_label, fontsize=ax_label_size)
    ax1.set_xlim((0, 6))
    ax1.tick_params(axis='both', labelsize=ax_label_size)
    # Y axis = RMSE axis
    ax1.set_ylabel(y_label, fontsize=ax_label_size)
    ax1.set_ylim((0, 0.13))

    ax2 = fig.add_subplot(122)
    # X axis = wind speed axis
    ax2.set_xlabel(x_label, fontsize=ax_label_size)
    ax2.set_xlim((0, 6))
    ax2.tick_params(axis='both', labelsize=ax_label_size, labelleft=False)
    # Y axis = RMSE axis
    # ax2.set_ylabel(y_label, fontsize=ax_label_size)
    ax2.set_ylim((0, 0.13))

    # Now, plotting the points
    for point in point_list_high_kpp:
        RMSE, index_w10, marker, color = point
        color_face = (color[0], color[1], color[2], 0.5)
        ax1.plot(index_w10, RMSE, mec=color, marker=marker, markersize=10, mfc=color_face)
    for point in point_list_high_SWB:
        RMSE, index_w10, marker, color = point
        color_face = (color[0], color[1], color[2], 0.5)
        ax1.plot(index_w10, RMSE, mec=color, marker=marker, markersize=10, mfc=color_face)

    for point in point_list_low_kpp:
        RMSE, index_w10, marker, color = point
        color_face = (color[0], color[1], color[2], 0.5)
        ax2.plot(index_w10, RMSE, mec=color, marker=marker, markersize=10, mfc=color_face)
    for point in point_list_low_SWB:
        RMSE, index_w10, marker, color = point
        color_face = (color[0], color[1], color[2], 0.5)
        ax2.plot(index_w10, RMSE, mec=color, marker=marker, markersize=10, mfc=color_face)

    # Now, altering the Y axis to list the wind speeds instead of the simple labels 1 - 5
    for ax in [ax1, ax2]:
        ax.set_xticks(range(7))
        ax.set_xticklabels(['', 0.85, 2.40, 4.35, 6.65, 9.30, ''])

    # Next, adding a legend to explain the color scheme and the marker type
    # The marker color indicates the theta value
    swb = [plt.plot([], [], c=utils_v.discrete_color_from_cmap(index=index_gamma, subdivisions=gamma_list.__len__()),
                    markersize=10, marker='X', label=label_gamma(gamma), linestyle='')[0]
           for index_gamma, gamma in enumerate(gamma_list)]
    kpp = [plt.plot([], [], c=utils_v.discrete_color_from_cmap(index=index_theta, subdivisions=theta_list.__len__()),
                    markersize=10, marker='o', label=label_theta(theta), linestyle='')[0]
           for index_theta, theta in enumerate(theta_list)]

    # Creating a legend
    ax2.legend(handles=swb + kpp, fontsize=legend_size, loc='upper right')

    plt.savefig(settings.figure_dir + 'model_evaluation_markov_0.png', bbox_inches='tight', dpi=600)


def label_theta(theta):
    return r'KPP, $\theta=$' + '{:.1f}'.format(theta)


def label_gamma(gamma):
    return r'SWB, $\gamma=$' + '{:.1f}'.format(gamma)

