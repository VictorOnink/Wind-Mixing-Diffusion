import matplotlib.pyplot as plt
import numpy as np
import analysis, settings
from visualization import utils_visualization as utils_v


def markov_1_RMSE_comparison(x_label=r'$u_{10}$ (m s$^{-1}$)', y_label=r'RMSE', fig_size=(16, 8),
                             ax_label_size=16, legend_size=11):
    """
    A figure showing the RMSE of M-0 and M-1 simulations relative to the field data for the various wind conditions
    :param x_label: x axis label
    :param y_label: y axis label
    :param fig_size: figure size
    :param ax_label_size: axis label fontsize
    :param legend_size: legend fontsize
    :return:
    """
    # Setting the wind speeds, rise velocities and alpha values that we want to plot
    w_10 = [0.85, 2.4, 4.35, 6.65, 9.3]
    w_r = [-0.03, -0.003]
    alpha = [0.0, 0.1, 0.3, 0.5, 0.7, 0.95]
    # Setting a small offset so that markers for different rise velocities aren't plotted on top of each other
    rise_offset = {-0.03: 0.15, -0.003: -0.15}
    # Selecting the marker type according to the rise velocity
    marker_type = {-0.03: 'o', -0.003: 'X'}

    # Looping through the KPP simulations, and retrieving the RMSE values for them. First, for the M-0 simulations
    point_list_KPP = []
    for index_w10, wind in enumerate(w_10):
        for rise in w_r:
            RMSE = analysis.determine_RMSE(wind, rise, 'KPP', 'Ceiling', alpha=0.0, output=True)
            plot_tuple = RMSE, index_w10 + 1, marker_type[rise], utils_v.return_color(0)
            point_list_KPP.append(plot_tuple)
    # Then for the M-1 simulations
    for index_w10, wind in enumerate(w_10):
        for rise in w_r:
            for index_a, a in enumerate(alpha):
                RMSE = analysis.determine_RMSE(wind, rise, 'KPP', 'Ceiling_Markov', alpha=a, output=True)
                plot_tuple = RMSE, index_w10 + 1 + rise_offset[rise], marker_type[rise], utils_v.return_color(index_a + 1)
                point_list_KPP.append(plot_tuple)

    # Looping through the SWB simulations, and retrieving the RMSE values for them, first for the M-0 simulations
    point_list_Kukulka = []
    for index_w10, wind in enumerate(w_10):
        for rise in w_r:
            RMSE = analysis.determine_RMSE(wind, rise, 'SWB', 'Ceiling', alpha=0.0, output=True)
            plot_tuple = RMSE, index_w10 + 1, marker_type[rise], utils_v.return_color(0)
            point_list_Kukulka.append(plot_tuple)
    # And then for M-1 simulations
    for index_w10, wind in enumerate(w_10):
        for rise in w_r:
            for index_a, a in enumerate(alpha):
                RMSE = analysis.determine_RMSE(wind, rise, 'SWB', 'Ceiling_Markov', alpha=a, output=True)
                plot_tuple = RMSE, index_w10 + 1 + rise_offset[rise], marker_type[rise], utils_v.return_color(index_a + 1)
                point_list_Kukulka.append(plot_tuple)

    # Creating the axis
    fig = plt.figure(figsize=fig_size)
    # Adding the axis for KPP
    ax = fig.add_subplot(121)
    ax.set_xlabel(x_label, fontsize=ax_label_size)
    ax.set_xlim((0, 6))
    ax.tick_params(axis='both', labelsize=ax_label_size)
    ax.set_ylabel(y_label, fontsize=ax_label_size)
    ax.set_ylim((0, 0.6))
    # Adding the axis for SWB
    ax2 = fig.add_subplot(122)
    ax2.set_xlim((0, 6))
    ax2.set_xlabel(x_label, fontsize=ax_label_size)
    ax2.tick_params(axis='both', labelsize=ax_label_size)
    ax2.tick_params(labelleft=False)
    # X axis = Concentration axis
    ax2.set_ylim((0, 0.6))

    ax.set_title(r'(a) KPP', fontsize=ax_label_size)
    ax2.set_title(r'(b) SWB', fontsize=ax_label_size)

    # Now, plotting the points
    for point in point_list_KPP:
        RMSE, index_w10, marker, color = point
        ax.plot(index_w10, RMSE, color=color, marker=marker, alpha=0.7, markersize=10, mfc=None)
    for point in point_list_Kukulka:
        RMSE, index_w10, marker, color = point
        ax2.plot(index_w10, RMSE, color=color, marker=marker, alpha=0.7, markersize=10, mfc=None)

    # Now, altering the Y axis to list the wind speeds instead of the simple labels 1 - 5
    ax.set_xticks(range(7))
    ax.set_xticklabels(['', 0.85, 2.40, 4.35, 6.65, 9.30, ''])
    ax2.set_xticks(range(7))
    ax2.set_xticklabels(['', 0.85, 2.40, 4.35, 6.65, 9.30, ''])

    # Next, adding a legend to explain the color scheme and the marker type
    # Showing the marker type according to the rise velocity
    marker = [plt.plot([], [], c='k', markersize=10, marker=marker_type[rise], label=label_marker(rise), linestyle='')[0] for rise in
              w_r]
    # Showing the color according to M-0/M-1 with alpha values
    markov0 = [plt.plot([], [], c=utils_v.return_color(0), markersize=10, marker='o', label='M0', linestyle='')[0]]
    markov1 = [plt.plot([], [], c=utils_v.return_color(ind + 1), markersize=10, marker='o',
                        label=r'M1 - $\alpha = $' + '{}'.format(a), linestyle='')[0] for ind, a in
               enumerate(alpha)]
    # Adding the legend
    ax2.legend(handles=marker + markov0 + markov1, fontsize=legend_size, loc='upper right')
    # Saving the figure
    plt.savefig(settings.figure_dir + 'model_evaluation_markov_1.png', bbox_inches='tight', dpi=600)


def label_marker(rise):
    """ Setting the figure label based on the rise velocity"""
    return r'$w_{rise}$' + ' = {}'.format(np.abs(rise)) + ' m s$^{-1}$'
