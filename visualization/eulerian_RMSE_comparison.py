import settings
import matplotlib.pyplot as plt
import numpy as np
import analysis


def eulerian_RMSE_comparison(x_label=r'$u_{10}$ (m s$^{-1}$)', y_label=r'RMSE', fig_size=(8, 8),
                             ax_label_size=16, legend_size=12):
    """
    plotting the RMSE differences between the field observations and the Eulerian model concentration profiles
    :param x_label: the x axis label
    :param y_label: the y axis label
    :param fig_size: the figure size
    :param ax_label_size:
    :param legend_size:
    :return:
    """
    # Setting the wind speeds and rise velocities for which we want to make the plot
    w_10 = [0.85, 2.4, 4.35, 6.65, 9.3]
    w_r = [0.03, 0.003]
    # Setting the diffusion types
    diffusion = ['KPP', 'Kukulka']
    # We add a slight offset to the x-axis positions so that the markers from the two diffusion types are not plotted
    # on top of each other
    diffusion_offset = {'KPP': 0.1, 'Kukulka': -0.1}
    # Selecting the marker type according which form of diffusion it is:
    marker_type = {'KPP': 'o', 'Kukulka': 'X', 'SWB': 'X'}
    # Selecting the marker color according to which rise velocity
    marker_color = {0.03: 'tab:blue', 0.003: 'tab:red', 0.0003: 'tab:green'}

    # Looping through the simulations, and retrieving the RMSE values for them
    point_list = []
    for index_w10, wind in enumerate(w_10):
        for rise in w_r:
            for diffusion_type in diffusion:
                RMSE = analysis.determine_RMSE_eulerian(wind, rise, diffusion_type)
                plot_tuple = RMSE, index_w10 + 1 + diffusion_offset[diffusion_type], marker_type[diffusion_type], marker_color[rise]
                point_list.append(plot_tuple)

    # Creating the axis
    fig = plt.figure(figsize=fig_size)
    ax = fig.add_subplot(111)
    # X axis = wind speed axis
    ax.set_xlabel(x_label, fontsize=ax_label_size)
    ax.set_xlim((0, 6))
    ax.tick_params(axis='both', labelsize=ax_label_size)
    # Y axis = RMSE axis
    ax.set_ylabel(y_label, fontsize=ax_label_size)
    ax.set_ylim((0, 0.6))

    # Now, plotting the points
    for point in point_list:
        RMSE, index_w10, marker, color = point
        ax.plot(index_w10, RMSE, color=color, marker=marker, alpha=0.7, markersize=10)

    # Now, altering the Y axis to list the wind speeds instead of the simple labels 1 - 5
    ax.set_xticks(range(7))
    ax.set_xticklabels(['', 0.85, 2.40, 4.35, 6.65, 9.30, ''])

    # Next, adding a legend to explain the color scheme and the marker type
    label_marker = ['KPP', 'SWB']
    # The marker type indicates the type of diffusion
    marker = [plt.plot([], [], c='k', markersize=10, marker=marker_type[label], label=label, linestyle='')[0] for label in
              label_marker]

    # The color indicates the rise velocity
    color = [plt.plot([], [], c=marker_color[rise], markersize=10, marker='o', label=label_color(rise), linestyle='')[0] for rise in
              w_r]

    ax.legend(handles=marker + color, fontsize=legend_size, loc='lower right')

    plt.savefig(settings.figure_dir + 'model_evaluation_eulerian.png', bbox_inches='tight', dpi=600)


def label_color(rise):
    """ defining the label for the given rise velocity"""
    return r'$w_r$ = {}'.format(np.abs(rise)) + ' m s$^{-1}$'
