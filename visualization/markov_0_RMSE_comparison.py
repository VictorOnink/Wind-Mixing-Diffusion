import settings
import matplotlib.pyplot as plt
import numpy as np
import analysis


def markov_0_RMSE_comparison(y_label=r'$u_{10}$ (m s$^{-1}$)', x_label=r'RMSE', fig_size=(8, 8),
                             ax_label_size=16, legend_size=12):
    w_10 = [0.85, 2.4, 4.35, 6.65, 9.3]
    w_r = [-0.03, -0.003]
    boundary = 'Reflect'
    diffusion = ['KPP', 'Kukulka']
    diffusion_offset = {'KPP': 0.1, 'Kukulka': -0.1}
    # Selecting the marker type according which form of diffusion it is:
    marker_type = {'KPP': 'o', 'Kukulka': 'X', 'SWB': 'X'}
    # Selecting the marker color according to which rise velocity
    marker_color = {-0.03: 'tab:blue', -0.003: 'tab:red', -0.0003: 'tab:green'}

    # Looping through the simulations, and retrieving the RMSE values for them
    point_list = []
    for index_w10, wind in enumerate(w_10):
        for rise in w_r:
            for diffusion_type in diffusion:
                RMSE = analysis.determine_RMSE(wind, rise, diffusion_type, boundary, alpha=0.0, output=True)
                plot_tuple = RMSE, index_w10 + 1 + diffusion_offset[diffusion_type], marker_type[diffusion_type], marker_color[rise]
                point_list.append(plot_tuple)

    # Creating the axis
    fig = plt.figure(figsize=fig_size)
    ax = fig.add_subplot(111)
    # Y axis = Depth axis
    ax.set_ylabel(y_label, fontsize=ax_label_size)
    ax.set_ylim((0, 6))
    ax.tick_params(axis='both', labelsize=ax_label_size)
    # X axis = Concentration axis
    ax.set_xlabel(x_label, fontsize=ax_label_size)
    ax.set_xlim((0, 1))

    # Now, plotting the points
    for point in point_list:
        RMSE, index_w10, marker, color = point
        ax.plot(RMSE, index_w10, color=color, marker=marker, alpha=0.7, markersize=10)

    # Now, altering the Y axis to list the wind speeds instead of the simple labels 1 - 5
    ax.set_yticks(range(7))
    ax.set_yticklabels(['', 0.85, 2.40, 4.35, 6.65, 9.30, ''])

    # Next, adding a legend to explain the color scheme and the marker type
    label_marker = ['KPP', 'SWB']
    marker = [plt.plot([], [], c='k', markersize=10, marker=marker_type[label], label=label, linestyle='')[0] for label in
              label_marker]

    def label_color(rise):
        return r'$w_r$ = {}'.format(np.abs(rise)) + ' m s$^{-1}$'

    color = [plt.plot([], [], c=marker_color[rise], markersize=10, marker='o', label=label_color(rise), linestyle='')[0] for rise in
              w_r]

    ax.legend(handles=marker + color, fontsize=legend_size, loc='lower right')

    plt.savefig(settings.figure_dir + 'model_evaluation_markov_0.png', bbox_inches='tight', dpi=600)