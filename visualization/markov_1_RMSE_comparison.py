import settings
import matplotlib.pyplot as plt
import utils
import numpy as np
import analysis
import visualization.utils_visualization


def markov_1_RMSE_comparison(x_label=r'$u_{10}$ (m s$^{-1}$)', y_label=r'RMSE', fig_size=(16, 8),
                             ax_label_size=16, legend_size=10):
    w_10 = [0.85, 2.4, 4.35, 6.65, 9.3]
    w_r = [-0.03, -0.003]
    alpha = [0.0, 0.1, 0.3, 0.5, 0.7, 0.95]
    rise_offset = {-0.03: 0.15, -0.003: -0.15}
    # Selecting the marker type according which form of diffusion it is:
    marker_type = {-0.03: 'o', -0.003: 'X'}

    # Looping through the KPP simulations, and retrieving the RMSE values for them
    point_list_KPP = []
    for index_w10, wind in enumerate(w_10):
        for rise in w_r:
            RMSE = analysis.determine_RMSE(wind, rise, 'KPP', 'Ceiling', alpha=0.0, output=True)
            plot_tuple = RMSE, index_w10 + 1, marker_type[rise], visualization.utils_visualization.return_color(0)
            point_list_KPP.append(plot_tuple)

    for index_w10, wind in enumerate(w_10):
        for rise in w_r:
            for index_a, a in enumerate(alpha):
                RMSE = analysis.determine_RMSE(wind, rise, 'KPP', 'Ceiling_Markov', alpha=a, output=True)
                plot_tuple = RMSE, index_w10 + 1 + rise_offset[rise], marker_type[rise], visualization.utils_visualization.return_color(index_a + 1)
                point_list_KPP.append(plot_tuple)

    # Looping through the SWB simulations, and retrieving the RMSE values for them
    point_list_Kukulka = []
    for index_w10, wind in enumerate(w_10):
        for rise in w_r:
            RMSE = analysis.determine_RMSE(wind, rise, 'Kukulka', 'Ceiling', alpha=0.0, output=True)
            plot_tuple = RMSE, index_w10 + 1, marker_type[rise], visualization.utils_visualization.return_color(0)
            point_list_Kukulka.append(plot_tuple)

    for index_w10, wind in enumerate(w_10):
        for rise in w_r:
            for index_a, a in enumerate(alpha):
                RMSE = analysis.determine_RMSE(wind, rise, 'Kukulka', 'Ceiling_Markov', alpha=a, output=True)
                plot_tuple = RMSE, index_w10 + 1 + rise_offset[rise], marker_type[rise], visualization.utils_visualization.return_color(index_a + 1)
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
    # Y axis = Depth axis
    ax2.set_xlim((0, 6))
    ax2.set_xlabel(x_label, fontsize=ax_label_size)
    ax2.tick_params(axis='both', labelsize=ax_label_size)
    ax2.tick_params(labelleft=False)
    # X axis = Concentration axis
    # ax2.set_ylabel(y_label, fontsize=ax_label_size)
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


    # # Next, adding a legend to explain the color scheme and the marker type
    def label_marker(rise):
        return r'$w_r$ = {}'.format(np.abs(rise)) + ' m s$^{-1}$'

    marker = [plt.plot([], [], c='k', markersize=10, marker=marker_type[rise], label=label_marker(rise), linestyle='')[0] for rise in
              w_r]
    markov0 = [plt.plot([], [], c=visualization.utils_visualization.return_color(0), markersize=10, marker='o', label='M0', linestyle='')[0]]
    markov1 = [plt.plot([], [], c=visualization.utils_visualization.return_color(ind + 1), markersize=10, marker='o', label=r'M1 - $\alpha = $' + '{}'.format(a), linestyle='')[0] for ind, a in
               enumerate(alpha)]

    ax2.legend(handles=marker + markov0 + markov1, fontsize=legend_size, loc='upper right')

    plt.savefig(settings.figure_dir + 'model_evaluation_markov_1.png', bbox_inches='tight', dpi=600)