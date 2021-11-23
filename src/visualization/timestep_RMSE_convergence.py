import matplotlib.pyplot as plt
import analysis, utils, settings
from visualization import utils_visualization as utils_v
import numpy as np


def timestep_RMSE_convergence(y_label='RMSE', x_label=r'$\Delta t$ (seconds)', fig_size=(16, 20), ax_label_size=16,
                              legend_size=12, boundary='Reflect'):
    """
    showing RMSE values of the concentration profiles with different integration timesteps relative to the concentration
    profile with dt = 1 second
    :param y_label: label of the y axis
    :param x_label: label of the x axis
    :param fig_size: size of the figure
    :param ax_label_size: fontsize of the axis labels
    :param legend_size: fontsize of the legend
    :param boundary: which boundary condition, and M-0 or M-1
    :return:
    """
    # Setting the axis limits
    xmax, xmin = 31, 0
    ymax, ymin = 0.2, 0
    ax_range = xmax, xmin, ymax, ymin

    # Get the base figure axis
    shape = (5, 2)
    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size, shape=shape, plot_num=10,
                             all_x_labels=True)

    # The w_rise of interest, and the dt of interest
    w_rise_list = [-0.03, -0.003, -0.0003]
    dt_list = [1, 5, 10, 15, 30]

    for row in range(shape[0]):
        # Setting the wind speed
        wind_range = utils.beaufort_limits()[row + 1]
        mean_wind = np.mean(wind_range)

        # Adding plot titles
        title_dict = {0: 'a', 1: 'b', 2: 'c', 3: 'd', 4: 'e', 5: 'f', 6: 'g', 7: 'h', 8: 'i', 9: 'j'}
        ax[2 * row + 0].set_title('({}) SWB, '.format(title_dict[2 * row + 0]) + r'$u_{10}$ = ' +
                                  '{:.2f}'.format(mean_wind) + r' m s$^{-1}$', fontsize=ax_label_size)
        ax[2 * row + 1].set_title('({}) KPP, '.format(title_dict[2 * row + 1]) + r'$u_{10}$ = ' +
                                  '{:.2f}'.format(mean_wind) + r' m s$^{-1}$', fontsize=ax_label_size)

        # Plotting the KPP RMSE trends
        for count, w_rise in enumerate(w_rise_list):
            RMSE = []
            for dt in dt_list:
                RMSE.append(analysis.reference_RMSE_difference(w_rise=w_rise, w_10=mean_wind, dt=dt,
                                                               diffusion_type='KPP'))
            ax[2 * row + 1].plot(dt_list, RMSE, linestyle='-', color=utils_v.return_color(count),
                                 label=r'$w_{rise}$ = ' + '{}'.format(w_rise) + r'm s$^{-1}$')

        # Plotting the SWB RMSE trends
        for count, w_rise in enumerate(w_rise_list):
            RMSE = []
            for dt in dt_list:
                RMSE.append(analysis.reference_RMSE_difference(w_rise=w_rise, w_10=mean_wind, dt=dt,
                                                               diffusion_type='SWB'))
            ax[2 * row + 0].plot(dt_list, RMSE, linestyle='-', color=utils_v.return_color(count),
                                 label=r'$w_{rise}$ = ' + '{}'.format(w_rise) + r'm s$^{-1}$')

    # Adding a legend
    ax[0].legend(fontsize=legend_size, loc='upper left')

    # Saving the figure
    plt.savefig(save_figure_name(boundary=boundary), bbox_inches='tight')


def save_figure_name(boundary):
    """" Setting the filename of the figure """
    boundary_dict = {'Reflect': 'M0', 'Reflect_Markov': 'M1'}
    if boundary_dict[boundary] is 'M0':
        filename = settings.figure_dir + 'timestep_RMSE_convergence_{}.png'.format(boundary_dict[boundary])
    else:
        filename = settings.figure_dir + 'timestep_RMSE_convergence_{}.png'.format(boundary_dict[boundary])
    return filename
