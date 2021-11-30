import matplotlib.pyplot as plt
import settings
import utils
from visualization import utils_visualization as utils_v
import string
import analysis


def plot_error_correlaton(diffusion, x_label=r'Correlation ($\rho$)', y_variable='RMSE', markov=False,
                          boundary='Ceiling', legend_size=12, label_size=14, title_size=12):
    """
    plot
    :return:
    """
    # Setting the y axis label
    if y_variable == 'RMSE':
        y_label = 'RMSE'

    # Setting some basic lists of parameters
    w_10_list = [0.85, 2.4, 4.35, 6.65, 9.3]
    w_10_marker = ['o', 'v', 's', 'P', 'X']
    alpha_list = [0.0, 0.1, 0.3, 0.5, 0.7, 0.95]
    if diffusion == 'KPP':
        theta_list = [1.0, 2.0, 3.0, 4.0, 5.0]
    else:
        theta_list = [1.0]
    wave_roughness_list = [False, False, True, True]
    w_rise = [-0.03, -0.003]

    # Creating the figure
    ax_range = 1.1, 0.5, 0.15, 0
    plot_num = 4
    ax = utils_v.base_figure(fig_size=(12, 10), ax_range=ax_range, y_label=y_label, x_label=x_label, ax_label_size=label_size,
                             all_x_labels=True, legend_axis=True, shape=(2, 2), plot_num=plot_num, log_xscale=False,
                             width_ratios=[1, 1, 0.4])

    # Setting the subplots
    for index in range(ax.__len__() - 2):
        ax[index].set_title(subplot_title(index), fontsize=title_size)

    # Plotting the correlation and RMSE
    for index_ax in range(plot_num):
        for index_w_10, w_10 in enumerate(w_10_list):
            for index_theta, theta in enumerate(theta_list):
                RMSE = analysis.determine_RMSE(w_10=w_10, w_rise=w_rise[index_ax % 2], diffusion_type=diffusion,
                                               alpha=0.0, boundary=boundary, theta=theta, output=True,
                                               wave_roughness=wave_roughness_list[index_ax])
                r, _ = analysis.correlation_field_model_data(w_10=w_10, w_rise=w_rise[index_ax % 2],
                                                             diffusion_type=diffusion, boundary=boundary, alpha=0.0,
                                                             theta=theta, to_print=False,
                                                             wave_roughness=wave_roughness_list[index_ax])
                color = utils_v.discrete_color_from_cmap(index=index_theta, subdivisions=theta_list.__len__())
                ax[index_ax].plot(r, RMSE, marker=w_10_marker[index_w_10], color=color, linestyle=None)

    # Creating the legend
    theta_lines = [plt.plot([], [], 'o', c=utils_v.discrete_color_from_cmap(index=index, subdivisions=theta_list.__len__()),
                            label=r'$\theta=$' + '{}'.format(theta))[0]
                   for index, theta in enumerate(theta_list)]
    wind_line = [plt.plot([], [], w_10_marker[ind], c='k', label='{:.2f}'.format(w_10) + r' m s$^{-1}$')[0] for ind, w_10 in enumerate(w_10_list)]
    ax[-2].legend(handles=wind_line + theta_lines, fontsize=legend_size, loc='upper left')
    ax[-2].axis('off')

    # Saving the figure
    file_name = settings.figure_dir + 'error_correlation_{}_markov={}.png'.format(diffusion, markov)
    plt.savefig(file_name, bbox_inches='tight', dpi=600)


def subplot_title(index):
    abc = string.ascii_lowercase
    w_rise = {0: 0.03, 1: 0.003}
    if index < 2:
        return '({}) {} m s'.format(abc[index], w_rise[index % 2]) + r'$^{-1}$, $z_0=$ Zhao & Li (2019)'
    else:
        return '({}) {} m s'.format(abc[index], w_rise[index % 2]) + r'$^{-1}$, $z_0=0.1\times H_s$'

