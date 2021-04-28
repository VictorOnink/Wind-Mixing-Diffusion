import settings
import matplotlib.pyplot as plt
import utils
from visualization import utils_visualization as utils_v
import numpy as np


def integration_timestep_control(w_rise_select, selection='w_10', close_up=None, y_label='Depth (m)', alpha=0.3,
                                x_label=r'Normalised Concentrations ($C/C_{max}$)', fig_size=(11, 6),
                                ax_label_size=16, legend_size=10, single_select=0,
                                output_step=-1, boundary='Reflect'):
    ax_range = utils_v.get_axes_range(close_up=close_up, norm_depth=False)

    # Selecting which model data we want to plot based on the diffusion type
    kukulka, kpp, artificial = utils_v.boolean_diff_type('all')
    # Selecting the timestep for which we want to plot the data
    # dt_list = [30, 30, 30, 15, 15, 15, 1, 1, 1]
    # line_style = ['-', '-', '-', '--', '--', '--', 'dotted', 'dotted', 'dotted']
    # w_rise = [-0.03, -0.003, -0.0003, -0.03, -0.003, -0.0003, -0.03, -0.003, -0.0003]
    dt_list = [30, 15, 10, 5, 1]
    line_style = ['-'] * len(dt_list)
    w_rise = [w_rise_select] * len(dt_list)

    # Creating the axis, one for KPP and one for Kukulka
    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size, shape=(1, 2), plot_num=2,
                             all_x_labels=True)

    mean_wind = np.mean(utils.utils_physics.beaufort_limits()[4])
    # for axis in ax:
    #     data_line, data_label = utils_v.add_observations(axis, norm_depth=False, alpha=alpha,
    #                                                      wind_range=utils.beaufort_limits()[4])

    # First, plotting the KPP data:
    for count, dt in enumerate(dt_list):
        if kpp:
            profile_dict = utils_v.get_concentration_list([mean_wind], [w_rise[count]], selection, single_select,
                                                          output_step=output_step, diffusion_type='KPP',
                                                          boundary=boundary, alpha_list=[0], dt=dt)
            for counter in range(len(profile_dict['concentration_list'])):
                ax[0].plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                           label=label(w_rise[count], dt),
                           linestyle=line_style[count], color=utils_v.discrete_color_from_cmap(count, len(dt_list)))
        if kukulka:
            profile_dict = utils_v.get_concentration_list([mean_wind], [w_rise[count]], selection, single_select,
                                                          output_step=output_step, diffusion_type='Kukulka',
                                                          boundary=boundary, alpha_list=[0], dt=dt)
            for counter in range(len(profile_dict['concentration_list'])):
                ax[1].plot(profile_dict['concentration_list'][counter], profile_dict['depth_bins'],
                           label=label(w_rise[count], dt),
                           linestyle=line_style[count], color=utils_v.discrete_color_from_cmap(count, len(dt_list)))

    lines, labels = ax[1].get_legend_handles_labels()
    # Adding the legends
    # ax[0].legend(data_line, data_label, fontsize=legend_size, loc='lower right')
    ax[1].legend(lines, labels, fontsize=legend_size, loc='lower right')

    # Adding subplot titles
    ax[0].set_title(r'(a) KPP', fontsize=ax_label_size)
    ax[1].set_title(r'(b) SWB', fontsize=ax_label_size)

    # Saving the figure
    plt.savefig(settings.figure_dir + 'm0_dt_int_check_w_rise={}_mld={}'.format(w_rise_select, settings.MLD) + '.png',
                bbox_inches='tight', dpi=600)


def label(w_rise, dt):
    return 'M-0, $w_r$ = ' + r'{}'.format(np.abs(w_rise)) + r' m s$^{-1}$, dt = ' + '{} s'.format(dt)
