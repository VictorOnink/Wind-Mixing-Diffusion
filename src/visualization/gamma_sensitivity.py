import matplotlib.pyplot as plt
import utils
import settings
from visualization import utils_visualization as utils_v
import numpy as np


class gamma_sensitivity:
    def __init__(self, w_rise, with_observations=True, close_up=None):
        # Setting the rise velocity that we are plotting
        self.w_rise = w_rise
        self.boundary = 'Ceiling'
        self.w_10_list = [0.85, 2.4, 4.35, 6.65, 9.3]
        self.gamma_list = [0.5, 1.0, 1.5, 2.0]
        # Data loading parameters
        self.output_step = -1
        self.single_select = 1
        # Setting plot parameters
        self.with_observations = with_observations
        self.add_variability = True
        self.y_label = 'Depth (m)'
        self.x_label = r'Normalised Concentrations'
        self.close_up = close_up
        self.fig_size = (16, 8)
        self.ax_label_size = 16
        self.legend_size = 12
        self.alpha = 0.15

    def plot(self):
        # Setting the axes ranges
        ax_range = utils_v.get_axes_range(close_up=self.close_up, norm_depth=False)

        # Titles for the subplots
        sub_titles = [r'(a) u$_{10}$=0.2-1.5 m s$^{-1}$', r'(b) u$_{10}$=1.5-3.3 m s$^{-1}$',
                      r'(c) u$_{10}$=3.3-5.4 m s$^{-1}$', r'(d) u$_{10}$=5.4-7.9 m s$^{-1}$',
                      r'(e) u$_{10}$=7.9-10.7 m s$^{-1}$']

        # Get the base figure axis
        plot_num = 6
        ax = utils_v.base_figure(self.fig_size, ax_range, self.y_label, self.x_label, self.ax_label_size,
                                 shape=(2, 3), plot_num=plot_num, legend_axis=False)

        # set subtitles
        for scale in range(sub_titles.__len__()):
            ax[scale].set_title(sub_titles[scale], fontsize=self.ax_label_size)

        # Load wind conditions
        beaufort = utils.beaufort_limits()

        for scale in range(self.w_10_list.__len__()):
            if self.with_observations:
                _, _ = utils_v.add_observations(ax[scale], norm_depth=False, wind_range=beaufort[scale + 1],
                                                alpha=self.alpha, mean_concentrations=True)
            for index_gamma, gamma in enumerate(self.gamma_list):
                plot_color = utils_v.discrete_color_from_cmap(index=index_gamma, subdivisions=self.gamma_list.__len__())
                profile_dict = utils_v.get_concentration_list(self.w_10_list, [self.w_rise], 'all', self.single_select,
                                                              output_step=self.output_step, diffusion_type='SWB',
                                                              boundary=self.boundary, alpha_list=[0], gamma=gamma)
                _, w_rise = profile_dict['parameter_concentrations'][scale]
                concentration = profile_dict['concentration_list'][scale]
                depth = profile_dict['depth_bins']
                ax[scale].plot(concentration, depth, label=self.label_gamma(gamma),
                               linestyle='-', color=plot_color)
                if self.add_variability:
                    std = profile_dict['std_list'][scale]
                    upper_limit, lower_limit = concentration + std, concentration - std
                    ax[scale].fill_betweenx(depth, lower_limit, upper_limit, alpha=0.2, color=plot_color)

            lines, labels = ax[scale].get_legend_handles_labels()

        # Labels for field data
        if self.with_observations:
            field_lines, field_labels = lines[-6:], labels[-6:]

        # Labels for gamma
        gamma_lines = [plt.plot([], [],
                                c=utils_v.discrete_color_from_cmap(index=index, subdivisions=self.gamma_list.__len__()),
                                label=r'$\gamma=$' + '{}'.format(gamma), linestyle='-')[0]
                       for index, gamma in enumerate(self.gamma_list)]

        # Adding the legend
        if self.with_observations:
            handles = field_lines + gamma_lines
        else:
            handles = gamma_lines
        ax[-1].legend(handles=handles, fontsize=self.legend_size, loc='upper left')
        ax[-1].axis('off')

        # Saving the figure
        save_name = settings.figure_dir + 'SWB_gamma_influence_w_r={}.png'.format(w_rise)
        plt.savefig(save_name, bbox_inches='tight')

    @staticmethod
    def label_gamma(gamma):
        label = r'$\gamma=$' + '{}'.format(gamma)
        return label
