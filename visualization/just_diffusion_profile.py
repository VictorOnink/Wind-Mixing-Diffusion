import settings
import matplotlib.pyplot as plt
import utils
import utils.utils_physics
from visualization import utils_visualization as utils_v
import numpy as np


def just_diffusion_profile(w_10_list, y_label='Depth (m)', x_label=r'$K_z$ ($10^{-2}$ m$^2$ s$^{-1}$)',
                           fig_size=(8, 8), ax_label_size=16, legend_size=12):
    ymax, ymin = 0, -1 * (settings.MLD + 10)
    ax_range = (1.55, 0, ymax, ymin)
    depth = np.linspace(ymax, np.abs(ymin), 1000)

    # Creating the axis
    ax = utils_v.base_figure(fig_size, ax_range, y_label, x_label, ax_label_size)

    # Plotting the diffusion profile according to the Kukulka approach
    for count, w_10 in enumerate(w_10_list):
        profile = utils.utils_physics.get_vertical_diffusion_profile(w_10, depth, 'Kukulka')
        ax.plot(profile * 100, -1 * depth, color=utils_v.discrete_color_from_cmap(count, len(w_10_list)), linestyle='-',
                label=utils_v.label_diffusivity_profile(w_10, 'Kukulka'))

    # Plotting the diffusion profile according the KPP approach
    for count, w_10 in enumerate(w_10_list):
        profile = utils.utils_physics.get_vertical_diffusion_profile(w_10, depth, 'KPP')
        ax.plot(profile * 100, -1 * depth, color=utils_v.discrete_color_from_cmap(count, len(w_10_list)),
                linestyle='--', label=utils_v.label_diffusivity_profile(w_10, 'KPP'))

    # Adding the legend
    ax.legend(fontsize=legend_size, loc='lower right')

    plt.savefig(settings.figure_dir + 'diffusion_profile_MLD=20.png', bbox_inches='tight')