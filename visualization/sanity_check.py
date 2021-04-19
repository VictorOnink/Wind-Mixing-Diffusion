import settings
import matplotlib.pyplot as plt
import utils
import numpy as np
from netCDF4 import Dataset

import utils.utils_filenames


def sanity_check(w_10, w_rise, diffusion_type, boundary, alpha_list):
    dataset = Dataset(
        utils.utils_filenames.get_parcels_output_name(w_10, w_rise, diffusion_type, boundary=boundary, mld=settings.MLD,
                                                      alpha=alpha_list[0]))
    start, end = 0, -1
    depth = dataset.variables['z'][0, start:end]
    ratio = dataset.variables['ratio'][0, start:end]
    sig2 = dataset.variables['sig2_store'][0, start:end]
    w_T = dataset.variables['w_T_store'][0, start:end]
    w_T2 = np.square(w_T)
    dWz = dataset.variables['dWz_store'][0, start:end]
    dHis = dataset.variables['d_History_store'][0, start:end]
    dGrad = dataset.variables['d_Gradient_store'][0, start:end]
    dDiff = dataset.variables['d_diff_store'][0, start:end]

    dim = 8
    fig = plt.figure(figsize=(dim, dim))

    var = [depth, sig2, w_T, w_T2, ratio, dHis, dGrad, dDiff]
    title = ['depth', 'sig2', 'w_T', 'w_T2', 'ratio', 'dHis', 'dGrad', 'dDiff']
    shape = (len(var), 1)
    count = 1

    for row in range(shape[0]):
        for column in range(shape[1]):
            ax_sub = fig.add_subplot(shape[0], shape[1], count)
            ax_sub.plot(var[count - 1])
            ax_sub.set_ylabel(title[count - 1])
            count += 1
    # ax_sub.set_ylim([0, 600])

    plt.tight_layout()
    plt.savefig(settings.figure_dir + '/sanity')