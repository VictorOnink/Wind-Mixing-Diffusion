# A number of analysis functions for determining the e-folding scale of the vertical distribution of particles
import settings as SET
from parcels_simulation_functions import get_parcels_output_name
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt


def depth_concentration(k_z, w_10, w_rise):
    dataset = Dataset(get_parcels_output_name(k_z, w_10, w_rise))
    depth = dataset.variables['z'][:,-1]
    # Bins, at 0.5 meter intervals
    bin = np.arange(0, 100, 0.1)
    # The concentration
    concentrations, bin_edges = np.histogram(depth, bins=bin)
    quick_plot(bin_edges,concentrations)


def quick_plot(x, y):
    plt.plot(x[:-1], y, '.', c='r')
    plt.xlim((x.min(), x.max()))
    plt.ylim(y.min(), y.max())
    plt.show()
