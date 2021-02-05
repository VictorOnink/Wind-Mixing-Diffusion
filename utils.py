import settings
import pickle
import math
import os
import numpy as np


def get_parcels_output_name(w_10, w_rise, diffusion_type, boundary, mld=settings.MLD):
    return settings.output_dir + diffusion_type + '_' + boundary + '_w10_{}_w_rise_{}_MLD_{}.nc'.format(w_10, w_rise,
                                                                                                        mld)


def get_concentration_output_name(w_10, w_rise, diffusion_type, boundary, mld=settings.MLD):
    return settings.conc_dir + diffusion_type + '_' + boundary + '_conc_w10_{}_w_rise_{}_MLD_{}'.format(w_10, w_rise,
                                                                                                        mld)


def get_data_output_name(prefix: str):
    return settings.data_dir + 'standardized_data_' + prefix


def save_obj(filename, object):
    with open(filename + '.pkl', 'wb') as f:
        pickle.dump(object, f, pickle.HIGHEST_PROTOCOL)


def load_obj(filename):
    with open(filename + '.pkl', 'rb') as f:
        return pickle.load(f)


def _check_file_exist(File: str):
    return os.path.isfile(File)


def remove_file(conduct: bool, File: str):
    if conduct:
        if _check_file_exist(File):
            os.remove(File)


def find_nearest_index(depth, z_ref):
    return (np.abs(depth - z_ref)).argmin()


def return_color(index):
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray',
              'tab:olive', 'tab:cyan']
    num = len(colors)
    return colors[index % num]


def beaufort_limits():
    return [(0, 0.2), (0.2, 1.5), (1.5, 3.3), (3.3, 5.4), (5.4, 7.9), (7.9, 10.7), (10.7, 13.8), (13.8, 17.1),
            (17.1, 20.7), (20.7, 24.4), (24.4, 28.4), (28.5, 32.6), (32.7, 100)]


def determine_Cd(w_10):
    return min(max(1.2E-3, 1.0E-3 * (0.49 + 0.065 * w_10)), 2.12E-3)


def determine_tau(w_10, rho_air):
    C_d = determine_Cd(w_10)
    return C_d * rho_air * w_10 ** 2


def determine_wave_height(w_10):
    tau_wind = determine_tau(w_10, settings.rho_a)  # wind stress at the ocean surface
    u_air = math.sqrt(tau_wind / settings.rho_a)  # friction velocity of air
    H_s = 0.96 * settings.g ** (-1) * settings.wave_age ** 1.5 * u_air ** 2  # significant wave height (m)
    return H_s


def determine_surface_roughness(w_10):
    # the surface roughness z_0, following Zhao & Li (2019)
    g, rho_a = settings.g, settings.rho_a
    u_s = math.sqrt(determine_tau(w_10, rho_a) / rho_a)  # frictional velocity
    # Rewriting the wave age in terms of the w_10 wind speed, instead of the frictional velocity used by
    # Kukulka et al. (2012)
    beta = settings.wave_age * u_s / w_10
    # Since beta > 0.4 (corresponding to relative mature wind waves), we have from eq. 36 that:
    z_0 = 3.5153e-5 * beta ** (-0.42) * w_10 **2 / g
    return z_0


def determine_kukulka_e_length(w_10, w_rise):
    # Getting environmental parameters
    tau_wind = determine_tau(w_10, settings.rho_a)  # wind stress at the ocean surface
    u_water = math.sqrt(tau_wind / settings.rho_w)  # friction velocity of water
    H_s = determine_wave_height(w_10)  # significant wave height (m)
    A0 = 1.5 * u_water * settings.vk * H_s  # A0 constant Kukulka et al. (2012)
    return math.fabs(w_rise) / A0


def get_vertical_diffusion_profile(w_10, depth: np.array, diffusion_type: str, mld: float = settings.MLD, H_s_frac=1.):
    rho_w = settings.rho_w  # density sea water (kg/m^3)
    rho_a = settings.rho_a  # density air (kg/m^3)
    u_s = math.sqrt(determine_tau(w_10, rho_a) / rho_w)  # shear velocity
    H_s = determine_wave_height(w_10)  # significant wave height (m)
    k = settings.vk  # von Karman constant
    phi = settings.phi
    z0 = determine_surface_roughness(w_10)

    if diffusion_type == 'Kukulka':
        profile = 1.5 * k * u_s * H_s * np.ones(depth.shape)
        profile[depth > (H_s * H_s_frac)] *= ((H_s * H_s_frac) ** 1.5 * np.power(depth[depth > (H_s * H_s_frac)], -1.5))
    elif diffusion_type == 'KPP':
        alpha = (k * u_s) / phi
        profile = alpha * (depth + z0) * np.power(1 - depth / mld, 2)
        profile[depth > mld] = 0
    return profile + settings.bulk_diff


def get_vertical_diffusion_gradient_profile(w_10, depth: np.array, diffusion_type: str, mld: float = settings.MLD,
                                            H_s_frac=1.):
    rho_w = settings.rho_w  # density sea water (kg/m^3)
    rho_a = settings.rho_a  # density air (kg/m^3)
    u_s = math.sqrt(determine_tau(w_10, rho_a) / rho_w)  # shear velocity
    k = settings.vk  # von Karman constant
    H_s = determine_wave_height(w_10)  # significant wave height (m)
    phi = settings.phi
    z0 = determine_wave_height(w_10)
    if diffusion_type == 'Kukulka':
        profile = -2.25 * k * u_s * H_s * (H_s * H_s_frac) ** 1.5 * np.power(depth, -2.5) * np.ones(depth.shape)
        profile[depth < (H_s * H_s_frac)] = 0
    elif diffusion_type == 'KPP':
        alpha = (k * u_s) / (phi * mld ** 2)
        profile = alpha * (mld - depth) * (mld - 3 * depth - 2 * z0)
        profile[depth > mld] = 0
    return profile


def exclude_field_data(exclude, sources):
    if exclude is not None:
        if exclude is not None:
            for source in exclude:
                sources.remove(source)
    return sources