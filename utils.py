import settings
import pickle
import math
import os
import numpy as np


def get_parcels_output_name(k_z, w_10, w_rise, diffusion_type, boundary):
    if diffusion_type == 'Rouse':
        return settings.output_dir + diffusion_type + '_' + boundary + '_k_z_{}_w10_{}_w_rise_{}.nc'.format(k_z, w_10,
                                                                                                            w_rise)
    else:
        return settings.output_dir + diffusion_type + '_' + boundary + '_w10_{}_w_rise_{}.nc'.format(w_10, w_rise)


def get_concentration_output_name(k_z, w_10, w_rise, diffusion_type, boundary):
    if diffusion_type == 'Rouse':
        return settings.conc_dir + diffusion_type + '_' + boundary + '_conc_k_z_{}_w10_{}_w_rise_{}'.format(k_z, w_10,
                                                                                                            w_rise)
    else:
        return settings.conc_dir + diffusion_type + '_' + boundary + '_conc_w10_{}_w_rise_{}'.format(w_10, w_rise)


def save_obj(filename, object):
    with open(filename + '.pkl', 'wb') as f:
        pickle.dump(object, f, pickle.HIGHEST_PROTOCOL)


def load_obj(filename):
    with open(filename + '.pkl', 'rb') as f:
        return pickle.load(f)


def determine_Cd(w_10):
    return min(max(1.2E-3, 1.0E-3*(0.49+0.065*w_10)), 2.12E-3)


def determine_tau(w_10, rho_air):
    C_d = determine_Cd(w_10)
    return C_d * rho_air * w_10**2


def determine_kukulka_e_length(w_10, w_rise):
    # basic constants
    rho_w = 1027                                                # density sea water (kg/m^3)
    rho_a = 1.22                                                # density air (kg/m^3)
    vk = 0.4                                                    # von Karman constant
    wave_age = 35                                               # assuming fully developed sea, Kukulka et al. (2012)
    g = 9.81                                                    # acceleration due to gravity (m s^-2)
    # Getting environmental parameters
    tau_wind = determine_tau(w_10, rho_a)                 # wind stress at the ocean surface
    u_water = math.sqrt(tau_wind / rho_w)                       # friction velocity of water
    u_air = math.sqrt(tau_wind / rho_a)                         # friction velocity of air
    H_s = 0.96 * g ** (-1) * wave_age ** 1.5 * u_air ** 2       # significant wave height (m)
    A0 = 1.5 * u_water * vk * H_s                               # A0 constant Kukulka et al. (2012)
    return math.fabs(w_rise) / A0


def _check_file_exist(File: str):
    return os.path.isfile(File)


def remove_file(conduct: bool, File:str):
    if conduct:
        if _check_file_exist(File):
            os.remove(File)


def get_vertical_diffusion_profile(w_10, k_z, depth: np.array, diffusion_type: str):
    rho_w = settings.rho_w # density sea water (kg/m^3)
    rho_a = settings.rho_a # density air (kg/m^3)
    wave_age = settings.wave_age # wave age of developed wave field (Kukulka et al., 2012)
    g = settings.g # Gravitational acceleration (
    u_s = math.sqrt(determine_tau(w_10, rho_a) / rho_w) # shear velocity
    u_air = math.sqrt(determine_tau(w_10, rho_a) / rho_a)
    H_s = 0.96 * g ** (-1) * wave_age ** 1.5 * u_air ** 2  # significant wave height (m)
    k = settings.vk # von Karman constant
    phi = settings.phi
    z0 = settings.z0
    MLD = settings.MLD
    if diffusion_type == 'Rouse':
        return k * u_s * depth + k_z
    if diffusion_type == 'Kukulka':
        return 1.5 * k * u_s * H_s * np.ones(depth.shape)
    if diffusion_type == 'KPP':
        alpha = (k * u_s) / phi
        return alpha * (depth + z0) * np.power(1 - depth / MLD, 2)


def get_vertical_diffusion_gradient_profile(w_10, depth: np.array, diffusion_type: str):
    rho_w = settings.rho_w # density sea water (kg/m^3)
    rho_a = settings.rho_a # density air (kg/m^3)
    u_s = math.sqrt(determine_tau(w_10, rho_a) / rho_w) # shear velocity
    k = settings.vk # von Karman constant
    phi = settings.phi
    z0 = settings.z0
    MLD = settings.MLD
    if diffusion_type == 'Rouse':
        return k * u_s * np.ones(depth.shape)
    if diffusion_type == 'Kukulka':
        return np.zeros(depth.shape)
    if diffusion_type == 'KPP':
        alpha = (k * u_s) / (phi * MLD**2)
        return alpha * (MLD - depth) * (MLD - 3*depth - 2 * z0)