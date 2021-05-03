import math

import numpy as np
import scipy.optimize
import settings
from utils.utils_files import check_file_exist, save_obj, load_obj, find_nearest_index


def beaufort_limits():
    return [(0, 0.2), (0.2, 1.5), (1.5, 3.3), (3.3, 5.4), (5.4, 7.9), (7.9, 10.7), (10.7, 13.8), (13.8, 17.1),
            (17.1, 20.7), (20.7, 24.4), (24.4, 28.4), (28.5, 32.6), (32.7, 100)]


def determine_Cd(w_10):
    """ Following Large & Pond (1981) https://doi.org/10.1175/1520-0485(1981)011%3C0324:OOMFMI%3E2.0.CO;2 """
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
    depth = np.abs(depth)
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
        profile += settings.bulk_diff
    elif diffusion_type == 'KPP':
        alpha = (k * u_s) / phi
        profile = alpha * (depth + z0) * np.power(1 - depth / mld, 2)
        profile[depth > mld] = 0
        profile += settings.bulk_diff
    elif diffusion_type == 'artificial':
        A = 1e-1
        z_f = 5
        profile = A * np.exp(-depth / z_f) * np.sin(np.pi * depth / mld)
        # profile[(depth > 5) & (depth < 7)] = 0.0002
        profile += math.fabs(profile.min()) + settings.bulk_diff
    return profile


def get_vertical_diffusion_gradient_profile(w_10, depth: np.array, diffusion_type: str, mld: float = settings.MLD,
                                            H_s_frac=1.):
    depth = np.abs(depth)
    rho_w = settings.rho_w  # density sea water (kg/m^3)
    rho_a = settings.rho_a  # density air (kg/m^3)
    u_s = math.sqrt(determine_tau(w_10, rho_a) / rho_w)  # shear velocity
    k = settings.vk  # von Karman constant
    H_s = determine_wave_height(w_10)  # significant wave height (m)
    phi = settings.phi
    z0 = determine_surface_roughness(w_10)
    if diffusion_type == 'Kukulka':
        profile = -2.25 * k * u_s * H_s * (H_s * H_s_frac) ** 1.5 * np.power(depth, -2.5) * np.ones(depth.shape)
        profile[depth < (H_s * H_s_frac)] = 0
    elif diffusion_type == 'KPP':
        alpha = (k * u_s) / (phi * mld ** 2)
        profile = alpha * (mld - depth) * (mld - 3 * depth - 2 * z0)
        profile[depth > mld] = 0
    elif diffusion_type == 'artificial':
        A = 1e-2
        z_f = 5
        profile = -A / z_f * np.exp(-depth / z_f) * (-10*np.pi*np.cos(np.pi * depth / mld) + mld * np.sin(np.pi * depth / mld)) / (z_f * mld)
        # profile[(depth > 5) & (depth < 7)] = 0
    return profile


def get_T_L_profile(T_L_amp: float, depth: np.array, T_L_min: float = 2 * settings.dt_int.seconds):
    profile = T_L_min + T_L_amp * np.sin(np.pi * depth / settings.MLD)
    profile[depth > settings.MLD] = T_L_min
    return profile


def lagrangian_integral_timescale(w_10):
    """
    Determining the Lagrangian integral timescale following Denman & Gargett (1983). If the mixed layer depth is deeper
    than the turbulent Ekman layer thickness, then the Lagrangian integral timescale is dependent only on the coriolis
    parameter. Otherwise, when the mixed layer depth is less than the Ekman layer thickness, then the timescale becomes
    inversely proportional to the surface wind speed

    The ekman timescale is approximately 30 minutes at latitude 45 (MLD > L_E)
    """
    # Friction velocity
    u_w = math.sqrt(determine_tau(w_10, settings.rho_a) / settings.rho_w)
    # RMS eddy velocity
    u_t = 2 * u_w
    # Coriolis parameter
    f = 2 * 7.2921e-5 * math.sin(settings.latitude)
    # The turbulent Ekman layer thickness
    L_e = 0.4 * u_w / f
    # Now determining the time scale
    if L_e < settings.MLD:
        T_L = 0.2 / f
    else:
        T_L = settings.MLD / u_t
    print("The lagrangian integral timescale is {} minutes".format(T_L / 60.))
    return T_L


def determine_particle_size(w_rise):
    """
    Determining the equivalent spherical particle size for a rise velocity based on Enders et al. (2015
    https://doi.org/10.1016/j.marpolbul.2015.09.027
    We do this by means of a pre-calculated lookup table for PE and PP particles with varying sizes of L, as for some
    reason I don't get correct results when
    """
    lookup_file = settings.data_dir + 'particle_size_rise_velocity_lookup'
    if not check_file_exist(lookup_file):
        # The optimization function
        def to_optimize(w_rise):
            if material is 'PE':
                rho_p = settings.rho_p_pe
            elif material is 'PP':
                rho_p = settings.rho_p_pp
            left = (1. - rho_p / settings.rho_w) * 8. / 3. * L * settings.g
            Re = 2. * L * np.abs(w_rise) / settings.nu
            right = np.square(w_rise) * (24. / Re + 5. / np.sqrt(Re) + 2. / 5.)
            return np.abs(left - right)
        # The range of particle sizes for which we will compute the rise velocities for
        L_range = np.logspace(0, -5, 2000)
        lookup_dict = {'L': L_range, 'w_rise_PP': np.zeros(L_range.shape, dtype=float),
                       'w_rise_PE': np.zeros(L_range.shape, dtype=float)}
        for ind, L in enumerate(L_range):
            material = 'PP'
            lookup_dict['w_rise_PP'][ind] = scipy.optimize.minimize_scalar(to_optimize, bounds=[-100, 0], method='bounded').x
            material = 'PE'
            lookup_dict['w_rise_PE'][ind] = scipy.optimize.minimize_scalar(to_optimize, bounds=[-100, 0], method='bounded').x

        save_obj(lookup_file, lookup_dict)
    # Loading the lookup table
    lookup_dict = load_obj(lookup_file)
    ind_PP = find_nearest_index(lookup_dict['w_rise_PP'], w_rise)
    ind_PE = find_nearest_index(lookup_dict['w_rise_PE'], w_rise)
    L_range = lookup_dict['L']
    str_format = w_rise, 2 * L_range[ind_PP], 2 * L_range[ind_PE]
    output = 'A rise velocity of {:.2E} is approximately a PP particle of diameter {:.2E} or a PE ' \
             'particle of diameter {:.2E}'.format(*str_format)
    print(output)


def determine_mixed_layer(w_10, w_rise, diffusion_type='KPP'):
    """
    Determining the depth of the surface layer within particles are assumed to be distributed homogeneously, following
    the boundary condition approach of Ross & Sharples (2004)
    """

    def to_optimize(z_t):
        dK = get_vertical_diffusion_gradient_profile(w_10, np.array([z_t]), diffusion_type)
        dt = settings.dt_int.seconds
        K = get_vertical_diffusion_profile(w_10, np.array([z_t]), diffusion_type)
        RHS = dK[0] * dt + np.sqrt(6 * K[0] * dt) - w_rise * dt
        return np.abs(z_t - RHS)

    mixing_depth = scipy.optimize.minimize_scalar(to_optimize, bounds=[0, 100], method='bounded').x
    print('The surface turbulent mixed layer depth {} m'.format(mixing_depth))
    return mixing_depth


def DeleteParticle(particle, fieldset, time):
    particle.delete()