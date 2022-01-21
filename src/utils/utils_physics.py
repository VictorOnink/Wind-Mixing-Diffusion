import math
import numpy as np
import scipy.optimize
import settings
from utils.utils_files import check_file_exist, save_obj, load_obj, find_nearest_index


def beaufort_limits():
    """
    Min and Max wind speed values corresponding to the various Beaufort sea states according to
    https://nl.wikipedia.org/wiki/Schaal_van_Beaufort
    :return:
    """
    return [(0, 0.2), (0.2, 1.5), (1.5, 3.3), (3.3, 5.4), (5.4, 7.9), (7.9, 10.7), (10.7, 13.8), (13.8, 17.1),
            (17.1, 20.7), (20.7, 24.4), (24.4, 28.4), (28.5, 32.6), (32.7, 100)]


def determine_Cd(w_10):
    """
    Determining the wind speed drag coefficient
    Following Large & Pond (1981) https://doi.org/10.1175/1520-0485(1981)011%3C0324:OOMFMI%3E2.0.CO;2
    """
    return min(max(1.2E-3, 1.0E-3 * (0.49 + 0.065 * w_10)), 2.12E-3)


def determine_tau(w_10, rho_air):
    """
    Determing the surface wind stress according to the 10m wind speed and the air density
    :param w_10:
    :param rho_air:
    :return:
    """
    # Determining the drag coefficient
    C_d = determine_Cd(w_10)
    return C_d * rho_air * w_10 ** 2


def determine_wave_height(w_10):
    """
    Given the 10m wind speed, determining the significant wave height following Kukulka et al. (2012)
    :param w_10:
    :return:
    """
    tau_wind = determine_tau(w_10, settings.rho_a)  # wind stress at the ocean surface
    u_air = math.sqrt(tau_wind / settings.rho_a)  # friction velocity of air
    H_s = 0.96 * settings.g ** (-1) * settings.wave_age ** 1.5 * u_air ** 2  # significant wave height (m)
    return H_s


def determine_surface_roughness(w_10):
    """
    the surface roughness z_0 for a given 10m wind speed, following Zhao & Li (2019)
    :param w_10:
    :return:
    """
    # determining the frictional velocity following Kukulka et al. (2012)
    u_s = math.sqrt(determine_tau(w_10, settings.rho_a) / settings.rho_a)
    # Rewriting the wave age in terms of the w_10 wind speed, instead of the frictional velocity used by
    # Kukulka et al. (2012)
    beta = settings.wave_age * u_s / w_10
    # Since beta > 0.4 (corresponding to relative mature wind waves), we have from eq. 36 that:
    z_0 = 3.5153e-5 * beta ** (-0.42) * w_10 ** 2 / settings.g
    return z_0


def get_vertical_diffusion_profile(w_10, depth: np.array, diffusion_type: str, mld: float = settings.MLD, gamma=1.,
                                   theta=1, wave_roughness=False):
    """
    Determining the vertical diffusion profile at the given depth levels for the given diffusion type
    :param w_10: 10m wind speeds
    :param depth: the depth array (can be negative in the case of the eulerian model setup)
    :param diffusion_type: KPP or SWB or artifial (where that was only a test case)
    :param mld: the mixed layer depth
    :param gamma: Relevent for SWB diffusion, setting until what depth we have constant diffusion as a fraction of
                     the significant wave height
    :param theta: Langmuir circulation amplification factor
    :param wave_roughness: if True, use the 0.1 * significant wave height as the surface roughness
    :return: array with Kz value at each depth level
    """
    depth = np.abs(depth)
    # Determining the surface friction velocity
    u_s = math.sqrt(determine_tau(w_10, settings.rho_a) / settings.rho_w)
    # Getting the significant wave height
    H_s = determine_wave_height(w_10)
    # Getting the surface roughness
    if wave_roughness:
        z0 = 0.1 * H_s
    else:
        z0 = determine_surface_roughness(w_10)
    if diffusion_type == 'SWB':
        profile = 1.5 * settings.vk * u_s * H_s * np.ones(depth.shape)
        profile[depth > (H_s * gamma)] *= ((H_s * gamma) ** 1.5 * np.power(depth[depth > (H_s * gamma)], -1.5))
        profile += settings.bulk_diff
    elif diffusion_type == 'KPP':
        alpha = (settings.vk * u_s * theta) / settings.phi
        profile = alpha * (depth + z0) * np.power(1 - depth / mld, 2)
        profile[depth > mld] = 0
        profile += settings.bulk_diff
    elif diffusion_type == 'artificial':
        A = 1e-1
        z_f = 5
        profile = A * np.exp(-depth / z_f) * np.sin(np.pi * depth / mld)
        profile += math.fabs(profile.min()) + settings.bulk_diff
    return profile


def get_vertical_diffusion_gradient_profile(w_10, depth: np.array, diffusion_type: str, mld: float = settings.MLD,
                                            gamma=1., theta=1, wave_roughness=False):
    """
    Analytically determining vertical gradient of the vertical diffusion profile at the given depth levels
    :param w_10: 10m wind speeds
    :param depth: the depth array (can be negative in the case of the eulerian model setup)
    :param diffusion_type: KPP or SWB or artifial (where that was only a test case)
    :param mld: the mixed layer depth
    :param gamma: Relevent for SWB diffusion, setting until what depth we have constant diffusion as a fraction of
                     the significant wave height
    :param theta: Langmuir circulation amplification factor
    :param wave_roughness: if True, use 0.1 * significant wave height as the surface roughness
    :return: array with dKz value at each depth level
    """
    depth = np.abs(depth)
    # Determining the surface friction velocity
    u_s = math.sqrt(determine_tau(w_10, settings.rho_a) / settings.rho_w)
    # Getting the significant wave height
    H_s = determine_wave_height(w_10)
    # Getting the surface roughness
    if wave_roughness:
        z0 = 0.1 * H_s
    else:
        z0 = determine_surface_roughness(w_10)
    if diffusion_type == 'SWB':
        profile = -2.25 * settings.vk * u_s * H_s * (H_s * gamma) ** 1.5 * np.power(depth, -2.5) * np.ones(depth.shape)
        profile[depth < (H_s * gamma)] = 0
    elif diffusion_type == 'KPP':
        alpha = (settings.vk * u_s * theta) / (settings.phi * mld ** 2)
        profile = alpha * (mld - depth) * (mld - 3 * depth - 2 * z0)
        profile[depth > mld] = 0
    elif diffusion_type == 'artificial':
        A = 1e-2
        z_f = 5
        profile = -A / z_f * np.exp(-depth / z_f) * (-10*np.pi*np.cos(np.pi * depth / mld) + mld * np.sin(np.pi * depth / mld)) / (z_f * mld)
    return profile


def determine_particle_size(w_rise, conduct=False):
    """
    Determining the equivalent spherical particle size for a rise velocity based on Enders et al. (2015)
    https://doi.org/10.1016/j.marpolbul.2015.09.027
    We do this by means of a pre-calculated lookup table for PE and PP particles with varying sizes of L, as for some
    reason I don't get correct results when I go from rise velocity to size
    """
    if conduct:
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
                lookup_dict['w_rise_PP'][ind] = scipy.optimize.minimize_scalar(to_optimize, bounds=[-100, 0],
                                                                               method='bounded').x
                material = 'PE'
                lookup_dict['w_rise_PE'][ind] = scipy.optimize.minimize_scalar(to_optimize, bounds=[-100, 0],
                                                                               method='bounded').x

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
    """
    Kernel for deleting a parcels particle when it goes out of bounds
    """
    particle.delete()