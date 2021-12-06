import settings


def get_parcels_output_name(w_10, w_rise, diffusion_type, boundary, alpha, mld=settings.MLD, dt=settings.dt_int.seconds,
                            theta=1.0, wave_roughness=False):
    """
    Setting the filename for the parcels output
    :param w_10: 10m wind speed
    :param w_rise: rise velocity
    :param diffusion_type: KPP or SWB diffusion
    :param boundary: the boundary type, and whether it is a M-0 or M-1 simulation
    :param alpha: memory term for M-1 simulations
    :param mld: mixed layer depth
    :param dt: integration timestep
    :param theta: Langmuir circulation amplification factor
    :param wave_roughness: if True, have surface roughness be wave height dependent
    :return:
    """
    name = settings.output_dir + '{}_{}_w10_{}_w_rise_{}_MLD_{}'.format(diffusion_type, boundary, w_10, w_rise, mld)
    if 'Markov' in boundary:
        name += '_alpha_list={}'.format(alpha)
    # Relic when I had already run the dt = 1s simulation and didn't want to redo them due to how long it took
    if dt != 1:
        name += '_dt={}'.format(dt)
    if diffusion_type == 'KPP':
        name += '_theta={}'.format(theta)
        if wave_roughness:
            name += '_wave_roughness'
    return name + '.nc'


def get_concentration_output_name(w_10, w_rise, diffusion_type, boundary, alpha=None, mld=settings.MLD,
                                  dt=settings.dt_int.seconds, theta=1.0, with_pkl=False, wave_roughness=False):
    """
    Setting the file name of the concentration output
    :param w_10:
    :param w_rise:
    :param diffusion_type:
    :param boundary:
    :param alpha:
    :param mld:
    :param dt:
    :param theta
    :param with_pkl: if True, add '.pkl' to the filename
    :param wave_roughness: if True, have surface roughness be wave height dependent
    :return:
    """
    name = settings.conc_dir + '{}_{}_conc_w10_{}_w_rise_{}_MLD_{}'.format(diffusion_type, boundary, w_10, w_rise, mld)
    if 'Markov' in boundary:
        name += '_alpha_list={}'.format(alpha)
    if dt != 1:
        name += '_dt={}'.format(dt)
    if diffusion_type == 'KPP':
        name += '_theta={}'.format(theta)
        if wave_roughness:
            name += '_wave_roughness'
    if with_pkl:
        name += '.pkl'
    return name


def get_data_output_name(prefix: str):
    """
    Setting the output name of the standardized field data
    :param prefix:
    :return:
    """
    return settings.data_dir + 'standardized_data_' + prefix