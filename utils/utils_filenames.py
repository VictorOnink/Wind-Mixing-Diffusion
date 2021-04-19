import settings


def get_parcels_output_name(w_10, w_rise, diffusion_type, boundary, alpha, mld=settings.MLD):
    name = settings.output_dir + diffusion_type + '_' + boundary + '_w10_{}_w_rise_{}_MLD_{}'.format(w_10, w_rise,
                                                                                                        mld)
    if 'Markov' in boundary:
        name += 'alpha_list={}'.format(alpha)
    return name + '.nc'


def get_eulerian_output_name(w_10, w_rise, diffusion_type, mld=settings.MLD):
    str_format = diffusion_type, w_rise, w_10, mld
    name = settings.eulout_dir + 'eulerian_{}_w_r={}_w_10={}_MLD={}'.format(*str_format)
    return name


def get_concentration_output_name(w_10, w_rise, diffusion_type, boundary, alpha, mld=settings.MLD):
    name = settings.conc_dir + diffusion_type + '_' + boundary + '_conc_w10_{}_w_rise_{}_MLD_{}'.format(w_10, w_rise,
                                                                                                        mld)
    if 'Markov' in boundary:
        name += 'alpha_list={}'.format(alpha)
    return name


def get_data_output_name(prefix: str):
    return settings.data_dir + 'standardized_data_' + prefix