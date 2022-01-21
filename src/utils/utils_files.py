import os
import pickle
import numpy as np


def save_obj(filename, item):
    """
    Saving a item using pickle
    :param filename:
    :param item:
    :return:
    """
    with open(filename + '.pkl', 'wb') as f:
        pickle.dump(item, f, pickle.HIGHEST_PROTOCOL)


def load_obj(filename):
    """
    Loading the data that had been previously pickled using save_obj()
    :param filename:
    :return:
    """
    with open(filename + '.pkl', 'rb') as f:
        return pickle.load(f)


def check_file_exist(file_name: str):
    """
    Check if the file called file_name already exists
    :param file_name:
    :return:
    """
    return os.path.isfile(file_name)


def remove_file(conduct: bool, file_name: str):
    """
    Removing a file, but only if conduct is True and the file exists
    :param conduct:
    :param file_name:
    :return:
    """
    if conduct:
        if check_file_exist(file_name):
            os.remove(file_name)


def find_nearest_index(depth, z_ref):
    """
    Return the index of the element of the deoth array that is closest to the value of z_ref
    :param depth:
    :param z_ref:
    :return:
    """
    return (np.abs(depth - z_ref)).argmin()


def exclude_field_data(exclude, sources):
    """
    The option to remove field data from analysis.
    :param exclude: a list containing the field data that is to be excluded
    :param sources: a list containing all field data sources
    :return:
    """
    if exclude is not None:
        for source in exclude:
            sources.remove(source)
    return sources


def simulation_description(wind, rise, alpha, theta, gamma, diffusion, wave_roughness, boundary):
    if diffusion == 'KPP':
        str_format = boundary, wind, rise, theta, wave_roughness
        statement = 'Diffusion KPP, boundary {}, w_10={}, rise={}, theta={}, wave_roughness={}'.format(*str_format)
    elif diffusion == 'SWB':
        str_format = boundary, wind, rise, gamma
        statement = 'Diffusion SWB, boundary {}, w_10={}, rise={}, gamma={}'.format(*str_format)
    if 'Markov' in boundary:
        statement += ', alpha={}'.format(alpha)
    print(statement)

