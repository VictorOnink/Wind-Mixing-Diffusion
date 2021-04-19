import os
import pickle

import numpy as np


def save_obj(filename, item):
    with open(filename + '.pkl', 'wb') as f:
        pickle.dump(item, f, pickle.HIGHEST_PROTOCOL)


def load_obj(filename):
    with open(filename + '.pkl', 'rb') as f:
        return pickle.load(f)


def check_file_exist(file_name: str):
    return os.path.isfile(file_name)


def remove_file(conduct: bool, file_name: str):
    if conduct:
        if check_file_exist(file_name):
            os.remove(file_name)


def find_nearest_index(depth, z_ref):
    return (np.abs(depth - z_ref)).argmin()


def exclude_field_data(exclude, sources):
    if exclude is not None:
        for source in exclude:
            sources.remove(source)
    return sources