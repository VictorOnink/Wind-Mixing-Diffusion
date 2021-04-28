import numpy as np
import settings
import utils
import progressbar


def eulerian_vertical_run(w_10, w_rise, diffusion_type):
    # This only works with buoyant particles
    w_rise = np.abs(w_rise)
    # Getting the relevant arrays for depth, Kz and concentration, where initially all concentration is at the surface
    model_arrays = get_grid_arrays(w_10, diffusion_type, w_rise)
    # Getting the grid and integration resolutions
    dZ, dt = get_grid_resolutions(model_arrays)
    # Looping through all the timesteps
    integration_number = int(settings.runtime.seconds / 2 // dt)
    pbar = progressbar.ProgressBar()
    for integration in pbar(range(integration_number)):
        model_arrays = integration_step(model_arrays, dZ, dt, w_rise)
    # Preparing the output
    output_arrays = get_output_arrays(model_arrays)
    output_name = utils.get_eulerian_output_name(w_10, w_rise, diffusion_type)
    utils.save_obj(output_name, output_arrays)


def get_grid_arrays(w_10, diffusion_type, w_rise):
    model_arrays = {}
    # Creating the arrays
    model_arrays['Z'] = np.zeros((settings.depth_levels + 2), dtype=np.float32)
    model_arrays['C'] = np.zeros((settings.depth_levels + 2), dtype=np.float32)
    model_arrays['Kz'] = np.zeros((settings.depth_levels + 2), dtype=np.float32)
    model_arrays['f_i'] = np.zeros((settings.depth_levels + 2), dtype=np.float32)
    model_arrays['g_i'] = np.zeros((settings.depth_levels + 2), dtype=np.float32)
    model_arrays['h_i'] = np.zeros((settings.depth_levels + 2), dtype=np.float32)
    model_arrays['w_r'] = np.zeros((settings.depth_levels + 2), dtype=np.float32)
    # Initializing the arrays
    model_arrays['Z'][1:-1] = np.linspace(-settings.max_depth, 0, num=settings.depth_levels)
    model_arrays['C'][-3] = 1
    model_arrays['Kz'][1:-1] = utils.get_vertical_diffusion_profile(w_10, model_arrays['Z'][1:-1], diffusion_type)
    model_arrays['w_r'][1:-1] = w_rise
    return model_arrays


def get_grid_resolutions(model_arrays):
    dZ = np.abs(model_arrays['Z'][2] - model_arrays['Z'][1])
    dt = 0.3
    return dZ, dt


def integration_step(model_arrays, dZ, dt, w_rise):
    # First, we calculate the diffusive flux at the bottom of each cell interface
    for i in range(1, settings.depth_levels + 1):
        model_arrays['f_i'][i] = model_arrays['Kz'][i] * (model_arrays['C'][i] - model_arrays['C'][i - 1]) / dZ
    # Next, we calculate the gradient of the flux
    for i in range(1, settings.depth_levels + 1):
        model_arrays['g_i'][i] = (model_arrays['f_i'][i + 1] - model_arrays['f_i'][i]) / dZ
    # Next, calculating the advective flux using an upstream scheme
    for i in range(1, settings.depth_levels + 1):
        w_r_i, w_r_ii = model_arrays['w_r'][i], model_arrays['w_r'][i + 1]
        model_arrays['h_i'][i] = 0.5 * (w_r_i + np.abs(w_r_i)) * model_arrays['C'][i - 1] + \
                                 0.5 * (w_r_i - np.abs(w_r_i)) * model_arrays['C'][i] + \
                                 - 0.5 * (w_r_ii + np.abs(w_r_ii)) * model_arrays['C'][i] + \
                                 - 0.5 * (w_r_ii - np.abs(w_r_ii)) * model_arrays['C'][i + 1]
    # # Adding the boundary condition that the flux at the bottom and surface are both equal to zero
    # model_arrays['h_i'][1], model_arrays['h_i'][-3] = 0, 0
    # Finally, we calculate the new concentration for each depth level
    model_arrays['C'] = model_arrays['C'] + dt * (model_arrays['g_i'] + model_arrays['h_i'])
    # Cheching mass conservation
    assert np.abs(sum(model_arrays['C']) - 1.0) < 0.001, "Mass isn't being adequately conserved"
    return model_arrays


def get_output_arrays(model_arrays):
    output_arrays = {'Z': model_arrays['Z'][1:-1], 'C': model_arrays['C'][1:-1]}
    # # Normalizing the concentrations
    # output_arrays['C'] = output_arrays['C'] / np.max(output_arrays['C'])
    return output_arrays
