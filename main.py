import utils
import parcels_simulation_functions
import analysis
import visualization
from progressbar import ProgressBar

def run():
    k_z = [1e-4]#[1e-6, 1e-5, 1e-4, 1e-3, 1e-2]
    w_10 = [15]#[0.1, 5, 10]
    w_rise = [-0.03]#[-0.03, -0.003, -0.0003]
    diffusion = 'KPP'
    boundary = 'Mixed'
    pbar=ProgressBar()
    for k in pbar(k_z):
        for wind in w_10:
            for rise in w_rise:
                # Option to remove a previous file if it exists in case I want to rerun a simulation. Setting conduct to
                # False deactivates the remove file function
                utils.remove_file(conduct=True, File=utils.get_parcels_output_name(k, wind, rise, diffusion, boundary))

                if not utils._check_file_exist(utils.get_parcels_output_name(k, wind, rise, diffusion, boundary)):
                    parcels_simulation_functions.vertical_diffusion_run(k, wind, rise, diffusion_type=diffusion,
                                                                        boundary=boundary)
                    analysis.depth_concentration(k, wind, rise, diffusion_type=diffusion, boundary=boundary)

                else:
                    print('This simulation has already been carried out')

    # Running visualization codes
    # visualization.basic_profile_figure(k_z_list=k_z, w_10_list=w_10, w_rise_list=w_rise, selection='k_z', kukulka=False,
    #                                    single_select=0, close_up=(0, -10), diffusion_type=diffusion)
    visualization.timestep_comparison(k_z_list=k_z, w_10_list=w_10, w_rise_list=w_rise, selection='k_z',
                                      close_up=(0, -100), diffusion_type=diffusion, interval=2, rouse=False,
                                      boundary=boundary)


if __name__ == '__main__':
    run()