import utils
import parcels_simulation_functions
import analysis
import visualization
import field_data
from progressbar import ProgressBar


def run():
    """
    Field data standardization
    """
    field_data.data_standardization()

    """
    Parcels simulations
    """
    # k_z = [1e-4]#[1e-4, 1e-4, 1e-4]
    # w_10 = [5]#[2, 5, 10]
    # w_rise = [-0.003]#[-0.03, -0.003, -0.0003]
    # diffusion = 'KPP'
    # boundary = 'Reflect_Markov'
    # pbar = ProgressBar()
    # for k in pbar(k_z):
    #     for wind in w_10:
    #         for rise in w_rise:
    #             # Option to remove a previous file if it exists in case I want to rerun a simulation. Setting conduct to
    #             # False deactivates the remove file function
    #             utils.remove_file(conduct=False, File=utils.get_parcels_output_name(k, wind, rise, diffusion, boundary))
    #
    #             if not utils._check_file_exist(utils.get_parcels_output_name(k, wind, rise, diffusion, boundary)):
    #                 parcels_simulation_functions.vertical_diffusion_run(k, wind, rise, diffusion_type=diffusion,
    #                                                                     boundary=boundary)
    #                 analysis.depth_concentration(k, wind, rise, diffusion_type=diffusion, boundary=boundary)
    #
    #             else:
    #                 print('This simulation has already been carried out')

    """
    Visualization of simulations
    """
    # visualization.basic_profile_figure(k_z_list=k_z, w_10_list=w_10, w_rise_list=w_rise, selection='w_10',
    #                                    kukulka=True, single_select=1, close_up=(0, -100), diffusion_type=diffusion,
    #                                    boundary=boundary)
    # visualization.basic_profile_figure(k_z_list=k_z, w_10_list=w_10, w_rise_list=w_rise, selection='w_10',
    #                                    kukulka=True, single_select=1, close_up=(0, -30), diffusion_type=diffusion,
    #                                    boundary=boundary)
    # visualization.basic_profile_figure(k_z_list=k_z, w_10_list=w_10, w_rise_list=w_rise, selection='w_rise',
    #                                    kukulka=True, single_select=2, close_up=(0, -100), diffusion_type=diffusion,
    #                                    boundary=boundary)
    # visualization.basic_profile_figure(k_z_list=k_z, w_10_list=w_10, w_rise_list=w_rise, selection='w_rise',
    #                                    kukulka=True, single_select=2, close_up=(0, -30), diffusion_type=diffusion,
    #                                    boundary=boundary)
    # visualization.timestep_comparison(k_z_list=[1e-4], w_10_list=[10], w_rise_list=[-0.003], selection='k_z',
    #                                   close_up=(0, -100), diffusion_type=diffusion, interval=2,
    #                                   boundary=boundary)
    # visualization.timestep_comparison(k_z_list=[1e-4], w_10_list=[10], w_rise_list=[-0.003], selection='k_z',
    #                                   close_up=(0, -30), diffusion_type=diffusion, interval=2,
    #                                   boundary=boundary)
    # visualization.boundary_condition_comparison(k_z_list=[1e-4], w_10_list=[5], w_rise_list=[-0.003], close_up=(0, -30),
    #                                             diffusion_type=diffusion, boundary=boundary)

if __name__ == '__main__':
    run()