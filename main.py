import utils
import utils_visualization as utils_v
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
    if not utils._check_file_exist(
            utils_v.field_data_figure_names(close_up=(0, -10), wind_sort=True, norm_depth=False)):
        visualization.plot_field_data_overview(wind_sort=True, close_up=(0, -10))
    if not utils._check_file_exist(
            utils_v.field_data_figure_names(close_up=(0, -35), wind_sort=True, norm_depth=False)):
        visualization.plot_field_data_overview(wind_sort=True, close_up=(0, -35))
    if not utils._check_file_exist(utils_v.field_data_figure_names(wind_sort=True, norm_depth=True)):
        visualization.plot_field_data_overview(wind_sort=True, norm_depth=True)
    if not utils._check_file_exist(
            utils_v.field_data_figure_names(wind_sort=False, norm_depth=True, close_up=(0, -10))):
        visualization.plot_field_data_overview(wind_sort=False, norm_depth=True, close_up=(0, -10))

    """
    Parcels simulations
    """
    w_10 = [0.85, 2.4, 4.35, 6.65, 9.3]
    w_rise = [-0.03, -0.003, -0.0003]
    diffusion = 'KPP'  # 'KPP'
    boundary = 'Reflect'  # 'Reflect_Markov'
    pbar = ProgressBar()
    for wind in pbar(w_10):
        for rise in w_rise:
            # Option to remove a previous file if it exists in case I want to rerun a simulation. Setting conduct to
            # False deactivates the remove file function
            utils.remove_file(conduct=False, File=utils.get_parcels_output_name(wind, rise, diffusion, boundary))

            if not utils._check_file_exist(utils.get_parcels_output_name(wind, rise, diffusion, boundary)):
                parcels_simulation_functions.vertical_diffusion_run(wind, rise, diffusion_type=diffusion,
                                                                    boundary=boundary)
                analysis.depth_concentration(wind, rise, diffusion_type=diffusion, boundary=boundary)

            else:
                print('This simulation has already been carried out')

    """
    Visualization of simulations
    """
    # visualization.basic_profile_figure(w_10_list=w_10, w_rise_list=w_rise, selection='w_10',
    #                                    single_select=1, close_up=(0, -20), diffusion_type=diffusion,
    #                                    boundary=boundary,diffusion_curve=True)
    # visualization.basic_profile_figure(w_10_list=w_10, w_rise_list=w_rise, selection='w_rise',
    #                                    single_select=2, close_up=(0, -20), diffusion_type=diffusion,
    #                                    boundary=boundary)
    # visualization.timestep_comparison(w_10_list=[10], w_rise_list=[-0.003], selection='k_z',
    #                                   close_up=(0, -100), diffusion_type=diffusion, interval=2,
    #                                   boundary=boundary)
    # visualization.boundary_condition_comparison(w_10_list=[5], w_rise_list=[-0.003], close_up=(0, -30),
    #                                             diffusion_type=diffusion, boundary=boundary)


if __name__ == '__main__':
    run()
