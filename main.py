import utils
import utils_visualization as utils_v
import parcels_simulation_functions
import analysis
import visualization
import field_data
from progressbar import ProgressBar
import settings
import ubelix_submission


def run():
    """
    Field data standardization and other processing
    """
    field_data_processing()
    """
    Parcels simulations
    """
    w_10 = [6.65]  # [0.85, 2.4, 4.35, 6.65, 9.3]
    w_rise = [-0.003]  # [-0.03, -0.003, -0.0003]
    diffusion = 'KPP'
    # boundary_options = ['Mixed', 'Reflect', 'Reduce_dt', 'Mixed_Markov', 'Reflect_Markov', 'Reduce_dt_Markov']
    boundary = 'Reflect_Markov'
    pbar = ProgressBar()
    for wind in pbar(w_10):
        for rise in w_rise:
            if settings.server is 'laptop':
                # Option to remove a previous file if it exists in case I want to rerun a simulation. Setting conduct to
                # False deactivates the remove file function
                conc_file = utils.get_concentration_output_name(w_10, w_rise, diffusion, boundary)
                utils.remove_file(conduct=True, file_name=conc_file)
                if not utils.check_file_exist(conc_file):
                    parcels_simulation_functions.vertical_diffusion_run(wind, rise, diffusion_type=diffusion,
                                                                        boundary=boundary)
                    analysis.depth_concentration(wind, rise, diffusion_type=diffusion, boundary=boundary)

                else:
                    print('This simulation has already been carried out')
            elif settings.server is 'ubelix':
                ubelix_submission.ubelix_submission(diffusion, boundary, wind, rise)


    """
    Visualization of simulations
    """
    plotting(w_10, w_rise, diffusion, boundary)


def field_data_processing():
    if settings.server is 'laptop':
        field_data.data_standardization()
        if not utils.check_file_exist(
                utils_v.field_data_figure_names(close_up=(0, -10), wind_sort=True, norm_depth=False)):
            visualization.plot_field_data_overview(wind_sort=True, close_up=(0, -10))
        if not utils.check_file_exist(
                utils_v.field_data_figure_names(close_up=(0, -35), wind_sort=True, norm_depth=False)):
            visualization.plot_field_data_overview(wind_sort=True, close_up=(0, -35))
        if not utils.check_file_exist(utils_v.field_data_figure_names(wind_sort=True, norm_depth=True, close_up=(0, -1))):
            visualization.plot_field_data_overview(wind_sort=True, norm_depth=True, close_up=(0, -1))
        if not utils.check_file_exist(
                utils_v.field_data_figure_names(wind_sort=False, norm_depth=True, close_up=(0, -10))):
            visualization.plot_field_data_overview(wind_sort=False, norm_depth=True, close_up=(0, -10))
        # Determining correlations for field-measured concentrations and depth
        analysis.correlation_depth_concentration()
        analysis.range_MLD_values()


def plotting(w_10, w_rise, diffusion, boundary):
    if settings.server is 'laptop':
        # visualization.sanity_check(wind, rise, diffusion_type=diffusion, boundary=boundary)
        # Just plotting the diffusion profiles over various wind conditions
        # visualization.just_diffusion_profile(w_10_list=[0.85, 2.4, 4.35, 6.65, 9.3])

        # Plotting the profiles with different boundary conditions
        # visualization.boundary_condition_comparison(w_rise_list=[-0.003], diffusion_type='KPP', close_up=(0, -10),)
        # visualization.boundary_condition_comparison(w_rise_list=[-0.003], diffusion_type='Kukulka', close_up=(0, -10))

        # Plotting the multi-wind condition figures
        # visualization.plot_model_field_data_comparison(w_10_list=w_10, w_rise_list=w_rise, selection='w_10',
        #                                                wind_sort=True, single_select=2, close_up=(0, -20),
        #                                                diffusion_type='KPP', boundary='all')
        # visualization.plot_model_field_data_comparison(w_10_list=w_10, w_rise_list=w_rise, selection='w_10',
        #                                                wind_sort=True, single_select=2, close_up=(0, -20),
        #                                                diffusion_type='Kukulka', boundary='all')
        # visualization.plot_model_field_data_comparison(w_10_list=w_10, w_rise_list=w_rise, selection='w_10',
        #                                                wind_sort=True, single_select=2, close_up=(0, -20),
        #                                                diffusion_type='all', boundary='Reflect')

        visualization.basic_profile_figure(w_10_list=w_10, w_rise_list=w_rise, diffusion_type=diffusion,
                                           boundary=boundary, selection='w_10', single_select=0, close_up=(0, -30),
                                           diffusion_curve=False)
        # visualization.basic_profile_figure(w_10_list=w_10, w_rise_list=w_rise, diffusion_type=diffusion,
        #                                    boundary=boundary, selection='w_10', single_select=0, close_up=(0, -100),
        #                                    diffusion_curve=False)

        # Plotting for the depth profiles for multiple MLD levels
        # for beaufort in range(1, 6):
        #     visualization.mld_depth_influence(w_rise_list=w_rise, MLD_list=[15.0, 20.0, 30.0], beaufort=beaufort)
        #     visualization.mld_depth_influence(w_rise_list=w_rise, MLD_list=[15.0, 20.0, 30.0], beaufort=beaufort,
        #                                       diffusion_type='Kukulka')

        # # Plotting the depth profile for multiple timesteps
        # visualization.timestep_comparison(w_10_list=w_10, w_rise_list=w_rise,
        #                                   close_up=(0, -30), diffusion_type=diffusion, interval=1,
        #                                   boundary=boundary, mld=20.0, diffusion_curve=False)
        # visualization.timestep_comparison(w_10_list=w_10, w_rise_list=w_rise,
        #                                   close_up=(0, -100), diffusion_type=diffusion, interval=1,
        #                                   boundary=boundary, mld=20.0, diffusion_curve=False)

        # Just the field data
        # visualization.plot_field_data_overview(wind_sort=True, close_up=(0,-20))

        # Testing different values of alpha
        visualization.Markov_TL_dependence(w_rise_list=w_rise, single_select=0, close_up=(0, -30),
                                           diffusion_type=diffusion)


if __name__ == '__main__':
    run()
