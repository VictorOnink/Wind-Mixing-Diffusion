import utils
import utils_visualization as utils_v
import parcels_simulation_functions
import analysis
import visualization
import field_data
from progressbar import ProgressBar
import settings
import ubelix_submission
import eulerian_simulation_functions

w_10 = [0.85, 2.4, 4.35, 6.65, 9.3]
w_rise = [-0.03, -0.003, -0.0003]
alpha = [0.0]  # [0.0, 0.1, 0.3, 0.5, 0.7, 0.95]
diffusion = 'KPP'
# boundary_options = ['Mixed', 'Reflect', 'Reduce_dt', 'Mixed_Markov', 'Reflect_Markov', 'Reduce_dt_Markov']
boundary = 'Reflect'


def parcels_simulations(wind, rise, alpha):
    if settings.server is 'laptop':
        # Option to remove a previous file if it exists in case I want to rerun a simulation. Setting
        # conduct to False deactivates the remove file function
        concentration_file = utils.get_concentration_output_name(wind, rise, diffusion, boundary, alpha=alpha) + '.pkl'
        utils.remove_file(conduct=False, file_name=concentration_file)
        if not utils.check_file_exist(concentration_file):
            parcels_simulation_functions.vertical_diffusion_run(wind, rise, diffusion_type=diffusion,
                                                                boundary=boundary, alpha=alpha)
            analysis.depth_concentration(wind, rise, diffusion_type=diffusion, boundary=boundary,
                                         alpha=alpha)
        else:
            print('This simulation has already been carried out')
        analysis.determine_RMSE(wind, rise, diffusion, boundary, alpha)
    elif settings.server is 'ubelix':
        # ubelix_submission.ubelix_submission(diffusion, boundary, wind, rise, alpha, submission='parcels')
        ubelix_submission.ubelix_submission(diffusion, boundary, wind, rise, alpha, submission='eulerian')


def field_data_processing():
    if settings.server is 'laptop':
        field_data.data_standardization()
        if not utils.check_file_exist(
                utils_v.field_data_figure_names(close_up=(0, -10), wind_sort=True, norm_depth=False)):
            visualization.plot_field_data_overview(wind_sort=True, close_up=(0, -10))
        if not utils.check_file_exist(
                utils_v.field_data_figure_names(close_up=(0, -35), wind_sort=True, norm_depth=False)):
            visualization.plot_field_data_overview(wind_sort=True, close_up=(0, -35))
        if not utils.check_file_exist(
                utils_v.field_data_figure_names(wind_sort=True, norm_depth=True, close_up=(0, -1))):
            visualization.plot_field_data_overview(wind_sort=True, norm_depth=True, close_up=(0, -1))
        if not utils.check_file_exist(
                utils_v.field_data_figure_names(wind_sort=False, norm_depth=True, close_up=(0, -10))):
            visualization.plot_field_data_overview(wind_sort=False, norm_depth=True, close_up=(0, -10))
        # Determining correlations for field-measured concentrations and depth
        analysis.correlation_depth_concentration()
        analysis.range_MLD_values()


def plotting():
    if settings.server is 'laptop':
        # visualization.sanity_check(wind, rise, diffusion_type=diffusion, boundary=boundary, alpha_list=alpha)
        # Just plotting the diffusion profiles over various wind conditions
        # visualization.just_diffusion_profile(w_10_list=[0.85, 2.4, 4.35, 6.65, 9.3])

        # Plotting the profiles with different boundary conditions
        # visualization.boundary_condition_comparison(w_rise_list=[-0.003], alpha_list=alpha, close_up=(0, -10),
        #                                             diffusion_type='KPP')
        # visualization.boundary_condition_comparison(w_rise_list=[-0.003], alpha_list=alpha, close_up=(0, -10),
        #                                             diffusion_type='Kukulka')

        # Plotting the multi-wind condition figures
        # visualization.plot_model_field_data_comparison(w_10_list=w_10, w_rise_list=w_rise, alpha_list=alpha,
        #                                                selection='w_10', single_select=2, wind_sort=True,
        #                                                close_up=(0, -20), diffusion_type='KPP', boundary='Reflect')
        # visualization.plot_model_field_data_comparison(w_10_list=w_10, w_rise_list=w_rise, alpha_list=alpha,
        #                                                selection='w_10', single_select=2, wind_sort=True,
        #                                                close_up=(0, -20), diffusion_type='Kukulka', boundary='all')
        # visualization.plot_model_field_data_comparison(w_10_list=w_10, w_rise_list=w_rise, alpha_list=alpha,
        #                                                selection='w_10', single_select=2, wind_sort=True,
        #                                                close_up=(0, -20), diffusion_type='all', boundary='Reflect')
        # visualization.plot_model_field_data_comparison(w_10_list=[6.65], w_rise_list=[-0.03, -0.003], alpha_list=alpha,
        #                                                selection='w_rise', single_select=0, wind_sort=False,
        #                                                close_up=(0, -20), diffusion_type='all', boundary='Reflect',
        #                                                fig_size=(8, 8))
        # visualization.plot_model_field_data_comparison(w_10_list=w_10, w_rise_list=w_rise, alpha_list=[0.95],
        #                                                selection='w_10', single_select=2, wind_sort=True,
        #                                                close_up=(0, -20), diffusion_type='all',
        #                                                boundary='Reflect_Markov')

        # visualization.basic_profile_figure(w_10_list=w_10, w_rise_list=w_rise, alpha_list=alpha, selection='w_10',
        #                                    close_up=(0, -30), single_select=0, diffusion_type=diffusion,
        #                                    boundary=boundary, diffusion_curve=False)
        # visualization.basic_profile_figure(w_10_list=w_10, w_rise_list=w_rise, alpha_list=alpha,
        #                                    diffusion_type=diffusion,
        #                                    boundary=boundary, selection='w_10', single_select=0, close_up=(0, -100),
        #                                    diffusion_curve=False)

        # Plotting for the depth profiles for multiple MLD levels
        # for beaufort in range(1, 6):
        #     visualization.mld_depth_influence(w_rise_list=w_rise, MLD_list=[15.0, 20.0, 30.0], alpha_list=alpha,
        #                                       beaufort=beaufort)
        #     visualization.mld_depth_influence(w_rise_list=w_rise, MLD_list=[15.0, 20.0, 30.0], alpha_list=alpha,
        #                                       beaufort=beaufort, diffusion_type='Kukulka')

        # # Plotting the depth profile for multiple timesteps
        # visualization.timestep_comparison(w_10_list=[0.85], w_rise_list=[-0.0003], alpha_list=[0.0], close_up=(0, -30),
        #                                   mld=20.0, diffusion_type=diffusion, interval=1, boundary=boundary,
        #                                   diffusion_curve=False)
        # visualization.timestep_comparison(w_10_list=w_10, w_rise_list=w_rise, alpha_list=alpha, close_up=(0, -100),
        #                                   mld=20.0, diffusion_type=diffusion, interval=1, boundary=boundary,
        #                                   diffusion_curve=False)

        # Just the field data
        # visualization.plot_field_data_overview(wind_sort=True, close_up=(0, -25))

        # Testing different values of alpha_list
        # visualization.Markov_alpha_dependence(w_rise_list=w_rise, single_select=0, close_up=(0, -30),
        #                                       diffusion_type=diffusion)

        # Comparing Markov with different alpha values, but then showing more than one type of diffusion
        # visualization.diffusion_markov_comparison(w_rise_list=[-0.03], single_select=0, close_up=(0, -25))

        # Creating a figure to compare the RMSE values for the Markov-0 runs
        # visualization.markov_0_RMSE_comparison()

        # Creating a similar figure to compare the RMSE values for the Markov-0 runs
        # visualization.markov_1_RMSE_comparison()
        pass


if __name__ == '__main__':
    """
    Field data standardization and other processing
    """
    field_data_processing()
    """
    Synchronize with the data stored on the ubelix server
    """
    ubelix_submission.ubelix_synchronization(update=False)
    """
    Parcels simulations
    """
    pbar = ProgressBar()
    for wind in pbar(w_10):
        for rise in w_rise:
            for alpha_val in alpha:
                parcels_simulations(wind=wind, rise=rise, alpha=alpha_val)
    """
    Visualization of simulations
    """
    plotting()
