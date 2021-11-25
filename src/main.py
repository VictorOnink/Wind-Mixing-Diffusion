import analysis, utils, visualization, field_data, parcels_simulation_functions
from progressbar import ProgressBar

w_10 = [0.85, 2.4, 4.35, 6.65, 9.3]  # [0.85, 2.4, 4.35, 6.65, 9.3]
w_rise = [-0.03, -0.003, -0.0003]
alpha = [0.0]  # [0.0, 0.1, 0.3, 0.5, 0.7, 0.95]
diffusion = 'SWB'
boundary = 'Ceiling'


def parcels_simulations(wind, rise, alpha):
    # Option to remove a previous file if it exists in case I want to rerun a simulation. Setting
    # conduct to False deactivates the remove file function
    concentration_file = utils.get_concentration_output_name(wind, rise, diffusion, boundary, alpha=alpha) + '.pkl'
    utils.remove_file(conduct=True, file_name=concentration_file)
    if not utils.check_file_exist(concentration_file):
        parcels_simulation_functions.vertical_diffusion_run(wind, rise, diffusion_type=diffusion,
                                                            boundary=boundary, alpha=alpha)
        analysis.depth_concentration(wind, rise, diffusion_type=diffusion, boundary=boundary,
                                     alpha=alpha, remove_file=True)
    else:
        print('This simulation has already been carried out')
    # Print the RMSE error between the calculated concentration profile and the field measurements
    analysis.determine_RMSE(wind, rise, diffusion, boundary, alpha, conduct=False)
    analysis.depth_bin_numbers(wind, rise, diffusion, boundary, alpha, conduct=False)
    analysis.compute_modelling_efficiency(wind, rise, diffusion, boundary, alpha, conduct=False)
    analysis.correlation_field_model_data(wind, rise, diffusion, boundary, alpha, conduct=False)


def field_data_processing():
    field_data.data_standardization()
    if not utils.check_file_exist(
            visualization.field_data_figure_names(close_up=(0, -10), wind_sort=True, norm_depth=False)):
        visualization.plot_field_data_overview(wind_sort=True, close_up=(0, -10))
    if not utils.check_file_exist(
            visualization.field_data_figure_names(close_up=(0, -35), wind_sort=True, norm_depth=False)):
        visualization.plot_field_data_overview(wind_sort=True, close_up=(0, -35))
    if not utils.check_file_exist(
            visualization.field_data_figure_names(wind_sort=True, norm_depth=True, close_up=(0, -1))):
        visualization.plot_field_data_overview(wind_sort=True, norm_depth=True, close_up=(0, -1))
    if not utils.check_file_exist(
            visualization.field_data_figure_names(wind_sort=False, norm_depth=True, close_up=(0, -10))):
        visualization.plot_field_data_overview(wind_sort=False, norm_depth=True, close_up=(0, -10))
    # Determining correlations for field-measured concentrations and depth
    analysis.correlation_depth_concentration(conduct=False)
    analysis.range_MLD_values(conduct=False)


def plotting():
    # Just plotting the diffusion profiles over various wind conditions
    # visualization.just_diffusion_profile(w_10_list=[0.85, 2.4, 4.35, 6.65, 9.3])

    # Plotting the profiles with different boundary conditions
    # visualization.boundary_condition_comparison(w_rise_list=[-0.003], alpha_list=alpha, close_up=(0, -10),
    #                                             diffusion_type='KPP')
    # visualization.boundary_condition_comparison(w_rise_list=[-0.003], alpha_list=alpha, close_up=(0, -10),
    #                                             diffusion_type='SWB')
    # visualization.boundary_condition_comparison(w_rise_list=[-0.03, -0.003, -0.0003], alpha_list=alpha,
    #                                             close_up=(0, -20), diffusion_type='KPP')
    # visualization.boundary_condition_comparison(w_rise_list=[-0.03, -0.003, -0.0003], alpha_list=alpha,
    #                                             close_up=(0, -20), diffusion_type='SWB')
    # for Bft in [1, 2, 3, 4, 5]:
    #     visualization.multiple_boundary_condition_comparison(close_up=(0, -20), beaufort=Bft)

    # Plotting the multi-wind condition figures
    # visualization.plot_model_field_data_comparison(w_10_list=w_10, w_rise_list=w_rise, alpha_list=alpha,
    #                                                selection='w_10', single_select=2, wind_sort=True,
    #                                                close_up=(0, -20), diffusion_type='all', boundary='Ceiling')

    # Just the field data
    # visualization.plot_field_data_overview(wind_sort=True, close_up=(0, -25))

    # Testing different values of alpha_list
    # visualization.markov_alpha_dependence(w_rise_list=w_rise, single_select=0, close_up=(0, -30),
    #                                       diffusion_type=diffusion)

    # Comparing Markov with different alpha values, but then showing more than one type of diffusion
    # visualization.diffusion_markov_comparison(w_rise_list=[-0.003], single_select=0, close_up=(0, -25))

    # Creating a figure to compare the RMSE values for the Markov-0 runs
    # visualization.markov_0_RMSE_comparison()

    # Creating a similar figure to compare the RMSE values for the Markov-0 runs
    # visualization.markov_1_RMSE_comparison()

    # The influence of the integration timestep for Markov 0
    # visualization.multiple_integration_timestep_control(beaufort=4, boundary='Reflect')
    # visualization.multiple_integration_timestep_control(beaufort=4, boundary='Ceiling')

    pass


if __name__ == '__main__':
    """
    Field data standardization and other processing
    """
    field_data_processing()
    """
    Parcels simulations
    """
    pbar = ProgressBar()
    for wind in pbar(w_10):
        for rise in w_rise:
            for alpha_val in alpha:
                parcels_simulations(wind=wind, rise=rise, alpha=alpha_val)
    analysis.timestep_dependent_RMSE(conduct=False, diffusion_type='SWB')
    analysis.timestep_dependent_RMSE(conduct=False, diffusion_type='KPP')
    """
    Visualization of simulations
    """
    plotting()
