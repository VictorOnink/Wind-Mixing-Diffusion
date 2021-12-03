import analysis
import utils
import visualization
import field_data
import parcels_simulation_functions
from progressbar import ProgressBar

w_10 = [0.85, 2.4, 4.35, 6.65, 9.3]  # [0.85, 2.4, 4.35, 6.65, 9.3]
w_rise = [-0.03, -0.003, -0.0003]  # [-0.03, -0.003, -0.0003]
alpha = [0.0]  # [0.0, 0.1, 0.3, 0.5, 0.7, 0.95]
theta = [1.0]  # [1.0, 2.0, 3.0, 4.0, 5.0]
diffusion = 'KPP'
wave_roughness = False
boundary = 'Ceiling'

if diffusion == 'SWB':
    assert theta.__len__() == 1, "Theta is not a relevant parameter for SWB diffusion, don't run multiple simulations!"


def parcels_simulations(wind, rise, alpha, theta):
    """
    Conducting the parcels simulation
    :param wind: surface wind speed
    :param rise: particle rise velocity
    :param alpha: for M-1 simulations, the value of the memory term alpha
    :param theta: for KPP diffusion, the langmuir circulation amplification factor
    :return:
    """
    utils.simulation_description(wind, rise, alpha, theta, diffusion, wave_roughness, boundary)
    # Option to remove a previous file if it exists in case I want to rerun a simulation. Setting
    # conduct to False deactivates the remove file function
    concentration_file = utils.get_concentration_output_name(wind, rise, diffusion, boundary, alpha=alpha, theta=theta,
                                                             with_pkl=True, wave_roughness=wave_roughness)
    utils.remove_file(conduct=False, file_name=concentration_file)
    if not utils.check_file_exist(concentration_file):
        parcels_simulation_functions.vertical_diffusion_run(wind, rise, diffusion_type=diffusion, boundary=boundary,
                                                            alpha=alpha, wave_roughness=wave_roughness, theta=theta)
        analysis.depth_concentration(wind, rise, diffusion_type=diffusion, boundary=boundary,
                                     alpha=alpha, theta=theta, remove_file=True, wave_roughness=wave_roughness)
    else:
        print('This simulation has already been carried out')
    # Calculate the RMSE error between the calculated concentration profile and the field measurements
    analysis.determine_RMSE(wind, rise, diffusion, boundary, alpha, theta=theta, conduct=False,
                            wave_roughness=wave_roughness)
    # Calculate the percentage of particles within the prescribed depth bins
    analysis.depth_bin_numbers(wind, rise, diffusion, boundary, alpha, theta=theta, conduct=False,
                               wave_roughness=wave_roughness)
    # Calculate the correlation between the field data and the modelled distribution
    analysis.correlation_field_model_data(wind, rise, diffusion, boundary, alpha, theta=theta, conduct=False,
                                          wave_roughness=wave_roughness)
    # Calculate the mean percentage deviation between the field data and modelled distribution
    analysis.mean_percentage_deviation(wind=wind, rise=rise, diffusion=diffusion, boundary=boundary, alpha=alpha,
                                       theta=theta, wave_roughness=wave_roughness, conduct=False)
    # Calculate the w_r/w'_max ratio for each wind condition and rise velocity
    analysis.rise_velocity_turbulence_ratio(wind, rise, diffusion, theta=theta, conduct=False,
                                            wave_roughness=wave_roughness)


def field_data_processing():
    # Carry out standardization functions
    field_data.data_standardization()
    # Determining correlations for field-measured concentrations and depth
    analysis.correlation_depth_concentration(conduct=False)
    analysis.range_MLD_values(conduct=False)


def plotting():
    # Figure 1: Plotting just the diffusion profiles for various wind conditions
    # visualization.just_diffusion_profile(w_10_list=[0.85, 2.4, 4.35, 6.65, 9.3])
    # visualization.just_diffusion_profile(w_10_list=[0.85, 2.4, 4.35, 6.65, 9.3], with_theta=False, swb=False,
    #                                      wave_roughness=True)

    # Figure 2: Plotting vertical profiles for KPP and SWB diffusion under various wind conditions
    visualization.plot_model_field_data_comparison(w_10_list=w_10, w_rise_list=w_rise, alpha_list=alpha,
                                                   selection='w_10', single_select=2, wind_sort=True,
                                                   close_up=(0, -20), diffusion_type='all', boundary='Ceiling',
                                                   theta=1.0, wave_roughness=False)

    # Figure 3: Influence of the Langmuir amplification factor and z_0 for KPP
    # visualization.theta_langmuir_wave_roughness_sensitivity(w_10_list=w_10, w_rise=-0.003,
    #                                                         theta_list=[1.0, 2.0, 3.0, 4.0, 5.0], close_up=(0, -20))

    # Figure 4: Root mean square error (RMSE) for the M-0 model runs compared with observations
    # visualization.markov_0_RMSE_comparison()

    # Figure 5: Comparing M-0 with M-1 models with various alpha values for both KPP and SWB diffusion
    # visualization.diffusion_markov_comparison(w_rise_list=[-0.003], single_select=0, close_up=(0, -25), theta=1.0,
    #                                           wave_roughness=False)

    # Figure 6: RMSE for the M-1 model runs compared with observations
    # visualization.markov_1_RMSE_comparison()

    # Supplementary Figures: Influence of dt for M-0 runs for Reflect and Ceiling conditions
    # visualization.multiple_integration_timestep_control(beaufort=4, boundary='Reflect')
    # visualization.multiple_integration_timestep_control(beaufort=4, boundary='Reflect', wave_roughness=True)
    # visualization.multiple_integration_timestep_control(beaufort=4, boundary='Ceiling')
    # visualization.multiple_integration_timestep_control(beaufort=4, boundary='Ceiling', wave_roughness=True)

    # Supplementary Figures: Influence of boundary conditions on overall profile
    # visualization.multiple_boundary_condition_comparison(close_up=(0, -20), beaufort=4)
    # visualization.multiple_boundary_condition_comparison(close_up=(0, -20), beaufort=4, wave_roughness=True)

    # Supplementary Figure: Time evolution of the vertical profiles for M-0 for both KPP and SWB diffusion
    # visualization.timestep_comparison(close_up=(0, -25))

    # Supplementary Figure: Influence of the theta Langmuir amplification factor
    # w_10 = [0.85, 2.4, 4.35, 6.65, 9.3]
    # visualization.theta_langmuir_wave_roughness_sensitivity(w_10_list=w_10, w_rise=-0.03,
    #                                                         theta_list=[1.0, 2.0, 3.0, 4.0, 5.0], close_up=(0, -20))
    # visualization.theta_langmuir_wave_roughness_sensitivity(w_10_list=w_10, w_rise=-0.0003,
    #                                                         theta_list=[1.0, 2.0, 3.0, 4.0, 5.0],
    #                                                         close_up=(0, -20), with_observations=False)

    # Supplementary figure: plotting the correlation and RMSE between field measurements and modelled profiles
    # visualization.plot_error_correlaton(markov=False, diffusion='KPP')
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
                for langmuir in theta:
                    parcels_simulations(wind=wind, rise=rise, alpha=alpha_val, theta=langmuir)
    analysis.timestep_dependent_RMSE(conduct=False, diffusion_type='SWB')
    analysis.timestep_dependent_RMSE(conduct=False, diffusion_type='KPP')
    """
    Visualization of simulations
    """
    plotting()
