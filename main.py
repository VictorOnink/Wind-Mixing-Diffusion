import settings as SET
import parcels_simulation_functions
import analysis
import visualization
from progressbar import ProgressBar

def run():
    k_z = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2]
    w_10 = [0.1, 10, 15, 20]
    w_rise = [-0.03, -0.003, -0.0003]
    pbar=ProgressBar()
    for k in pbar(k_z):
        for wind in w_10:
            for rise in w_rise:
                print('Going through all parameter combinations')
                parcels_simulation_functions.vertical_diffusion_run(k, wind, rise)
                analysis.depth_concentration(k, wind, rise)
    visualization.basic_profile_figure(k_z_list=k_z, w_10_list=w_10, w_rise_list=w_rise, selection='k_z')
    visualization.basic_profile_figure(k_z_list=k_z, w_10_list=w_10, w_rise_list=w_rise, selection='k_z',
                                       close_up=(0, -10))
    visualization.basic_profile_figure(k_z_list=k_z, w_10_list=w_10, w_rise_list=w_rise, selection='w_10')
    visualization.basic_profile_figure(k_z_list=k_z, w_10_list=w_10, w_rise_list=w_rise, selection='w_10',
                                       close_up=(0, -10))
    visualization.basic_profile_figure(k_z_list=k_z, w_10_list=w_10, w_rise_list=w_rise, selection='w_rise')
    visualization.basic_profile_figure(k_z_list=k_z, w_10_list=w_10, w_rise_list=w_rise, selection='w_rise',
                                       close_up=(0, -10))


if __name__ == '__main__':
    run()