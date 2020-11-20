import settings as SET
from parcels_simulation_functions import vertical_diffusion_run
from analysis import depth_concentration
from progressbar import ProgressBar

def run():
    k_z = [1e-6, 1e-5, 1e-4]
    w_10 = [0.1, 10, 15, 20]
    w_rise = [-0.03, -0.003, -0.0003]
    pbar=ProgressBar()
    for k in pbar(k_z):
        for wind in w_10:
            for rise in w_rise:
                vertical_diffusion_run(k, wind, rise)
    #depth_concentration(k_z, w_10, w_rise)


if __name__ == '__main__':
    run()