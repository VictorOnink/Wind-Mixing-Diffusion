import settings as SET
from parcels_simulation_functions import vertical_diffusion_run

def run():
    k_z = 1e-5
    vertical_diffusion_run(k_z)

if __name__ == '__main__':
    run()