import settings
import os
import parcels_simulation_functions
import analysis
import utils


def ubelix_submission(diffusion, boundary, wind, rise):
    """
    Creating a bash script that can be used to submit a simulation on the ubelix server
    """
    # Creating the submission file
    file_name = settings.bin_dir + r'ubelix_submission.sh'
    print(file_name)
    file = open(file_name, "w+")
    file.write("#!/bin/sh\n")
    file.write("#SBATCH --mail-type=begin,end,fail\n")
    file.write("#SBATCH --mail-user=victor.onink@climate.unibe.ch\n")
    file.write("#SBATCH --job-name={}\n".format(run_name(diffusion, boundary, wind, rise)))
    file.write("#SBATCH --output=bin/{}.o%j\n".format(run_name(diffusion, boundary, wind, rise)))
    file.write("#SBATCH --mem-per-cpu=6G\n")
    file.write("#SBATCH --time=00:10:00\n")
    file.write('#SBATCH --partition=debug\n')
    file.write('source /home/ubelix/climate/vo18e689/.bash_profile\n')
    file.write('source /home/ubelix/climate/vo18e689/anaconda3/bin/activate py3_parcels_v2_2\n')
    file.write('cd "{}"\n'.format(settings.code_dir))
    var_dict = {'diffusion': diffusion, 'boundary': boundary, 'wind': wind, 'rise': rise}
    export_var(file, var_dict)
    file.write('python ubelix_submission.py -p 10 -v')
    file.close()
    # Submitting the actual job
    os.system('sbatch {}'.format(file_name))
    # Removing the submission file
    utils.remove_file(conduct=True, file_name=file_name)


def export_var(file, var_dict):
    for var in var_dict:
        file.write("{}={}\n".format(var, var_dict[var]))
        file.write("export {}\n".format(var))


def run_name(diffusion, boundary, wind, rise):
    name = '{}_{}_wind={}_rise={}'.format(diffusion, boundary, wind, rise)
    return name


if __name__ == '__main__':
    diffusion, boundary = os.getenv('diffusion'), os.getenv('boundary')
    wind, rise = float(os.getenv('wind')), float(os.getenv('rise'))
    parcels_simulation_functions.vertical_diffusion_run(wind, rise, diffusion_type=diffusion, boundary=boundary)
    analysis.depth_concentration(wind, rise, diffusion_type=diffusion, boundary=boundary)