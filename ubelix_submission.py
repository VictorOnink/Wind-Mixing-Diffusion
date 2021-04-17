import settings
import os
import parcels_simulation_functions
import analysis
import utils
import eulerian_simulation_functions


def ubelix_submission(diffusion, boundary, wind, rise, alpha, submission):
    """
    Creating a bash script that can be used to submit a simulation on the ubelix server
    """
    # Creating the submission file
    file_name = settings.bin_dir + r'ubelix_submission.sh'
    print(file_name)
    file = open(file_name, "w")
    file.write("#!/bin/sh\n")
    file.write("#SBATCH --mail-type=begin,end,fail\n")
    file.write("#SBATCH --mail-user=victor.onink@climate.unibe.ch\n")
    file.write("#SBATCH --job-name={}\n".format(run_name(diffusion, boundary, wind, rise, alpha)))
    file.write("#SBATCH --output=bin/{}.o%j\n".format(run_name(diffusion, boundary, wind, rise, alpha)))
    file.write("#SBATCH --mem-per-cpu=6G\n")
    if submission == 'parcels':
        file.write("#SBATCH --time=00:30:00\n")
    elif submission == 'eulerian':
        file.write("#SBATCH --time=01:00:00\n")
    file.write('#SBATCH --partition=all\n')
    file.write('source /storage/homefs/vo18e689/.bash_profile\n')
    file.write('source /storage/homefs/vo18e689/anaconda3/bin/activate py3_parcels\n')
    s
    file.write('cd "{}"\n'.format(settings.code_dir))
    var_dict = {'diffusion': diffusion, 'boundary': boundary, 'wind': wind, 'rise': rise, 'alpha_list': alpha,
                'submission': submission}
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


def run_name(diffusion, boundary, wind, rise, alpha):
    if 'Markov' in boundary:
        name = '{}_{}_wind={}_rise={}_alpha={}'.format(diffusion, boundary, wind, rise, alpha)
    else:
        name = '{}_{}_wind={}_rise={}'.format(diffusion, boundary, wind, rise)
    return name


def ubelix_synchronization(update: bool = False):
    """ Download any concentration files computed on the ubelix server to my laptop """
    update_dict = {'parcels': [settings.root_direc['ubelix'] + 'concentration_output/', settings.conc_dir],
                   'eulerian': [settings.root_direc['ubelix'] + 'eulerian_output/', settings.eulout_dir]}
    if settings.server is 'laptop':
        if update:
            for ubelix_folder in update_dict:
                current_folder = os.getcwd()
                os.chdir(update_dict[ubelix_folder][1])
                command = 'rsync -av ubelix:{}* .'.format(update_dict[ubelix_folder][0])
                os.system(command)
                os.chdir(current_folder)


if __name__ == '__main__':
    submission = os.getenv('submission')
    diffusion, boundary = os.getenv('diffusion'), os.getenv('boundary')
    wind, rise, alpha = float(os.getenv('wind')), float(os.getenv('rise')), float(os.getenv('alpha_list'))
    if submission == 'parcels':
        parcels_simulation_functions.vertical_diffusion_run(wind, rise, alpha=alpha, diffusion_type=diffusion,
                                                            boundary=boundary)
        analysis.depth_concentration(w_10=wind, w_rise=rise, alpha=alpha, diffusion_type=diffusion, boundary=boundary)
    if submission == 'eulerian':
        eulerian_simulation_functions.eulerian_vertical_run(w_10=wind, w_rise=rise, diffusion_type=diffusion)
