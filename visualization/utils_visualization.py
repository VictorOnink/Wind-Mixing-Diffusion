import settings
import numpy as np
import matplotlib.pyplot as plt
import utils
from matplotlib.gridspec import GridSpec


def diffusion_curve_axis(ax, ax_label_size, w_10, profile_dict, diffusion_type, color, linestyle='dotted',
                         gradient=False):
    """
    this creates and returns the axis object on which we plot the diffusion curve
    :param ax: axis object to which we add the diffusion axis
    :param ax_label_size: fontsize of the axis labels
    :param w_10: 10m wind speed
    :param profile_dict: dictionary containing the depth array that we use to calculate the diffusion profile
    :param diffusion_type: SWB or KPP
    :param color: color of the diffusion profile line
    :param linestyle: the linestyle of the diffusion profile
    :param gradient: boolean statement, where if gradient == True, we also plot the diffusion profile vertical gradient
    :return:
    """
    # Creating the axis for the diffusion curve, which will have the same x axis as ax but have an opposing y axis
    ax2 = ax.twiny()
    # Get the proper diffusion curve
    depth = profile_dict['depth_bins']
    profile = utils.get_vertical_diffusion_profile(w_10, depth * -1, diffusion_type)
    # X axis = Concentration axis
    ax2.set_xlabel(r'$K_z$ (m$^2$ s$^{-1}$)', fontsize=ax_label_size)
    ax2.set_xlim((0, 0.05))
    # The actual plotting
    ax2.plot(profile, depth, color=color, linestyle=linestyle, label=label_diffusivity_profile(w_10, diffusion_type))
    # This is if we want to include the diffusion profile gradient too, but this tends to look a lot messier
    if gradient:
        ax2.set_xlim((-0.005, 0.05))
        profile = utils.get_vertical_diffusion_gradient_profile(w_10, depth * -1, diffusion_type)
        ax2.plot(profile, depth, color='red', linestyle='-.',
                 label=label_diffusivity_profile(w_10, diffusion_type))
    return ax2


def label_diffusivity_profile(w_10, diffusion_type):
    """
    Label for the diffusion profile
    :param w_10:
    :param diffusion_type:
    :return:
    """
    if diffusion_type == 'SWB':
        return r'SWB, u$_{10}$' + ' = {:.2f}'.format(w_10) + ' m s$^{-1}$'
    elif diffusion_type == 'KPP':
        return r'KPP, u$_{10}$ '+'= {:.2f}'.format(w_10) + 'm s$^{-1}$, MLD = ' + '{} m'.format(settings.MLD)


def base_figure(fig_size, ax_range, y_label, x_label, ax_label_size, shape=(1, 1), plot_num=1, all_x_labels=False,
                legend_axis=False):
    """
    Function creating the base figure that we use as a foundation for almost all figures
    :param fig_size: size of the figure
    :param ax_range: the limits of the x and y axes
    :param y_label: the y label
    :param x_label: the x label
    :param ax_label_size: the fontsize of the axis labels
    :param shape: the shape of the array (rows, columns)
    :param plot_num: how many subplots we want to create (e.g. in case we have a 2x3 figure but only
    :param all_x_labels: if True, all subplots in the bottom row of teh figure will have x labels, otherwise just the
                         middle one
    :param legend_axis: if true, we add an additional column in which we can add the legend (in case it is too big to
                        fit within a subplot)
    :return:
    """
    # Loading the axis limits
    xmax, xmin, ymax, ymin = ax_range
    # Creating the figure
    fig = plt.figure(figsize=fig_size)
    # If we have just a single subplot in the figure
    if shape == (1, 1):
        if legend_axis:
            grid = GridSpec(nrows=shape[0], ncols=shape[1] + 1, figure=fig)
        else:
            grid = GridSpec(nrows=shape[0], ncols=shape[1], figure=fig)
        ax_sub = fig.add_subplot(grid[0, 0])
        # Y axis = Depth axis
        ax_sub.set_ylabel(y_label, fontsize=ax_label_size)
        ax_sub.set_ylim((ymin, ymax))
        ax_sub.tick_params(axis='both', labelsize=ax_label_size)
        # X axis = Concentration axis
        ax_sub.set_xlabel(x_label, fontsize=ax_label_size)
        ax_sub.set_xscale('log')
        ax_sub.set_xlim((xmin, xmax))
        if not legend_axis:
            return ax_sub
        else:
            ax = [ax_sub]
            ax_legend = fig.add_subplot(grid[0, -1])
            ax_legend.set_axis_off()
            ax.append(ax_legend)
            return ax
    # if we have a more complex figure than just one subplot
    else:
        if legend_axis:
            grid = GridSpec(nrows=shape[0], ncols=shape[1] + 1, figure=fig)
        else:
            grid = GridSpec(nrows=shape[0], ncols=shape[1], figure=fig)
        ax, fig_index = [], 1
        for row in range(shape[0]):
            for column in range(shape[1]):
                ax_sub = fig.add_subplot(grid[row, column])
                ax_sub.set_ylim((ymin, ymax))
                ax_sub.tick_params(axis='both', labelsize=ax_label_size)
                # Only add y labels if we are in the first column
                if column == 0:
                    ax_sub.set_ylabel(y_label, fontsize=ax_label_size)
                else:
                    ax_sub.tick_params(labelleft=False)
                # Only add x labels if we are in the bottom row, and only to the middle one unless all_x_labels == True
                ax_sub.set_xscale('log')
                ax_sub.set_xlim((xmin, xmax))
                if row == (shape[0] - 1):
                    if not all_x_labels and column % 2 is 1:
                        ax_sub.set_xlabel(x_label, fontsize=ax_label_size)
                    elif all_x_labels:
                        ax_sub.set_xlabel(x_label, fontsize=ax_label_size)
                    else:
                        ax_sub.tick_params(labelbottom=False)
                else:
                    ax_sub.tick_params(labelbottom=False)
                # Add the axis to the list
                ax.append(ax_sub)
                # Continuing on to the next plot
                fig_index += 1
                if fig_index > plot_num:
                    break
        # If we want a seperate axis just for adding a legend, then here we create that axis (without any axis lines_
        if legend_axis:
            for rows in range(shape[0]):
                ax_legend = fig.add_subplot(grid[rows, -1])
                ax_legend.set_axis_off()
                ax.append(ax_legend)
        return tuple(ax)


def label_SWB(selection, parameters):
    """ Basic label for a plot showing SWB diffusion results """
    w_10, w_rise = parameters
    w_rise = np.abs(w_rise)
    if selection is 'w_rise':
        return r'SWB, w$_{rise}$ '+'= {} m s'.format(w_rise) + r'$^{-1}$'
    elif selection is 'w_10':
        return r'SWB, u$_{10}$ '+'= {} m s'.format(w_10) + r'$^{-1}$'


def label_KPP(selection, parameters, mld=settings.MLD):
    """ Basic label for a plot showing KPP diffusion results """
    w_10, w_rise = parameters
    w_rise = np.abs(w_rise)
    if selection is 'w_rise':
        return r'KPP, w$_{rise}$ '+'= {}'.format(w_rise) + 'm s$^{-1}$, MLD = ' + '{:.1f} m'.format(mld)
    elif selection is 'w_10':
        return r'KPP, u$_{10}$ '+'= {}'.format(w_10) + 'm s$^{-1}$, MLD = ' + '{:.1f} m'.format(mld)


def label_alpha_comparison(boundary, alpha):
    """ Basic label for when we compare an M-0 with M-1 runs with varying values of alpha"""
    if boundary is 'Reflect':
        return 'M-0'
    else:
        return r'M-1, $\alpha$ = ' + '{}'.format(alpha)


def get_axes_range(close_up, norm_depth, delta_y=0.5):
    """
    A function to specify the axes ranges
    :param close_up: if not None, then it indicates the x-axis limits as (max, min)
    :param norm_depth: if True, then we plot the depths normalized by the MLD
    :param delta_y: this slightly extends the limits on the y axis so that markers right on the border are not cut off
    :return: (xmax, xmin, ymax, ymin)
    """
    # The x-axis shows the relative concentration
    xmax, xmin = 1e0, 1e-4
    if norm_depth:
        delta_y = 0.1
    # If we have a close_up plot
    if close_up is not None:
        ymax, ymin = close_up
        ymax += delta_y
    else:
        # If we have normalized depths
        if norm_depth:
            ymax, ymin = 0 + delta_y, -1
        else:
            ymax, ymin = 0 + delta_y, -1 * settings.max_depth - delta_y
    return xmax, xmin, ymax, ymin


def determine_linestyle(boundary, boundary_list, kpp, swb, diffusion_type):
    """
    Given the boundary and diffusion of a run we are plotting, specifying which linestyle to use
    :return:
    """
    if len(boundary_list) is 1:
        if swb and kpp:
            if diffusion_type is 'SWB':
                return '-'
            elif diffusion_type is 'KPP':
                return '--'
        elif swb:
            return '-'
        elif kpp:
            return '--'
    # For when we are trying to plot multiple boundary conditions in one plot
    else:
        line_style = {'Reflect': '-', 'Reflect_Markov': '--'}
        return line_style[boundary]


def boolean_diff_type(diffusion_type):
    """
    A basic function, that indicates which diffusion types are to be included in a plot
    Note: by all, we just mean KPP and SWB, since the artificial profile was just a test and not generally used for any
    results
    :param diffusion_type:
    :return:
    """
    if diffusion_type is 'SWB':
        swb, kpp, artificial = True, False, False
    elif diffusion_type is 'KPP':
        swb, kpp, artificial = False, True, False
    elif diffusion_type is 'all':
        swb, kpp, artificial = True, True, False
    elif diffusion_type is 'artificial':
        swb, kpp, artificial = False, False, True
    return swb, kpp, artificial


def get_concentration_list(w_10_list, w_rise_list, selection, single_select, diffusion_type, alpha_list, output_step=-1,
                           all_timesteps=False, boundary='Ceiling', mld=settings.MLD, dt=settings.dt_int.seconds):
    """
    A function that returns a dictionary, which in turns contains lists containing all concentration profiles for the
    runs specified in w_10_list, w_rise_list and alpha_list, for the given boundary condition and diffusion type
    Note: This is could be written much more efficiently now (completely with dictionaries), but given how all other
    code is written to account for this format I won't be rewriting it at this time
    :param w_10_list: list containing all required w_10 values
    :param w_rise_list: list containing all required w_rise values
    :param selection: are we plotting w_10 or w_r values?
    :param single_select: the index that we take of the parameter list that we aren't plotting (so if selection = w_10,
                          then single_select would be the index of the value of w_rise_list that we use to select the
                          data)
    :param diffusion_type: SWB or Kukulka
    :param alpha_list: either a list of alpha values or just a single value (only relevant for M-1 simulations)
    :param output_step: which output step do we want to return? default is -1, so the last output of the simulation
    :param all_timesteps: if True, return all the concentrations for all timesteps and not just the one specified in
                          output_step
    :param boundary: which boundary condition, and M-0 or M-1
    :param mld: the mixing layer depth, with teh default just being taken from the settings.py file
    :param dt: the integration timestep, with the default taken from settings.py
    :return: dictionary containing lists with the relevant outputs
    """
    output_dic = {'concentration_list': [],         # List containing the concentration profiles
                  'parameter_concentrations': []}   # List containing tuples with the w-10 and w_rise values for each
                                                    # profile in 'concentration_list'

    # Using the selection parameter to adapt w_rise_list or w_10_list to contain just one value. selection == 'all' will
    # just return all simulations for a particular diffusion_type and boundary
    if selection == 'w_10':
        w_rise_list = [w_rise_list[single_select]]
    elif selection == 'w_rise':
        w_10_list = [w_10_list[single_select]]
    # Specifying the alpha values
    if type(alpha_list) is list:
        alpha = alpha_list[0]
    else:
        alpha = alpha_list

    # Looping through the w_rise and w_10 values in w_rise_list and w_10_list
    for w_rise in w_rise_list:
        for w_10 in w_10_list:
            # Loading the dictionary containing the concentrations
            input_dir = utils.load_obj(utils.get_concentration_output_name(w_10, w_rise, diffusion_type, boundary,
                                                                           mld=mld, alpha=alpha, dt=dt))
            # Selecting the timeslice of interest
            if output_step == -1:
                concentration = [input_dir['last_time_slice']]
            else:
                concentration = [output_step]
            if all_timesteps:
                concentration = range(input_dir['last_time_slice'] + 1)
            # Appending the concentration profiles to 'concentration_list', where the concentrations are now normalized
            # by the total number of particles in the simulation
            for c in concentration:
                output_dic['concentration_list'].append(input_dir[c] / input_dir[c].sum())
            # Get the depth bins of the concentrations
            output_dic['depth_bins'] = input_dir['bin_edges'][:-1] * -1
            # Saving the physical parameters of the concentration distribution we have just loaded
            output_dic['parameter_concentrations'].append((w_10, w_rise))
    return output_dic


def add_observations(ax, sources=None, wind_range=None, norm_depth=False, alpha=0.5):
    """
    Function that adds field data to the axis object ax
    :param ax: the axis object to which we add the field data
    :param sources: unless we specify a specific field data source (so source is not None), we just add all field data
    :param wind_range: specifying the particular wind range for which want to add field data, otherwise we just add all
                       the field data
    :param norm_depth: if True, normalize the field data by the MLD
    :param alpha: setting the opacity of the markers
    :return:
    """
    # If no individual sources as specified, just assume I want to include all of them
    if sources is None:
        sources = ['Kooi', 'Pieper', 'Zettler', 'Kukulka', 'Egger']
    data_markers = ['v', 'P', '*', 'o', 's']

    data_labels = {'Kooi': 'Kooi et al. (2016)', 'Pieper': 'Pieper et al. (2020)',
                   'Zettler': 'Amaral-Zettler (Unpublished data)',
                   'Kukulka': 'Kukulka et al. (2012)', 'Egger': 'Egger et al. (2020)'}

    # Adding the actual plotting
    for count, source in enumerate(sources):
        # Load the correct data dictionary
        data_dict = utils.load_obj(utils.get_data_output_name(source))
        source_conc = data_dict['concentration']
        # Normalized depths vs un-normalized
        if norm_depth:
            source_depth = data_dict['depth_norm']
        else:
            source_depth = data_dict['depth']
        # Selecting the data within the wind_range
        if wind_range is not None:
            w_min, w_max = wind_range
            wind = data_dict['wind_speed']
            wind_select = (w_min < wind) & (wind <= w_max)
            if sum(wind_select) is 0:
                break
            else:
                source_depth = source_depth[wind_select]
                source_conc = source_conc[wind_select]

        # Making the actual plot
        ax.scatter(source_conc, -1 * source_depth, marker=data_markers[count], c='black',
                   alpha=alpha, label=data_labels[source])
    # Getting the lines and labels for the data, so we can add these later to the legend of the figure
    lines, labels = ax.get_legend_handles_labels()
    return lines, labels


def discrete_color_from_cmap(index, subdivisions, cmap='viridis_r'):
    """
    For a given colormap, we split it up into a number of sections equal to subdivisions, and then return the RGB values
    for the section corresponding to index
    :param index: which section we want the color for
    :param subdivisions: into how many sections do we want to split the colormap
    :param cmap: which colormap to use
    :return: the RGB values for the given index
    """
    cmap_steps = plt.cm.get_cmap(cmap, subdivisions)
    return cmap_steps(index)


def rise_velocity_selector(size_class, w_rise_list):
    """
    This removes the smallest rise velocity from w_rise_list
    :param size_class:
    :param w_rise_list:
    :return:
    """
    if size_class is 'large':
        if -0.0003 in w_rise_list:
            w_rise_list.remove(-0.0003)
    return w_rise_list


def return_color(index):
    """
    this returns the color from the list colors
    :param index:
    :return:
    """
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray',
              'tab:olive', 'tab:cyan']
    num = len(colors)
    return colors[index % num]