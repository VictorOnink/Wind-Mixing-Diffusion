import utils
import settings as SET
import numpy as np
import pandas as pd
from seabird.cnv import fCNV


def data_standardization():
    """
    Running all the data standardization functions. Each standardization function contains a check if the standardized
    file hasn't already been created, and so we only have to rerun this if we lose one of the standardized files for
    some reason
    """
    standardization_kukulka()
    standardization_kooi()
    standardization_Pieper()
    standardization_Zettler()


def standardization_kukulka():
    """
    Data from Kukulka et al., 2012, shared by Tobias Kukulka. The data is for 47 different stations/depths, and is in
    the format of:
    Station/ Depth/ Tow length/ Tow volume/ Plastic/ Volumetric Conc/ Avg. Wind Speed (kts)/ Mix Layer depth (m)
    """
    prefix = 'Kukulka'
    file_name = utils.get_data_output_name(prefix)
    if not utils._check_file_exist(file_name + '.pkl'):
        data = np.genfromtxt(SET.data_dir + 'atlantic_prf.dat')
        # Convert wind data from knots to m/s
        data[:, -2] *= 0.514

        # Normalise all the data, such that for each station I have the concentration as a fraction of the max
        # concentration recorded at that station
        station_numbers = np.unique(data[:, 0])
        for station in station_numbers:
            data[data[:, 0] == station, 5] /= np.max(data[data[:, 0] == station, 5])

        # Get a normalised depth array, where the depth is a fraction of the mixing layer depth
        depth_norm = data[:, 1] / data[:, -1]

        # Create a dictionary containing the normalized concentration, the depth, the normalized depth and wind speed
        output_dic = {'concentration': data[:, 5], 'depth': data[:, 1], 'depth_norm': depth_norm,
                      'wind_speed': data[:, -2]}

        # Pickling the array
        utils.save_obj(filename=file_name, object=output_dic)


def standardization_kooi():
    """
    Data from Kooi et al. (2016), located on figshare at https://dx.doi.org/10.6084/m9.figshare.3427862.v1
    In short, a multinet was used where all the nets are at 0.5m intervals. The datafile contains extensive data on the
    meteographic/oceanographic conditions, and the type of debris inside the net. There are also estimates of rise
    velocities of some of the fragments/line items.

    I don't have estimates of the mixing layer depth, but I do have CTD data for each station, and I could see about
    what exact approach Kukulka used to have the same measure the mixing layer depth (depends on the station how far
    down the CTD data was collected
    """
    prefix = 'Kooi'
    file_name = utils.get_data_output_name(prefix)
    if not utils._check_file_exist(file_name + '.pkl'):
        data_trawl = pd.read_excel(SET.data_dir + 'Data_KooiEtAl.xlsx', sheet_name='trawls')
        data_plastic = pd.read_excel(SET.data_dir + 'Data_KooiEtAl.xlsx', sheet_name='nets')


        # Determine the mixed layer depth
        MLD = determine_MLD(prefix)

        # Convert wind data from knots to m/s
        wind_data = 0.514 * data_trawl['WindSpeedKnots']

        # The net number corresponds to the depth at which a measurement was done
        data_plastic['Net'] = (data_plastic['Net'] - 1) * 0.5
        depth_levels = np.unique(data_plastic['Net'])

        # The station numbers (set such that the first station has index 0, also for the station column in the data
        # sheet)
        station_numbers = data_trawl['Station'] - 1
        data_plastic['Station'] -= 1

        # Creating an array that will contain all the concentrations for each station and depth
        concentration = pd.DataFrame(columns=station_numbers.values, index=depth_levels).fillna(0.0)

        # Now, cycling through all the stations to get the concentration as counts per depth level (essentially
        # concentrations)
        for station in station_numbers:
            for depth in depth_levels:
                concentration[station][depth] = data_plastic['NumberParticles'][(data_plastic['Station'] == station) &
                                                                                (data_plastic['Net'] == depth)].sum()

        # Normalizing everything by the max concentration
        concentration = concentration.apply(lambda x: x/x.max(), axis=0).values

        # Getting the wind_data and depth_levels into the same shape as the concentration array
        wind_data = pd.concat([wind_data]*depth_levels.shape[0], axis=1).transpose().values
        depth_levels = np.array([depth_levels, ]*station_numbers.shape[0]).transpose()
        MLD = np.array([MLD]*depth_levels.shape[0]).reshape(depth_levels.shape)

        # Normalising the depth array
        depth_norm = np.divide(depth_levels, MLD)

        # Create a dictionary containing the normalized concentration, the depth, the normalized depth and wind speed
        output_dic = {'concentration': concentration.flatten(), 'depth': depth_levels.flatten(),
                      'depth_norm': depth_norm.flatten(), 'wind_speed': wind_data.flatten()}

        # Pickling the array
        utils.save_obj(filename=file_name, object=output_dic)


def standardization_Pieper():
    """
    Data from Pieper et al. (2020), collected during the PE442 cruise from the Azores to Sicily. The data file was
    shared by Catharina Pieper.

    Unlike the previous data, these samples were not collected with manta-trawls, but instead using niskin bottles on
    the CTD carriage. Each sample is the number of microplastics in a liter of sample, where most of these were fibers,
    and as such they are good examples of particles with very low rising velocities.

    Sample concentrations and depths are from the data shared by Catharina, while the wind speed comes from the Casino
    data file from the R.V. Pelagia. CTD data is currently unavailable
    """
    prefix = 'Pieper'
    file_name = utils.get_data_output_name(prefix)
    if not utils._check_file_exist(file_name + '.pkl'):
        data_bottle = pd.read_excel(SET.data_dir + '2020_PE442_MPs.xlsx')

        # Get the station indicators, sample concentrations, sample depths, and sample types
        station_sample = data_bottle['Station Number']
        station_numbers = np.unique(station_sample)
        counts = data_bottle['MP/Liter']
        depths = data_bottle.Depth
        type = data_bottle['Filter Type']
        type_list = np.unique(type)

        # Creating a dataframe in which all the concentrations will go
        concentration = pd.DataFrame(columns=station_numbers, index=type_list).fillna(0.0)
        depth_dataframe = pd.DataFrame(columns=station_numbers, index=type_list).fillna(0.0)

        # Go through all the stations, and get the concentration at each depth, where we subtract the blanks and air
        # contamination counts
        for station in station_numbers:
            control = counts[(station == station_sample) & (type == 'Blank')].values + \
                      counts[(station == station_sample) & (type == 'AirContamination')].values
            # Go through all the different types of samples that were collected
            for sample in type_list:
                # There were always two replicas at each depth for each station, and I want the average
                replicas = counts[(station == station_sample) & (type == sample)].values
                # If a particular sample type was present at a station, then compute the mean concentration,
                # and subtracting the number of particles in the control
                if len(replicas) > 0:
                    concentration[station][sample] = max(0, np.mean(replicas) - control)
                    depth_dataframe[station][sample] = np.mean(depths[(station == station_sample) & (type == sample)].values)

        # Normalizing everything by the max concentration
        concentration = concentration.apply(lambda x: x / x.max(), axis=0).fillna(0.0)

        # Wind speed from the PE442 Casino file
        wind_data = pd.DataFrame(casino_wind(device='CTD with samples', cruise='PE442'))
        wind_data = pd.concat([wind_data] * type_list.shape[0], axis=1).transpose().values


        # Getting the mixing layer depth from the CTD data
        MLD = determine_MLD(prefix=prefix, station_numbers=station_numbers).values
        MLD = np.array([MLD, ]*type_list.shape[0]).reshape(type_list.shape[0], station_numbers.shape[0])

        # Normalising the depth according to MLD depth
        depth_norm = np.divide(depth_dataframe.values, MLD)

        # Saving everything into a dictionary
        output_dic = {'concentration': concentration.values.flatten(), 'depth': depth_dataframe.values.flatten(),
                      'depth_norm': depth_norm.flatten(), 'wind_speed': wind_data.flatten()}

        # Pickling the array
        utils.save_obj(filename=file_name, object=output_dic)


def standardization_Zettler():
    """
    Data collected during the PE448 South Atlantic cruise in January 2019. The data has been shared by Erik Zettler and
    is currently not published yet. Samples reflect sampled microplastic concentrations using a multinet, and samples
    using manta trawl data at the ocean surface.
    """
    prefix = 'Zettler'
    file_name = utils.get_data_output_name(prefix)
    if not utils._check_file_exist(file_name + '.pkl'):
        data_multi = pd.read_excel(SET.data_dir + 'PE448_multinet_data.xlsx')
        data_surf = pd.read_excel(SET.data_dir + 'Sample Log-PE448b-20190121.xlsx', sheet_name='MT')

        # Get the sample depths, counts, and volumes for the multi-net, and then the concentration (counts/volume)
        depth = data_multi.depth.dropna().reset_index(drop=True)
        counts_multi = data_multi['#pieces plastic']
        volume_multi = data_multi['volume m^3 (per flow meter)']
        concentration_multi = pd.DataFrame(np.divide(counts_multi.values,
                                                     volume_multi.values)).dropna().reset_index(drop=True)

        # Get the concentration for the surface trawls at the same point as the multi-net trawls
        counts_surf = data_surf.loc[[4, 13, 17], ['Total # of plastic pieces in tow']]
        volume_surf = data_surf.loc[[4, 13, 17], ['volume filtered (m^3; net mouth is 15cm high)']]
        concentration_surf = pd.DataFrame(np.divide(counts_surf.values, volume_surf.values))

        # Get a nice array with three columns corresponding to the three 'super station' sites
        concentrations = pd.DataFrame(columns=[0, 1, 2], index=[0, 1, 2, 3, 4, 5]).fillna(0.0)
        depths = pd.DataFrame(columns=[0, 1, 2], index=[0, 1, 2, 3, 4, 5]).fillna(0.0)
        for station in [0, 1, 2]:
            for rows in [0, 1, 2, 3, 4, 5]:
                if rows is 0:
                    concentrations[station][rows] = concentration_surf[0][station]
                else:
                    concentrations[station][rows] = concentration_multi[0][(rows - 1) + station * 5]
                    depths[station][rows] = depth[(rows - 1) + station * 5]

        # I don't have CTD data at the moment for these stations, so we'll just have an empty array for this for now
        MLD = pd.DataFrame(columns=[0, 1, 2], index=[0, 1, 2, 3, 4, 5])

        # Normalizing the concentrations and depths
        depth_norm = np.divide(depths.values, MLD.values)
        concentrations = concentrations.apply(lambda x: x / x.max(), axis=0).fillna(0.0)

        # Getting the wind data
        wind_data = pd.DataFrame(casino_wind(device='MultiNet', cruise='PE448'))
        wind_data = pd.concat([wind_data] * depths.shape[0], axis=1).transpose().values

        # Getting just the concentrations that are greater than 0
        greater_than = concentrations.values.flatten() > 0

        # Saving everything into a dictionary
        output_dic = {'concentration': concentrations.values.flatten()[greater_than],
                      'depth': depths.values.flatten()[greater_than],
                      'depth_norm': depth_norm.flatten()[greater_than],
                      'wind_speed': wind_data.flatten()[greater_than]}

        # Pickling the array
        utils.save_obj(filename=file_name, object=output_dic)



def determine_MLD(prefix: str, station_numbers=None):
    """
    Determine the mixing layer depth according to de Boyer Montegut et al. (2004), which is the MLD calculation approach
    used by Kukulka et al. (2004).

    MLD = depth at which there is a 0.2 degree temperature difference relative to the temperature at 10m depth
    """
    z_ref = 10 # reference depth in meters
    dif_ref = 0.2 # temperature difference relative to reference depth (degrees celcius)

    if prefix is 'Kooi':
        data_ctd = pd.read_excel(SET.data_dir + 'Data_KooiEtAl.xlsx', sheet_name='CTD')
        data_ctd.station -= 1
        station_numbers = np.sort(np.append(data_ctd.station.unique(), 24))
        # Array in which to store the determined MLD values
        MLD = np.zeros((1, station_numbers.shape[0]))
        for station in station_numbers:
            if station == 24:
                # CTD data for station 24 is missing from the datasheet
                MLD[0, station] = np.nan
            else:
                # Load the depth and temperature data for that particular station, where we also include a check that all
                # of the depth files are sorted correctly
                depth = data_ctd['Depth'][data_ctd.station == station].values
                depth_sort = depth.argsort()
                depth = depth[depth_sort]
                temp = data_ctd['Temperature'][data_ctd.station == station].values[depth_sort]

                # Determine the index that corresponds to a depth of 10m
                ind_10 = utils.find_nearest_index(depth=depth, z_ref=z_ref)
                temp_10 = temp[ind_10]

                # Determine the depth at which the temperature difference is equal to dif_ref with respect to z_ref
                depth, temp = depth[ind_10:], temp[ind_10:]
                MLD[0, station] = depth[np.where(np.abs(temp - temp_10) > dif_ref)[0][0]]

    if prefix is 'Pieper':
        MLD = pd.DataFrame(columns=station_numbers, index=[0]).fillna(0.0)
        # Check if there is a CTD file for the station in question
        for station in station_numbers:
            file_name = SET.data_dir + 'CTD_PE442/PE442_' + station + 'avg.cnv'
            if utils._check_file_exist(file_name):
                # Load the depth and temperature data for the particular station
                temperature = fCNV(file_name)['TEMP']
                depth = fCNV(file_name)['DEPTH']

                # Determine the index that corresponds to a depth of 10m
                ind_10 = utils.find_nearest_index(depth=depth, z_ref=z_ref)
                temp_10 = temperature[ind_10]

                # Determine the depth at which the temperature difference is equal to dif_ref with respect to z_ref
                depth, temp = depth[ind_10:], temperature[ind_10:]
                MLD[station] = depth[np.where(np.abs(temp - temp_10) > dif_ref)[0][0]]
            else:
                MLD[station] = np.nan
    return MLD


def casino_wind(device: str, cruise:str):
    if cruise is 'PE442':
        data = pd.read_csv(SET.data_dir + 'casino_{}.csv'.format(cruise), delimiter='\t')
    elif cruise is 'PE448':
        data = pd.read_excel(SET.data_dir + 'casino_{}.xlsx'.format(cruise))

    # Now, we want the wind data for the points in time when the device in question starts it's deployment
    wind_data = data.loc[(data['Device name'] == device) &
                         (data['Action name'] == 'Begin')]['PE_WEATHER_01_weather_trueairspeed'].reset_index(drop=True)

    if cruise is 'PE448':
        # There is an issue with the casino file. Based on my own recollection and examining a number of data points
        # that appeared to escape corruption, these corrections should now yield the true wind speed
        wind_data = pd.DataFrame(np.divide(wind_data.values, [1000, 1, 10000]))
    # For PE442, there was one occasion where the CTD malfunctioned, so we have three measurements that actually
    # correspond to just one station, namely points 4 and 5 are actually for station 5
    if cruise is 'PE442':
        wind_data = wind_data.drop(labels=[4, 5])
    return wind_data.values
