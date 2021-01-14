# Data Processing functions
## Luke Colosi | lcolosi@ucsd.edu | December 22nd, 2020

########### Shifting Longitude Orientation Function ###########
def shift_grid(data, lon, dlon):
    """
    shift_grid(data, lon, dlon)

        Function to shift global data horizontally to the right in order to center the geographic grid around a
        desired longitudinal location.

        Parameters
        ----------
        data : Numpy 2d masked array that will be shifted.
        lon : Numpy array longitude coordinates.
        dlon : Integer or float value corresponding to how many degrees of longitude are to be shifted. Positive dlon
               corresponds to a shift right and a negative dlon corresponds to a dhift left.

        Returns
        -------
        data_shift : numpy 2d array of shifted data
        lon_shift : numpy array of shifted longitude

        Libraries necessary to run function
        -----------------------------------
        Numpy : import numpy as np
    """

    # import libraries:
    import numpy as np

    # initialize variables:
    nlon = len(lon)
    lon_min = np.min(lon)
    lon_max = np.max(lon)

    # Set the resolution of longitude:
    dl = lon[1] - lon[0]

    # Determine how much to shift over by:
    nshift = int(round(dlon / dl))

    # set indices
    indxs_1, indxs_2 = np.arange(nshift, nlon, 1), np.arange(0, nshift, 1)
    indxs = np.concatenate((indxs_1, indxs_2), axis=0)

    # Apply to data and longitude:
    lon_shift = np.arange(lon_min + dlon, lon_max + dlon + 1, dl)
    data_shift = data[:, indxs]

    return data_shift, lon_shift


########### Importing processed Data Function ###########
def import_data(data_name, path):

    """
    Function to import IFREMER binned, CCMP2 deresolved, and WW3 deresolved data.

        import_data(data_name, path)

        Parameters
        ----------
        data_name : Specifies the data which will be imported as a string. Options include: data_name = 'IFREMER_swh',
                    'CCMP2_wsp', 'WW3_swh', 'WW3_wsp'.
        path : Path, as a string, to the directory where raw data is stored. Path can be relative or absolute.
            ex: path = '../data/'

        Returns
        -------
        data : Numpy 3d masked array of specified data with dimensions 8400, 133, 360.
        time : List of daily time steps with 8400 elements.
        lat : Numpy 1d array of degrees of Latitude
        lon : Numpy 1d array of degrees of Longitude

        Libraries necessary to run function
        -----------------------------------
        import numpy as np
        from netCDF4 import Dataset, num2date
        import glob
    """

    # Import libraries:
    import numpy as np
    from netCDF4 import Dataset, num2date
    import glob

    # set dimensions:
    nt, nlat, nlon = 8400, 133, 360

    # Set file and variable names
    if data_name == "IFREMER_swh":
        filenames = sorted(glob.glob(path + "IFREMER_binned_alt_swh_*.nc"))
        variable = "swh"
    elif data_name == "CCMP2_wsp":
        filenames = sorted(glob.glob(path + "CCMP2_deresolved_wsp_*.nc"))
        variable = "wsp"
    elif data_name == "WW3_swh":
        filenames = sorted(glob.glob(path + "ww3_deresolved_hs_*.nc"))
        variable = "hs"
    elif data_name == "WW3_wsp":
        filenames = sorted(glob.glob(path + "ww3_deresolved_wsp_*.nc"))
        variable = "wsp"

    # Call data and concatinate:

    # Initialize variables:
    data = np.ma.masked_all([nt, nlat, nlon])
    time = []
    year_c = np.array(
        [
            365,
            365,
            365,
            366,
            365,
            365,
            365,
            366,
            365,
            365,
            365,
            366,
            365,
            365,
            365,
            366,
            365,
            365,
            365,
            366,
            365,
            365,
            365,
        ]
    )
    iday, yc = 0, 0

    # Loop through filenames:
    for f in filenames:

        # set nc variable:
        nc = Dataset(f, "r")

        # call latitude and longitude:
        if yc == 0:
            lon = nc.variables["lon"][:]
            lat = nc.variables["lat"][:]

        # call data and time
        idata = nc.variables[variable][:]
        itime = num2date(nc.variables["time"][:], nc.variables["time"].units)

        # place the hs and time data into the 3D arrays
        data[iday : iday + year_c[yc], :, :] = idata
        time.extend(itime)

        # year counters:
        iday += year_c[yc]
        yc += 1

    return data, time, lat, lon
