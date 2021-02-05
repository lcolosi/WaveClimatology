# Functions for saving netCDF files for data repository
## Luke Colosi | lcolosi@ucsd.edu | October 18th, 2020

########### Global Attibute Function ###########
def add_global_atrributes(nc, summary):

    """
    add_global_atrributes(nc, attributes)

        Adds global attributes to an nc file from an ordered dictionary of attributes.

        Parameters
        ----------
        nc : The netCDF variable intitalized with nc = Dataset(filename, 'w', format='NETCDF4') which is used to
             write global attributes, variables, and variable attributes to the netCDF file.
        summary : A string containing the summary of what is contained in the netCDF file.

        Returns
        -------

        Libraries necessary to run function
        -----------------------------------
        import time
        from collections import OrderedDict


    """

    # Import libraries:
    import time
    from collections import OrderedDict

    attributes = OrderedDict()
    attributes[
        "title"
    ] = "Data from: The Seasonal Cycle of Significant Wave Height in the Ocean:  Local vs Remote Forcing"
    attributes["summary"] = summary
    attributes[
        "keywords"
    ] = "wave and wind climate, Surface gravity waves, Satellite altimetry, Annual cycle, Probability of swell"
    attributes[
        "references"
    ] = "Colosi, L. V., Villas Boas, A. B., and Gille, S. T. (2020). The Seasonal Cycle of Significant Wave Height in the Ocean: Local vs Remote Forcing. Journal of Geophysical Research: Oceans, 1(1), 1st ser., 1-2. doi:"
    attributes[
        "acknowledgement"
    ] = "This work was supported by the NASA SWOT (awards NNX16AH67G and 80NSSC20K1136) and Ocean Vector Winds Science Teams (award 80NSSC19K0059), by a NASA Earth and Space Science Fellowship awarded to Ana Villas Boas, and by the Philip and Elizabeth Hiestand Scholars program."
    attributes["date_created"] = time.ctime()
    attributes["creator_name"] = "Luke Vincent Colosi"
    attributes["creator_email"] = "lcolosi@ucsd.edu"
    attributes[
        "creator_institution"
    ] = "University of California - San Diego; Scripps Institution of Oceanography"
    attributes["contributor_name"] = "Bia Villas Boas"
    attributes["contributor_role"] = "Methodology"
    attributes["standard_name_vocabulary"] = "CF Standard Name Table v72"
    attributes["doi"] = ""

    # loop through global attributes and write to nc file
    for att in attributes.keys():
        nc.setncattr(att, attributes[att])


########### IFREMER Binned Along Track SWH Function ###########
def save_netcdf_binned_swh(swh, nobs, lon, lat, time, output, summary):

    """
    save_netcdf_binned_swh(swh, nobs lon, lat, time, output, summary)

        Function to save corrected binned satellite altimeter swh and wsp with longitude, latitude, and time
        variables into a netCDF file in the current directory or directory specified in the output variable.

        Parameters
        ----------
        swh : numpy 3D masked array of IFREMER correct swh
        nobs : numpy 3D masked array of number of observations averaged in bins
        lon : numpy array column vector of longitude coordinates
        lat : numpy array column vector of longitude coordinates
        time : numpy array of date2num values of time for days
        output: filename (path to file and file's name)
            ex: output = '/zdata/downloads/colosi_data_bk/binned_data/ifremer_p1_daily_data/
                          my_daily_binned_ifremer_data/IFREMER_binned_alt_swh_wsp_93_16.nc'
        summary : A string containing the summary of what is contained in the netCDF file.

        Returns
        -------

        Libraries necessary to run function
        -----------------------------------
        from netCDF4 import Dataset, num2date, date2num
        from datetime import datetime
        import netCDF4

        Important Note
        --------------
        The NetCDF file cannot be saved to a directory that has a file with the same name as the file being saved
        (permission will be denied to write over that file. Therefore, make sure all data files in the directory one
        is saving to have different names).
    """

    # Set path to my functions:
    import sys

    sys.path.append("../tools/")

    # import libraries:
    from netCDF4 import Dataset, num2date, date2num
    from datetime import datetime
    import netCDF4

    # import functions
    from save_netcdf_fields import add_global_atrributes

    # set dimensions for swh, latitude and longitude
    Nt, Ny, Nx = swh.shape

    # Initiate netCDF file in directory specified by output variable
    nc = Dataset(output, "w", format="NETCDF4")

    # set global attributes:
    add_global_atrributes(nc, summary)

    # set time, longitude and latitude dimensions
    time_dim = nc.createDimension("time", Nt)
    lon_dim = nc.createDimension("lon", Nx)
    lat_dim = nc.createDimension("lat", Ny)

    # Initiate dictionary for attribute for all variables:
    input_vars = {}

    # Set attributes for each variable
    ######## time ########
    input_vars["time"] = {}
    input_vars["time"]["calendar"] = "gregorian"
    input_vars["time"]["units"] = "days since 1900-01-01 00:00:00"
    input_vars["time"]["long_name"] = "Time"
    input_vars["time"]["standard_name"] = "time"
    input_vars["time"]["coverage_content_type"] = "coordinate"

    ######## latitude ########
    input_vars["lat"] = {}
    input_vars["lat"]["units"] = "degrees_north"
    input_vars["lat"]["long_name"] = "Latitude"
    input_vars["lat"]["standard_name"] = "latitude"
    input_vars["lat"]["coverage_content_type"] = "coordinate"

    ######## longitude ########
    input_vars["lon"] = {}
    input_vars["lon"]["units"] = "degrees_east"
    input_vars["lon"]["long_name"] = "Longitude"
    input_vars["lon"]["standard_name"] = "longitude"
    input_vars["lon"]["coverage_content_type"] = "coordinate"

    ######## Number of Observations ########
    input_vars["nobs"] = {}
    input_vars["nobs"]["long_name"] = "Number of Observations"
    input_vars["nobs"]["standard_name"] = "number_of_observations"
    input_vars["nobs"]["ancillary_variables"] = "swh"
    input_vars["nobs"]["coverage_content_type"] = "physicalMeasurement"
    input_vars["nobs"]["coordinates"] = "time lat lon"
    input_vars["nobs"][
        "comment"
    ] = "Provides a 3-dimensional array of the number of data points average in 1 by 1 degree bins."

    ######## SWH ########
    input_vars["swh"] = {}
    input_vars["swh"]["units"] = "m"
    input_vars["swh"][
        "long_name"
    ] = "binned corrected altimeter significant wave height"
    input_vars["swh"]["standard_name"] = "sea_surface_wave_significant_height"
    input_vars["swh"]["coverage_content_type"] = "physicalMeasurement"
    input_vars["swh"]["coordinates"] = "time lat lon"

    # Initialize variables:
    vars = {}

    # Create variable, set fill_value to the netCDF4 default_fillvals with 32 or 64-bit floating point:
    for var in input_vars.keys():

        if var == "time":
            vars["time"] = nc.createVariable(
                "time", "<f4", ("time",), fill_value=netCDF4.default_fillvals["f4"]
            )
        elif var == "lon":
            vars["lon"] = nc.createVariable(
                "lon", "<f4", ("lon",), fill_value=netCDF4.default_fillvals["f4"]
            )
        elif var == "lat":
            vars["lat"] = nc.createVariable(
                "lat", "<f4", ("lat",), fill_value=netCDF4.default_fillvals["f4"]
            )
        else:
            vars[var] = nc.createVariable(
                var,
                "<f8",
                ("time", "lat", "lon"),
                fill_value=netCDF4.default_fillvals["f8"],
            )

        # Set variable attributes
        for a in input_vars[var].keys():
            setattr(vars[var], a, input_vars[var][a])

    # Force netCDF file to convert variables to masked arrays when read
    for var_all in vars.keys():
        vars[var_all].set_auto_maskandscale(True)

    # assign data to netCDF variables
    vars["lat"][:] = lat
    vars["lon"][:] = lon
    vars["time"][:] = date2num(time, nc.variables["time"].units)
    vars["nobs"][:] = nobs
    vars["swh"][:] = swh

    # Close netCDF file modifications
    nc.close()


########### CCMP2 and WW3 Deresolved WSP Function ###########
def save_netcdf_deresolved_wsp(wsp, lon, lat, time, output, summary):

    """
    save_netcdf_deresolved_wsp(wsp, lon, lat, time, output, summary)

        Function to save CCMP2 and WW3 deresolved wsp with longitude, latitude, and time
        variables into a netCDF file in the current directory or directory specified in the output variable.

        Parameters
        ----------
        wsp : numpy 3D masked array of WSP
        lon : numpy array column vector of longitude coordinates
        lat : numpy array column vector of longitude coordinates
        time : numpy array of date2num values of time for days
        output: filename (path to file and file's name)
            ex: output = '/zdata/downloads/colosi_data_bk/binned_data/ifremer_p1_daily_data/
                          my_daily_binned_ifremer_data/IFREMER_binned_alt_swh_wsp_93_16.nc'
        summary : A string containing the summary of what is contained in the netCDF file.

        Returns
        -------

        Libraries necessary to run function
        -----------------------------------
        from netCDF4 import Dataset, num2date, date2num
        from datetime import datetime
        import netCDF4

        Important Note
        --------------
        The NetCDF file cannot be saved to a directory that has a file with the same name as the file being saved
        (permission will be denied to write over that file. Therefore, make sure all data files in the directory one
        is saving to have different names).
    """

    # Set path to my functions:
    import sys

    sys.path.append("../tools/")

    # import libraries:
    from netCDF4 import Dataset, num2date, date2num
    from datetime import datetime
    import netCDF4

    # import functions
    from save_netcdf_fields import add_global_atrributes

    # set dimensions for swh, latitude and longitude
    Nt, Ny, Nx = wsp.shape

    # Initiate netCDF file in directory specified by output variable
    nc = Dataset(output, "w", format="NETCDF4")

    # set global attributes:
    add_global_atrributes(nc, summary)

    # set time, longitude and latitude dimensions
    time_dim = nc.createDimension("time", Nt)
    lon_dim = nc.createDimension("lon", Nx)
    lat_dim = nc.createDimension("lat", Ny)

    # Initiate dictionary for attribute for all variables:
    input_vars = {}

    # Set attributes for each variable
    ######## time ########
    input_vars["time"] = {}
    input_vars["time"]["calendar"] = "gregorian"
    input_vars["time"]["units"] = "days since 1900-01-01 00:00:00"
    input_vars["time"]["long_name"] = "Time"
    input_vars["time"]["standard_name"] = "time"
    input_vars["time"]["coverage_content_type"] = "coordinate"

    ######## latitude ########
    input_vars["lat"] = {}
    input_vars["lat"]["units"] = "degrees_north"
    input_vars["lat"]["long_name"] = "Latitude"
    input_vars["lat"]["standard_name"] = "latitude"
    input_vars["lat"]["coverage_content_type"] = "coordinate"

    ######## longitude ########
    input_vars["lon"] = {}
    input_vars["lon"]["units"] = "degrees_east"
    input_vars["lon"]["long_name"] = "Longitude"
    input_vars["lon"]["standard_name"] = "longitude"
    input_vars["lon"]["coverage_content_type"] = "coordinate"

    ######## WSP ########
    input_vars["wsp"] = {}
    input_vars["wsp"]["units"] = "m s-1"
    input_vars["wsp"]["long_name"] = "wind speed"
    input_vars["wsp"]["standard_name"] = "wind_speed"
    input_vars["wsp"]["coverage_content_type"] = "physicalMeasurement"
    input_vars["wsp"]["coordinates"] = "time lat lon"

    # Initialize variables:
    vars = {}

    # Create variable, set fill_value to the netCDF4 default_fillvals with 32 or 64-bit floating point:
    for var in input_vars.keys():

        if var == "time":
            vars["time"] = nc.createVariable(
                "time", "<f4", ("time",), fill_value=netCDF4.default_fillvals["f4"]
            )
        elif var == "lon":
            vars["lon"] = nc.createVariable(
                "lon", "<f4", ("lon",), fill_value=netCDF4.default_fillvals["f4"]
            )
        elif var == "lat":
            vars["lat"] = nc.createVariable(
                "lat", "<f4", ("lat",), fill_value=netCDF4.default_fillvals["f4"]
            )
        else:
            vars[var] = nc.createVariable(
                var,
                "<f8",
                ("time", "lat", "lon"),
                fill_value=netCDF4.default_fillvals["f8"],
            )

        # Set variable attributes
        for a in input_vars[var].keys():
            setattr(vars[var], a, input_vars[var][a])

    # Force netCDF file to convert variables to masked arrays when read
    for var_all in vars.keys():
        vars[var_all].set_auto_maskandscale(True)

    # assign data to netCDF variables
    vars["lat"][:] = lat
    vars["lon"][:] = lon
    vars["time"][:] = date2num(time, nc.variables["time"].units)
    vars["wsp"][:] = wsp

    # Close netCDF file modifications
    nc.close()


########### WW3 Deresolved Hs Function ###########
def save_netcdf_deresolved_hs(hs, lon, lat, time, output, summary):

    """
    save_netcdf_deresolved_hs(hs, lon, lat, time, output, summary)

        Function to save WW3 deresolved hs with longitude, latitude, and time
        variables into a netCDF file in the current directory or directory specified in the output variable.

        Parameters
        ----------
        hs : numpy 3D masked array of WW3 hs
        lon : numpy array column vector of longitude coordinates
        lat : numpy array column vector of longitude coordinates
        time : numpy array of date2num values of time for days
        output: filename (path to file and file's name)
            ex: output = '/zdata/downloads/colosi_data_bk/binned_data/ifremer_p1_daily_data/
                          my_daily_binned_ifremer_data/IFREMER_binned_alt_swh_wsp_93_16.nc'
        summary : A string containing the summary of what is contained in the netCDF file.

        Returns
        -------

        Libraries necessary to run function
        -----------------------------------
        from netCDF4 import Dataset, num2date, date2num
        from datetime import datetime
        import netCDF4

        Important Note
        --------------
        The NetCDF file cannot be saved to a directory that has a file with the same name as the file being saved
        (permission will be denied to write over that file. Therefore, make sure all data files in the directory one
        is saving to have different names).
    """

    # Set path to my functions:
    import sys

    sys.path.append("../tools/")

    # import libraries:
    from netCDF4 import Dataset, num2date, date2num
    from datetime import datetime
    import netCDF4

    # import functions
    from save_netcdf_fields import add_global_atrributes

    # set dimensions for swh, latitude and longitude
    Nt, Ny, Nx = hs.shape

    # Initiate netCDF file in directory specified by output variable
    nc = Dataset(output, "w", format="NETCDF4")

    # set global attributes:
    add_global_atrributes(nc, summary)

    # set time, longitude and latitude dimensions
    time_dim = nc.createDimension("time", Nt)
    lon_dim = nc.createDimension("lon", Nx)
    lat_dim = nc.createDimension("lat", Ny)

    # Initiate dictionary for attribute for all variables:
    input_vars = {}

    # Set attributes for each variable
    ######## time ########
    input_vars["time"] = {}
    input_vars["time"]["calendar"] = "gregorian"
    input_vars["time"]["units"] = "days since 1900-01-01 00:00:00"
    input_vars["time"]["long_name"] = "Time"
    input_vars["time"]["standard_name"] = "time"
    input_vars["time"]["coverage_content_type"] = "coordinate"

    ######## latitude ########
    input_vars["lat"] = {}
    input_vars["lat"]["units"] = "degrees_north"
    input_vars["lat"]["long_name"] = "Latitude"
    input_vars["lat"]["standard_name"] = "latitude"
    input_vars["lat"]["coverage_content_type"] = "coordinate"

    ######## longitude ########
    input_vars["lon"] = {}
    input_vars["lon"]["units"] = "degrees_east"
    input_vars["lon"]["long_name"] = "Longitude"
    input_vars["lon"]["standard_name"] = "longitude"
    input_vars["lon"]["coverage_content_type"] = "coordinate"

    ######## WSP ########
    input_vars["hs"] = {}
    input_vars["hs"]["units"] = "m"
    input_vars["hs"]["long_name"] = "significant wave height"
    input_vars["hs"]["standard_name"] = "sea_surface_wave_significant_height"
    input_vars["hs"]["coverage_content_type"] = "physicalMeasurement"
    input_vars["hs"]["coordinates"] = "time lat lon"

    # Initialize variables:
    vars = {}

    # Create variable, set fill_value to the netCDF4 default_fillvals with 32 or 64-bit floating point:
    for var in input_vars.keys():

        if var == "time":
            vars["time"] = nc.createVariable(
                "time", "<f4", ("time",), fill_value=netCDF4.default_fillvals["f4"]
            )
        elif var == "lon":
            vars["lon"] = nc.createVariable(
                "lon", "<f4", ("lon",), fill_value=netCDF4.default_fillvals["f4"]
            )
        elif var == "lat":
            vars["lat"] = nc.createVariable(
                "lat", "<f4", ("lat",), fill_value=netCDF4.default_fillvals["f4"]
            )
        else:
            vars[var] = nc.createVariable(
                var,
                "<f8",
                ("time", "lat", "lon"),
                fill_value=netCDF4.default_fillvals["f8"],
            )

        # Set variable attributes
        for a in input_vars[var].keys():
            setattr(vars[var], a, input_vars[var][a])

    # Force netCDF file to convert variables to masked arrays when read
    for var_all in vars.keys():
        vars[var_all].set_auto_maskandscale(True)

    # assign data to netCDF variables
    vars["lat"][:] = lat
    vars["lon"][:] = lon
    vars["time"][:] = date2num(time, nc.variables["time"].units)
    vars["hs"][:] = hs

    # Close netCDF file modifications
    nc.close()


########### WW3 Deresolved fp Function ###########
def save_netcdf_deresolved_fp(fp, lon, lat, time, output, summary):

    """
    save_netcdf_deresolved_fp(fp, lon, lat, time, output, summary)

        Function to save WW3 deresolved fp with longitude, latitude, and time
        variables into a netCDF file in the current directory or directory specified in the output variable.

        Parameters
        ----------
        fp : numpy 3D masked array of WW3 fp
        lon : numpy array column vector of longitude coordinates
        lat : numpy array column vector of longitude coordinates
        time : numpy array of date2num values of time for days
        output: filename (path to file and file's name)
            ex: output = '/zdata/downloads/colosi_data_bk/binned_data/ifremer_p1_daily_data/
                          my_daily_binned_ifremer_data/IFREMER_binned_alt_swh_wsp_93_16.nc'
        summary : A string containing the summary of what is contained in the netCDF file.

        Returns
        -------

        Libraries necessary to run function
        -----------------------------------
        from netCDF4 import Dataset, num2date, date2num
        from datetime import datetime
        import netCDF4

        Important Note
        --------------
        The NetCDF file cannot be saved to a directory that has a file with the same name as the file being saved
        (permission will be denied to write over that file. Therefore, make sure all data files in the directory one
        is saving to have different names).
    """

    # Set path to my functions:
    import sys

    sys.path.append("../tools/")

    # import libraries:
    from netCDF4 import Dataset, num2date, date2num
    from datetime import datetime
    import netCDF4

    # import functions
    from save_netcdf_fields import add_global_atrributes

    # set dimensions for swh, latitude and longitude
    Nt, Ny, Nx = fp.shape

    # Initiate netCDF file in directory specified by output variable
    nc = Dataset(output, "w", format="NETCDF4")

    # set global attributes:
    add_global_atrributes(nc, summary)

    # set time, longitude and latitude dimensions
    time_dim = nc.createDimension("time", Nt)
    lon_dim = nc.createDimension("lon", Nx)
    lat_dim = nc.createDimension("lat", Ny)

    # Initiate dictionary for attribute for all variables:
    input_vars = {}

    # Set attributes for each variable
    ######## time ########
    input_vars["time"] = {}
    input_vars["time"]["calendar"] = "gregorian"
    input_vars["time"]["units"] = "days since 1900-01-01 00:00:00"
    input_vars["time"]["long_name"] = "Time"
    input_vars["time"]["standard_name"] = "time"
    input_vars["time"]["coverage_content_type"] = "coordinate"

    ######## latitude ########
    input_vars["lat"] = {}
    input_vars["lat"]["units"] = "degrees_north"
    input_vars["lat"]["long_name"] = "Latitude"
    input_vars["lat"]["standard_name"] = "latitude"
    input_vars["lat"]["coverage_content_type"] = "coordinate"

    ######## longitude ########
    input_vars["lon"] = {}
    input_vars["lon"]["units"] = "degrees_east"
    input_vars["lon"]["long_name"] = "Longitude"
    input_vars["lon"]["standard_name"] = "longitude"
    input_vars["lon"]["coverage_content_type"] = "coordinate"

    ######## WSP ########
    input_vars["fp"] = {}
    input_vars["fp"]["units"] = "s-1"
    input_vars["fp"]["long_name"] = "wave peak frequency"
    input_vars["fp"]["standard_name"] = "sea_surface_wave_peak_frequency"
    input_vars["fp"]["coverage_content_type"] = "physicalMeasurement"
    input_vars["fp"]["coordinates"] = "time lat lon"

    # Initialize variables:
    vars = {}

    # Create variable, set fill_value to the netCDF4 default_fillvals with 32 or 64-bit floating point:
    for var in input_vars.keys():

        if var == "time":
            vars["time"] = nc.createVariable(
                "time", "<f4", ("time",), fill_value=netCDF4.default_fillvals["f4"]
            )
        elif var == "lon":
            vars["lon"] = nc.createVariable(
                "lon", "<f4", ("lon",), fill_value=netCDF4.default_fillvals["f4"]
            )
        elif var == "lat":
            vars["lat"] = nc.createVariable(
                "lat", "<f4", ("lat",), fill_value=netCDF4.default_fillvals["f4"]
            )
        else:
            vars[var] = nc.createVariable(
                var,
                "<f8",
                ("time", "lat", "lon"),
                fill_value=netCDF4.default_fillvals["f8"],
            )

        # Set variable attributes
        for a in input_vars[var].keys():
            setattr(vars[var], a, input_vars[var][a])

    # Force netCDF file to convert variables to masked arrays when read
    for var_all in vars.keys():
        vars[var_all].set_auto_maskandscale(True)

    # assign data to netCDF variables
    vars["lat"][:] = lat
    vars["lon"][:] = lon
    vars["time"][:] = date2num(time, nc.variables["time"].units)
    vars["fp"][:] = fp

    # Close netCDF file modifications
    nc.close()


########### Weighted Least-Squares Parameters and Uncertainty Function ###########
def save_netcdf_lsf_parameters(
    a_amp,
    a_phase,
    s_amp,
    s_phase,
    fve,
    a_amp_unc,
    a_phase_unc,
    s_amp_unc,
    s_phase_unc,
    lon,
    lat,
    output,
    summary,
):

    """
    save_netcdf_lsf_parameters(a_amp, a_phase, s_amp, s_phase, fve, a_amp_unc, a_phase_unc, s_amp_unc, s_phase_unc,
                               lon, lat, output, summary)

        Function to save weighted least-squares fit parameters with longitude and latitude
        variables into a netCDF file in the current directory or directory specified in the output variable.

        Parameters
        ----------
        a_amp : numpy 2D masked array of annual cycle amplitude
        a_phase : numpy 2D masked array of annual cycle phase
        s_amp : numpy 2D masked array of semi-annual cycle amplitude
        s_phase : numpy 2D masked array of semi-annual cycle phase
        fve : numpy 2D masked array of fraction of variance explained by model
        a_amp_unc : numpy 2D masked array of annual cycle amplitude uncertainty
        a_phase_unc : numpy 2D masked array of annual cycle phase uncertainty
        s_amp_unc : numpy 2D masked array of semi-annual cycle amplitude uncertainty
        s_phase_unc : numpy 2D masked array of semi-annual cycle phase uncertainty
        lon : numpy array column vector of longitude coordinates
        lat : numpy array column vector of longitude coordinates
        output: filename (path to file and file's name)
            ex: output = '/zdata/downloads/colosi_data_bk/binned_data/ifremer_p1_daily_data/
                          my_daily_binned_ifremer_data/IFREMER_binned_alt_swh_wsp_93_16.nc'
        summary : A string containing the summary of what is contained in the netCDF file.

        Returns
        -------

        Libraries necessary to run function
        -----------------------------------
        from netCDF4 import Dataset, num2date, date2num
        from datetime import datetime
        import netCDF4

        Important Note
        --------------
        The NetCDF file cannot be saved to a directory that has a file with the same name as the file being saved
        (permission will be denied to write over that file. Therefore, make sure all data files in the directory one
        is saving to have different names).
    """

    # Set path to my functions:
    import sys

    sys.path.append("../tools/")

    # import libraries:
    from netCDF4 import Dataset, num2date, date2num
    from datetime import datetime
    import netCDF4

    # import functions
    from save_netcdf_fields import add_global_atrributes

    # set dimensions for parameters, latitude and longitude
    Ny, Nx = a_amp.shape

    # Initiate netCDF file in directory specified by output variable
    nc = Dataset(output, "w", format="NETCDF4")

    # set global attributes:
    add_global_atrributes(nc, summary)

    # set longitude and latitude dimensions
    lon_dim = nc.createDimension("lon", Nx)
    lat_dim = nc.createDimension("lat", Ny)

    # Initiate dictionary for attribute for all variables:
    input_vars = {}

    # Set attributes for each variable
    ######## latitude ########
    input_vars["lat"] = {}
    input_vars["lat"]["units"] = "degrees_north"
    input_vars["lat"]["long_name"] = "Latitude"
    input_vars["lat"]["standard_name"] = "latitude"
    input_vars["lat"]["coverage_content_type"] = "coordinate"

    ######## longitude ########
    input_vars["lon"] = {}
    input_vars["lon"]["units"] = "degrees_east"
    input_vars["lon"]["long_name"] = "Longitude"
    input_vars["lon"]["standard_name"] = "longitude"
    input_vars["lon"]["coverage_content_type"] = "coordinate"

    ######## Annual cycle Amplitude ########
    input_vars["a_amp"] = {}
    input_vars["a_amp"]["units"] = "m"
    input_vars["a_amp"]["long_name"] = "Annual Cycle Amplitude"
    input_vars["a_amp"]["coverage_content_type"] = "physicalMeasurement"
    input_vars["a_amp"]["coordinates"] = "lat lon"

    ######## Annual cycle Phase ########
    input_vars["a_phase"] = {}
    input_vars["a_phase"]["units"] = "months"
    input_vars["a_phase"]["long_name"] = "Annual Cycle Phase"
    input_vars["a_phase"]["coverage_content_type"] = "physicalMeasurement"
    input_vars["a_phase"]["coordinates"] = "lat lon"

    ######## Semi-Annual cycle Amplitude ########
    input_vars["s_amp"] = {}
    input_vars["s_amp"]["units"] = "m"
    input_vars["s_amp"]["long_name"] = "Semi-Annual Cycle Amplitude"
    input_vars["s_amp"]["coverage_content_type"] = "physicalMeasurement"
    input_vars["s_amp"]["coordinates"] = "lat lon"

    ######## Semi-Annual cycle Phase ########
    input_vars["s_phase"] = {}
    input_vars["s_phase"]["units"] = "months"
    input_vars["s_phase"]["long_name"] = "Semi-Annual Cycle Phase"
    input_vars["s_phase"]["coverage_content_type"] = "physicalMeasurement"
    input_vars["s_phase"]["coordinates"] = "lat lon"

    ######## FVE ########
    input_vars["fve"] = {}
    input_vars["fve"]["units"] = "%"
    input_vars["fve"]["long_name"] = "Fraction of Variance Explained by model"
    input_vars["fve"]["coverage_content_type"] = "physicalMeasurement"
    input_vars["fve"]["coordinates"] = "time lat lon"

    ######## Annual cycle Amplitude Uncertainty ########
    input_vars["a_amp_unc"] = {}
    input_vars["a_amp_unc"]["units"] = "m"
    input_vars["a_amp_unc"]["long_name"] = "Annual Cycle Amplitude Uncertainty"
    input_vars["a_amp_unc"]["coverage_content_type"] = "physicalMeasurement"
    input_vars["a_amp_unc"]["coordinates"] = "lat lon"

    ######## Annual cycle Phase Uncertainty ########
    input_vars["a_phase_unc"] = {}
    input_vars["a_phase_unc"]["units"] = "months"
    input_vars["a_phase_unc"]["long_name"] = "Annual Cycle Phase Uncertainty"
    input_vars["a_phase_unc"]["coverage_content_type"] = "physicalMeasurement"
    input_vars["a_phase_unc"]["coordinates"] = "lat lon"

    ######## Semi-Annual cycle Amplitude Uncertainty ########
    input_vars["s_amp_unc"] = {}
    input_vars["s_amp_unc"]["units"] = "m"
    input_vars["s_amp_unc"]["long_name"] = "Semi-Annual Cycle Amplitude Uncertainty"
    input_vars["s_amp_unc"]["coverage_content_type"] = "physicalMeasurement"
    input_vars["s_amp_unc"]["coordinates"] = "lat lon"

    ######## Semi-Annual cycle Phase Uncertainty ########
    input_vars["s_phase_unc"] = {}
    input_vars["s_phase_unc"]["units"] = "months"
    input_vars["s_phase_unc"]["long_name"] = "Semi-Annual Cycle Phase Uncertainty"
    input_vars["s_phase_unc"]["coverage_content_type"] = "physicalMeasurement"
    input_vars["s_phase_unc"]["coordinates"] = "lat lon"

    # Initialize variables:
    vars = {}

    # Create variable, set fill_value to the netCDF4 default_fillvals with 32 or 64-bit floating point:
    for var in input_vars.keys():

        if var == "lon":
            vars["lon"] = nc.createVariable(
                "lon", "<f4", ("lon",), fill_value=netCDF4.default_fillvals["f4"]
            )
        elif var == "lat":
            vars["lat"] = nc.createVariable(
                "lat", "<f4", ("lat",), fill_value=netCDF4.default_fillvals["f4"]
            )
        else:
            vars[var] = nc.createVariable(
                var, "<f8", ("lat", "lon"), fill_value=netCDF4.default_fillvals["f8"]
            )

        # Set variable attributes
        for a in input_vars[var].keys():
            setattr(vars[var], a, input_vars[var][a])

    # Force netCDF file to convert variables to masked arrays when read
    for var_all in vars.keys():
        vars[var_all].set_auto_maskandscale(True)

    # assign data to netCDF variables
    vars["lat"][:] = lat
    vars["lon"][:] = lon
    vars["a_amp"][:] = a_amp
    vars["a_phase"][:] = a_phase
    vars["s_amp"][:] = s_amp
    vars["s_phase"][:] = s_phase
    vars["fve"][:] = fve
    vars["a_amp_unc"][:] = a_amp_unc
    vars["a_phase_unc"][:] = a_phase_unc
    vars["s_amp_unc"][:] = s_amp_unc
    vars["s_phase_unc"][:] = s_phase_unc

    # Close netCDF file modifications
    nc.close()


########### Decorrelation Time Scale Function ###########
def save_netcdf_decor_scale(decor_scale, lon, lat, time, output, summary):

    """
    save_netcdf_decor_scale(decor_scale, lon, lat, time, output, summary)

        Function to save decorrelation time scale with longitude, latitude, and time
        variables into a netCDF file in the current directory or directory specified in the output variable.

        Parameters
        ----------
        decor_scale : numpy 3D masked array of monthly decorrelation time scales
        lon : numpy array column vector of longitude coordinates
        lat : numpy array column vector of longitude coordinates
        time : numpy array of date2num values of time for days
        output: filename (path to file and file's name)
            ex: output = '/zdata/downloads/colosi_data_bk/binned_data/ifremer_p1_daily_data/
                          my_daily_binned_ifremer_data/IFREMER_binned_alt_swh_wsp_93_16.nc'
        summary : A string containing the summary of what is contained in the netCDF file.

        Returns
        -------

        Libraries necessary to run function
        -----------------------------------
        from netCDF4 import Dataset, num2date, date2num
        from datetime import datetime
        import netCDF4

        Important Note
        --------------
        The NetCDF file cannot be saved to a directory that has a file with the same name as the file being saved
        (permission will be denied to write over that file. Therefore, make sure all data files in the directory one
        is saving to have different names).
    """

    # Set path to my functions:
    import sys

    sys.path.append("../tools/")

    # import libraries:
    from netCDF4 import Dataset, num2date, date2num
    from datetime import datetime
    import netCDF4

    # import functions
    from save_netcdf_fields import add_global_atrributes

    # set dimensions for swh, latitude and longitude
    Nt, Ny, Nx = decor_scale.shape

    # Initiate netCDF file in directory specified by output variable
    nc = Dataset(output, "w", format="NETCDF4")

    # set global attributes:
    add_global_atrributes(nc, summary)

    # set time, longitude and latitude dimensions
    time_dim = nc.createDimension("time", Nt)
    lon_dim = nc.createDimension("lon", Nx)
    lat_dim = nc.createDimension("lat", Ny)

    # Initiate dictionary for attribute for all variables:
    input_vars = {}

    # Set attributes for each variable
    ######## time ########
    input_vars["time"] = {}
    input_vars["time"]["calendar"] = "gregorian"
    input_vars["time"]["units"] = "days since 1900-01-01 00:00:00"
    input_vars["time"]["long_name"] = "Time"
    input_vars["time"]["standard_name"] = "time"
    input_vars["time"]["coverage_content_type"] = "coordinate"

    ######## latitude ########
    input_vars["lat"] = {}
    input_vars["lat"]["units"] = "degrees_north"
    input_vars["lat"]["long_name"] = "Latitude"
    input_vars["lat"]["standard_name"] = "latitude"
    input_vars["lat"]["coverage_content_type"] = "coordinate"

    ######## longitude ########
    input_vars["lon"] = {}
    input_vars["lon"]["units"] = "degrees_east"
    input_vars["lon"]["long_name"] = "Longitude"
    input_vars["lon"]["standard_name"] = "longitude"
    input_vars["lon"]["coverage_content_type"] = "coordinate"

    ######## Monthly Decorrelation time scales ########
    input_vars["decor_scale"] = {}
    input_vars["decor_scale"]["units"] = "days"
    input_vars["decor_scale"]["long_name"] = "Monthly Decorrelation time scales"
    input_vars["decor_scale"]["coverage_content_type"] = "physicalMeasurement"
    input_vars["decor_scale"]["coordinates"] = "time lat lon"

    # Initialize variables:
    vars = {}

    # Create variable, set fill_value to the netCDF4 default_fillvals with 32 or 64-bit floating point:
    for var in input_vars.keys():

        if var == "time":
            vars["time"] = nc.createVariable(
                "time", "<f4", ("time",), fill_value=netCDF4.default_fillvals["f4"]
            )
        elif var == "lon":
            vars["lon"] = nc.createVariable(
                "lon", "<f4", ("lon",), fill_value=netCDF4.default_fillvals["f4"]
            )
        elif var == "lat":
            vars["lat"] = nc.createVariable(
                "lat", "<f4", ("lat",), fill_value=netCDF4.default_fillvals["f4"]
            )
        else:
            vars[var] = nc.createVariable(
                var,
                "<f8",
                ("time", "lat", "lon"),
                fill_value=netCDF4.default_fillvals["f8"],
            )

        # Set variable attributes
        for a in input_vars[var].keys():
            setattr(vars[var], a, input_vars[var][a])

    # Force netCDF file to convert variables to masked arrays when read
    for var_all in vars.keys():
        vars[var_all].set_auto_maskandscale(True)

    # assign data to netCDF variables
    vars["lat"][:] = lat
    vars["lon"][:] = lon
    vars["time"][:] = date2num(time, nc.variables["time"].units)
    vars["decor_scale"][:] = decor_scale

    # Close netCDF file modifications
    nc.close()


########### Probability of Swell Function ###########
def save_netcdf_prob_swell(
    monthly_prob_swell,
    seasonal_prob_swell,
    lon,
    lat,
    monthly_time,
    seasonal_time,
    output,
    summary,
):

    """
    save_netcdf_prob_swell(monthly_prob_swell, seasonal_prob_swell, lon, lat, monthly_time, seasonal_time, output,
                           summary):

        Function to save monthly and seasonal probability of swell with longitude, latitude, monthly time and seasonal time
        variables into a netCDF file in the current directory or directory specified in the output variable.

        Parameters
        ----------
        monthly_prob_swell : numpy 3D masked array of monthly progression of probability of swell
        seasonal_prob_swell : numpy 3D masked array of seasonal progression of probability of swell
        lon : numpy array column vector of longitude coordinates
        lat : numpy array column vector of longitude coordinates
        monthly_time : numpy array of date2num values of time for months
        seasonal_time : numpy array of date2num values of time for seasons
        output: filename (path to file and file's name)
            ex: output = '/zdata/downloads/colosi_data_bk/binned_data/ifremer_p1_daily_data/
                          my_daily_binned_ifremer_data/IFREMER_binned_alt_swh_wsp_93_16.nc'
        summary : A string containing the summary of what is contained in the netCDF file.

        Returns
        -------

        Libraries necessary to run function
        -----------------------------------
        from netCDF4 import Dataset, num2date, date2num
        from datetime import datetime
        import netCDF4

        Important Note
        --------------
        The NetCDF file cannot be saved to a directory that has a file with the same name as the file being saved
        (permission will be denied to write over that file. Therefore, make sure all data files in the directory one
        is saving to have different names).
    """

    # Set path to my functions:
    import sys

    sys.path.append("../tools/")

    # import libraries:
    from netCDF4 import Dataset, num2date, date2num
    from datetime import datetime
    import netCDF4

    # import functions
    from save_netcdf_fields import add_global_atrributes

    # set dimensions for swh, latitude and longitude
    Nt_mon, Ny, Nx = monthly_prob_swell.shape
    Nt_sea, Ny, Nx = seasonal_prob_swell.shape

    # Initiate netCDF file in directory specified by output variable
    nc = Dataset(output, "w", format="NETCDF4")

    # set global attributes:
    add_global_atrributes(nc, summary)

    # set time, longitude and latitude dimensions
    time_dim_mon = nc.createDimension("time_mon", Nt_mon)
    time_dim_sea = nc.createDimension("time_sea", Nt_sea)
    lon_dim = nc.createDimension("lon", Nx)
    lat_dim = nc.createDimension("lat", Ny)

    # Initiate dictionary for attribute for all variables:
    input_vars = {}

    # Set attributes for each variable
    ######## monthly time ########
    input_vars["time_mon"] = {}
    input_vars["time_mon"]["calendar"] = "gregorian"
    input_vars["time_mon"]["units"] = "month of the year"
    input_vars["time_mon"]["long_name"] = "Time in months of the year"
    input_vars["time_mon"]["ancillary_variables"] = "monthly_prob_swell"
    input_vars["time_mon"]["standard_name"] = "time"
    input_vars["time_mon"]["coverage_content_type"] = "coordinate"
    input_vars["time_mon"][
        "comment"
    ] = "Monthly time 1-dimensional array contains integer values that identify the month which the monthly probability of swell is computed. For example, the 1 refers to the month of January, 2 refers to the month of February, and so on."

    ######## seasonal time ########
    input_vars["time_sea"] = {}
    input_vars["time_sea"]["calendar"] = "gregorian"
    input_vars["time_sea"]["units"] = "season of the year"
    input_vars["time_sea"]["long_name"] = "Time in seasons of the year"
    input_vars["time_sea"]["ancillary_variables"] = "seasonal_prob_swell"
    input_vars["time_sea"]["standard_name"] = "time"
    input_vars["time_sea"]["coverage_content_type"] = "coordinate"
    input_vars["time_sea"][
        "comment"
    ] = "Seasonal time 1-dimensional array contains interger values that identify the season which the seasonal probability of swell is computed. Integers 1, 2, 3, and 4 correspond to the seasons of winter, spring, summer, and fall respectively."

    ######## latitude ########
    input_vars["lat"] = {}
    input_vars["lat"]["units"] = "degrees_north"
    input_vars["lat"]["long_name"] = "Latitude"
    input_vars["lat"]["standard_name"] = "latitude"
    input_vars["lat"]["coverage_content_type"] = "coordinate"

    ######## longitude ########
    input_vars["lon"] = {}
    input_vars["lon"]["units"] = "degrees_east"
    input_vars["lon"]["long_name"] = "Longitude"
    input_vars["lon"]["standard_name"] = "longitude"
    input_vars["lon"]["coverage_content_type"] = "coordinate"

    ######## monthly Probability of Swell ########
    input_vars["monthly_prob_swell"] = {}
    input_vars["monthly_prob_swell"]["units"] = "%"
    input_vars["monthly_prob_swell"][
        "long_name"
    ] = "Monthly Progression of Probability of Swell"
    input_vars["monthly_prob_swell"]["coverage_content_type"] = "physicalMeasurement"
    input_vars["monthly_prob_swell"]["coordinates"] = "time lat lon"

    ######## monthly Probability of Swell ########
    input_vars["seasonal_prob_swell"] = {}
    input_vars["seasonal_prob_swell"]["units"] = "%"
    input_vars["seasonal_prob_swell"][
        "long_name"
    ] = "Seasonal Progression of Probability of Swell"
    input_vars["seasonal_prob_swell"]["coverage_content_type"] = "physicalMeasurement"
    input_vars["seasonal_prob_swell"]["coordinates"] = "time lat lon"

    # Initialize variables:
    vars = {}

    # Create variable, set fill_value to the netCDF4 default_fillvals with 32 or 64-bit floating point:
    for var in input_vars.keys():

        if var == "time_mon":
            vars["time_mon"] = nc.createVariable(
                "time_mon",
                "<f4",
                ("time_mon",),
                fill_value=netCDF4.default_fillvals["f4"],
            )
        elif var == "time_sea":
            vars["time_sea"] = nc.createVariable(
                "time_sea",
                "<f4",
                ("time_sea",),
                fill_value=netCDF4.default_fillvals["f4"],
            )
        elif var == "lon":
            vars["lon"] = nc.createVariable(
                "lon", "<f4", ("lon",), fill_value=netCDF4.default_fillvals["f4"]
            )
        elif var == "lat":
            vars["lat"] = nc.createVariable(
                "lat", "<f4", ("lat",), fill_value=netCDF4.default_fillvals["f4"]
            )
        elif var == "monthly_prob_swell":
            vars[var] = nc.createVariable(
                var,
                "<f8",
                ("time_mon", "lat", "lon"),
                fill_value=netCDF4.default_fillvals["f8"],
            )
        elif var == "seasonal_prob_swell":
            vars[var] = nc.createVariable(
                var,
                "<f8",
                ("time_sea", "lat", "lon"),
                fill_value=netCDF4.default_fillvals["f8"],
            )

        # Set variable attributes
        for a in input_vars[var].keys():
            setattr(vars[var], a, input_vars[var][a])

    # Force netCDF file to convert variables to masked arrays when read
    for var_all in vars.keys():
        vars[var_all].set_auto_maskandscale(True)

    # assign data to netCDF variables
    vars["lat"][:] = lat
    vars["lon"][:] = lon
    vars["time_mon"][:] = monthly_time
    vars["time_sea"][:] = seasonal_time
    vars["monthly_prob_swell"][:] = monthly_prob_swell
    vars["seasonal_prob_swell"][:] = seasonal_prob_swell

    # Close netCDF file modifications
    nc.close()


########### IFREMER Binned Along Track SWH and WSP Function ###########
def save_netcdf_binned_swh_wsp(swh, wsp, nobs, lon, lat, time, output, summary):

    """
    save_netcdf_binned_swh(swh, wsp, nobs lon, lat, time, output, summary)

        Function to save corrected binned satellite altimeter swh and wsp with longitude, latitude, and time
        variables into a netCDF file in the current directory or directory specified in the output variable.

        Parameters
        ----------
        swh : numpy 3D masked array of IFREMER correct swh
        wsp : numpy 3D masked array of IFREMER correct wsp
        nobs : numpy 3D masked array of number of observations averaged in bins
        lon : numpy array column vector of longitude coordinates
        lat : numpy array column vector of longitude coordinates
        time : numpy array of date2num values of time for days
        output: filename (path to file and file's name)
            ex: output = '/zdata/downloads/colosi_data_bk/binned_data/ifremer_p1_daily_data/
                          my_daily_binned_ifremer_data/IFREMER_binned_alt_swh_wsp_93_16.nc'
        summary : A string containing the summary of what is contained in the netCDF file.

        Returns
        -------

        Libraries necessary to run function
        -----------------------------------
        from netCDF4 import Dataset, num2date, date2num
        from datetime import datetime
        import netCDF4

        Important Note
        --------------
        The NetCDF file cannot be saved to a directory that has a file with the same name as the file being saved
        (permission will be denied to write over that file. Therefore, make sure all data files in the directory one
        is saving to have different names).
    """

    # Set path to my functions:
    import sys

    sys.path.append("../tools/")

    # import libraries:
    from netCDF4 import Dataset, num2date, date2num
    from datetime import datetime
    import netCDF4

    # import functions
    from save_netcdf_fields import add_global_atrributes

    # set dimensions for swh, latitude and longitude
    Nt, Ny, Nx = swh.shape

    # Initiate netCDF file in directory specified by output variable
    nc = Dataset(output, "w", format="NETCDF4")

    # set global attributes:
    add_global_atrributes(nc, summary)

    # set time, longitude and latitude dimensions
    time_dim = nc.createDimension("time", Nt)
    lon_dim = nc.createDimension("lon", Nx)
    lat_dim = nc.createDimension("lat", Ny)

    # Initiate dictionary for attribute for all variables:
    input_vars = {}

    # Set attributes for each variable
    ######## time ########
    input_vars["time"] = {}
    input_vars["time"]["calendar"] = "gregorian"
    input_vars["time"]["units"] = "days since 1900-01-01 00:00:00"
    input_vars["time"]["long_name"] = "Time"
    input_vars["time"]["standard_name"] = "time"
    input_vars["time"]["coverage_content_type"] = "coordinate"

    ######## latitude ########
    input_vars["lat"] = {}
    input_vars["lat"]["units"] = "degrees_north"
    input_vars["lat"]["long_name"] = "Latitude"
    input_vars["lat"]["standard_name"] = "latitude"
    input_vars["lat"]["coverage_content_type"] = "coordinate"

    ######## longitude ########
    input_vars["lon"] = {}
    input_vars["lon"]["units"] = "degrees_east"
    input_vars["lon"]["long_name"] = "Longitude"
    input_vars["lon"]["standard_name"] = "longitude"
    input_vars["lon"]["coverage_content_type"] = "coordinate"

    ######## Number of Observations ########
    input_vars["nobs"] = {}
    input_vars["nobs"]["long_name"] = "Number of Observations"
    input_vars["nobs"]["standard_name"] = "number_of_observations"
    input_vars["nobs"]["ancillary_variables"] = "swh"
    input_vars["nobs"]["coverage_content_type"] = "physicalMeasurement"
    input_vars["nobs"]["coordinates"] = "time lat lon"
    input_vars["nobs"][
        "comment"
    ] = "Provides a 3-dimensional array of the number of data points average in 1 by 1 degree bins."

    ######## SWH ########
    input_vars["swh"] = {}
    input_vars["swh"]["units"] = "m"
    input_vars["swh"][
        "long_name"
    ] = "binned corrected altimeter significant wave height"
    input_vars["swh"]["standard_name"] = "sea_surface_wave_significant_height"
    input_vars["swh"]["coverage_content_type"] = "physicalMeasurement"
    input_vars["swh"]["coordinates"] = "time lat lon"

    ######## WSP ########
    input_vars["wsp"] = {}
    input_vars["wsp"]["units"] = "m s-1"
    input_vars["wsp"]["long_name"] = "binned corrected altimeter wind speed"
    input_vars["wsp"]["standard_name"] = "wind_speed"
    input_vars["wsp"]["coverage_content_type"] = "physicalMeasurement"
    input_vars["wsp"]["coordinates"] = "time lat lon"

    # Initialize variables:
    vars = {}

    # Create variable, set fill_value to the netCDF4 default_fillvals with 32 or 64-bit floating point:
    for var in input_vars.keys():

        if var == "time":
            vars["time"] = nc.createVariable(
                "time", "<f4", ("time",), fill_value=netCDF4.default_fillvals["f4"]
            )
        elif var == "lon":
            vars["lon"] = nc.createVariable(
                "lon", "<f4", ("lon",), fill_value=netCDF4.default_fillvals["f4"]
            )
        elif var == "lat":
            vars["lat"] = nc.createVariable(
                "lat", "<f4", ("lat",), fill_value=netCDF4.default_fillvals["f4"]
            )
        else:
            vars[var] = nc.createVariable(
                var,
                "<f8",
                ("time", "lat", "lon"),
                fill_value=netCDF4.default_fillvals["f8"],
            )

        # Set variable attributes
        for a in input_vars[var].keys():
            setattr(vars[var], a, input_vars[var][a])

    # Force netCDF file to convert variables to masked arrays when read
    for var_all in vars.keys():
        vars[var_all].set_auto_maskandscale(True)

    # assign data to netCDF variables
    vars["lat"][:] = lat
    vars["lon"][:] = lon
    vars["time"][:] = date2num(time, nc.variables["time"].units)
    vars["nobs"][:] = nobs
    vars["swh"][:] = swh
    vars["wsp"][:] = wsp

    # Close netCDF file modifications
    nc.close()


########### ETOPO1 Deresolved Elevation Function ###########
def save_netcdf_deresolved_topo(topo, lon, lat, output, summary):

    """
    save_netcdf_deresolved_topo(topo, land_mask, lon, lat, output, summary)

        Function to save ETOPO1 integrated topography and bathymetry with longitude, latitude, and time
        variables into a netCDF file in the current directory or directory specified in the output variable.

        Parameters
        ----------
        topo : numpy 3D masked array of elevation
        land_mask : numpy 3D masked array of designating coordinates of land
        lon : numpy array column vector of longitude coordinates
        lat : numpy array column vector of longitude coordinates
        output: filename (path to file and file's name)
            ex: output = '/zdata/downloads/colosi_data_bk/binned_data/ifremer_p1_daily_data/
                          my_daily_binned_ifremer_data/IFREMER_binned_alt_swh_wsp_93_16.nc'
        summary : A string containing the summary of what is contained in the netCDF file.

        Returns
        -------

        Libraries necessary to run function
        -----------------------------------
        from netCDF4 import Dataset, num2date, date2num
        from datetime import datetime
        import netCDF4

        Important Note
        --------------
        The NetCDF file cannot be saved to a directory that has a file with the same name as the file being saved
        (permission will be denied to write over that file. Therefore, make sure all data files in the directory one
        is saving to have different names).
    """

    # Set path to my functions:
    import sys

    sys.path.append("../tools/")

    # import libraries:
    from netCDF4 import Dataset, num2date, date2num
    from datetime import datetime
    import netCDF4

    # import functions
    from save_netcdf_fields import add_global_atrributes

    # set dimensions for swh, latitude and longitude
    Nt, Ny, Nx = topo.shape

    # Initiate netCDF file in directory specified by output variable
    nc = Dataset(output, "w", format="NETCDF4")

    # set global attributes:
    add_global_atrributes(nc, summary)

    # set longitude and latitude dimensions
    lon_dim = nc.createDimension("lon", Nx)
    lat_dim = nc.createDimension("lat", Ny)

    # Initiate dictionary for attribute for all variables:
    input_vars = {}

    # Set attributes for each variable

    ######## latitude ########
    input_vars["lat"] = {}
    input_vars["lat"]["units"] = "degrees_north"
    input_vars["lat"]["long_name"] = "Latitude"
    input_vars["lat"]["standard_name"] = "latitude"
    input_vars["lat"]["coverage_content_type"] = "coordinate"

    ######## longitude ########
    input_vars["lon"] = {}
    input_vars["lon"]["units"] = "degrees_east"
    input_vars["lon"]["long_name"] = "Longitude"
    input_vars["lon"]["standard_name"] = "longitude"
    input_vars["lon"]["coverage_content_type"] = "coordinate"

    ######## Topography ########
    input_vars["topo"] = {}
    input_vars["topo"]["units"] = "m"
    input_vars["topo"]["long_name"] = "Earth surface elevation"
    input_vars["topo"]["coverage_content_type"] = "physicalMeasurement"
    input_vars["topo"]["coordinates"] = "lat lon"

    ######## land_mask ########
    input_vars["land_mask"] = {}
    input_vars["land_mask"]["units"] = "boolean"
    input_vars["land_mask"]["long_name"] = "coordinates of land"
    input_vars["land_mask"]["coverage_content_type"] = "physicalMeasurement"
    input_vars["land_mask"]["coordinates"] = "lat lon"

    # Initialize variables:
    vars = {}

    # Create variable, set fill_value to the netCDF4 default_fillvals with 32 or 64-bit floating point:
    for var in input_vars.keys():

        if var == "lon":
            vars["lon"] = nc.createVariable(
                "lon", "<f4", ("lon",), fill_value=netCDF4.default_fillvals["f4"]
            )
        elif var == "lat":
            vars["lat"] = nc.createVariable(
                "lat", "<f4", ("lat",), fill_value=netCDF4.default_fillvals["f4"]
            )
        else:
            vars[var] = nc.createVariable(
                var, "<f8", ("lat", "lon"), fill_value=netCDF4.default_fillvals["f8"]
            )

        # Set variable attributes
        for a in input_vars[var].keys():
            setattr(vars[var], a, input_vars[var][a])

    # Force netCDF file to convert variables to masked arrays when read
    for var_all in vars.keys():
        vars[var_all].set_auto_maskandscale(True)

    # assign data to netCDF variables
    vars["lat"][:] = lat
    vars["lon"][:] = lon
    vars["topo"][:] = topo
    vars["land_mask"][:] = land_mask

    # Close netCDF file modifications
    nc.close()
