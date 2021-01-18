# Ifremer binning SWH and WSP program

# Path to access python functions
import sys

sys.path.append("../tools/")

# Path to access raw data (CHANGE PATH IF DATA IS NOT in /data DIRECTORY)
data_path = "../data/"

# Import Libraries
import numpy as np
from datetime import datetime, timedelta
import glob
from netCDF4 import Dataset, num2date

# Import functions
from save_netcdf_fields import (
    add_global_atrributes,
    save_netcdf_binned_swh,
    save_netcdf_binned_swh_wsp,
)

# Set filenames for IFREMER product 1 along-track altimeter SWH and WSP daily data

# Initialize variables:
initial_year = 1993
final_year = 2015
filenames = []

# create year loop:
for iyear in range(initial_year, final_year + 1, 1):

    # obtain filenames from all days within the directory:
    files = sorted(glob.glob(data_path + "wm_%s" % iyear + "*.nc"))

    # concatinate all filenames into a list
    filenames += files

# Initialize datetime time series:
time = np.array([datetime(1993, 1, 1) + timedelta(days=i) for i in range(8400)])

# Set dimensions for data of space and time
regular_grid = 1.0
nlon, nlat = int(360.0 / regular_grid), int(180.0 / regular_grid)
grid_lat, grid_lon = np.arange(-90, 90, regular_grid), np.arange(0, 360, regular_grid)
nt = len(time)

# Initialize 3D masked arrays for uncorrected SWH, corrected SWH and WSP, time and number of observations
swh_array = np.ma.masked_all([nt, nlat, nlon])
swhcor_array = np.ma.masked_all([nt, nlat, nlon])
wspcor_array = np.ma.masked_all([nt, nlat, nlon])
N_array = np.ma.masked_all([nt, nlat, nlon])
time_array = []

# set day counter
iday = 0

# Loop through daily files
for f in filenames:

    # initialize dictionary for each file as masked arrays:
    data_base = {}
    data_base["swh"] = np.ma.zeros([nlat, nlon])
    data_base["swhcor"] = np.ma.zeros([nlat, nlon])
    data_base["wspcor"] = np.ma.zeros([nlat, nlon])
    data_base["time"] = np.ma.array(nt)
    data_base["lat"] = grid_lat
    data_base["lon"] = grid_lon
    data_base["N"] = np.ma.zeros([nlat, nlon])

    # Read data from ith file:
    nc = Dataset(f, "r")

    # call SWH, WSP, latitude, and longitude along-track data
    swh = nc.variables["swh"][:]
    swhcor = nc.variables["swhcor"][:]
    wspcor = nc.variables["wind_speed_cor"][:]
    lat = nc.variables["lat"][:]
    lon = nc.variables["lon"][:]

    # Center longitude around the pacific by adding 360 to all negative longitude values
    lon[lon < 0] = lon[lon < 0] + 360

    # Call time data
    time = num2date(nc.variables["time"][:], nc.variables["time"].units)

    # Create an array containing the data (year, month, and day) that each observation was recorded
    days = np.array([t.date() for t in time])

    # Obtain the unique year, month, and day
    day = np.unique(days)[0]

    # Loop through each along track data point:
    for idata in range(len(swh)):

        # set lon and lat geographic indices:
        indlat = np.abs(grid_lat - lat[idata]).argmin()
        indlon = np.abs(grid_lon - lon[idata]).argmin()
        # Sum data point at specified index:
        data_base["swh"][indlat, indlon] += swh[idata]
        data_base["swhcor"][indlat, indlon] += swhcor[idata]
        data_base["wspcor"][indlat, indlon] += wspcor[idata]
        data_base["N"][indlat, indlon] += 1

    # Mask all grid point where no data point are assigned
    N_mask = np.ma.masked_where(data_base["N"] == 0, data_base["N"])
    swh_mask = np.ma.masked_where(data_base["swh"] == 0, data_base["swh"])
    swhcor_mask = np.ma.masked_where(data_base["swhcor"] == 0, data_base["swhcor"])
    wspcor_mask = np.ma.masked_where(data_base["wspcor"] == 0, data_base["wspcor"])

    # Average of data points that have been summed together on the same grid point and save in 3D array
    swh_array[iday, :, :] = swh_mask / N_mask
    swhcor_array[iday, :, :] = swhcor_mask / N_mask
    wspcor_array[iday, :, :] = wspcor_mask / N_mask
    N_array[iday, :, :] = N_mask

    # Save time as a datetime object:
    time_array.append(datetime(day.year, day.month, day.day, 0, 0))

    # set counter
    iday += 1

    # print file name to keep track of the progresso of the function
    print(f)

# Extract data from between -66 to 66 degrees lattude:
ii = np.where(abs(grid_lat) <= 66)[0]
swh_array = swh_array[:, ii, :]
swhcor_array = swhcor_array[:, ii, :]
wspcor_array = wspcor_array[:, ii, :]
N_array = N_array[:, ii, :]
grid_lat = grid_lat[ii]

# Loop through years to save binned swhcor data:

# set year time arrays that correspond to the year at which swh data point was collected
years = np.array([y.year for y in time_array])

# Convert time_array to an array:
time_array = np.ma.array(time_array)

# loop through years and intialize yearly variables for binned SWH:
for iyear in np.unique(years):

    # Initialize index for ith year
    ind_year = years == iyear

    # Call ith year for binned swh, nobs and time
    swh_iyear = swhcor_array[ind_year, :, :]
    nobs_iyear = N_array[ind_year, :, :]
    time_iyear = time_array[ind_year]

    # Set summary variable:
    summary = (
        "Data contained in this netCDF file is derived from the French Research Institute for Exploitation of the Sea (IFREMER) cross-calibrated along-track satellite altimetry significant wave height (SWH) product (ftp://ftp.ifremer.fr/ifremer/cersat/products/swath/altimeters/waves). Thus, this data is an intermediate product. Here, IFREMER corrected along-track SWH (swhcor) is placed in 1 degree longitude by 1 degree latitude bins from 66N to 66S over the time period of 1 January %s"
        % iyear
        + " to 31 December %s" % iyear
        + " at daily intervals. Data points in bins are averaged and number of observations averaged in bins is recorded. Binned SWH data is stored in a 3-dimensional (time, latitude, longitude) masked array."
    )

    # Save ith year binned swh data
    save_netcdf_binned_swh(
        swh=swh_iyear,
        nobs=nobs_iyear,
        lon=grid_lon,
        lat=grid_lat,
        time=time_iyear,
        output="../data/ifremer_swh/IFREMER_binned_alt_swh_%s" % iyear + ".nc",
        summary=summary,
    )


# Save SWH and WSP into one large file:

# Set summary variable:
summary = "Data contained in this netCDF file is derived from the French Research Institute for Exploitation of the Sea (IFREMER) cross-calibrated along-track satellite altimetry significant wave height (SWH) and wind speed (WSP) product (ftp://ftp.ifremer.fr/ifremer/cersat/products/swath/altimeters/waves). Thus, this data is an intermediate product. Here, IFREMER corrected along-track SWH (swhcor) and WSP are placed in 1 degree longitude by 1 degree latitude bins from 66N to 66S over the time period of 1 January 1993 to 31 December 2015 at daily intervals. Data points in bins are averaged and number of observations averaged in bins is recorded. Binned SWH and WSP data are stored in a 3-dimensional (time, latitude, longitude) masked arrays."

# Save data
save_netcdf_binned_swh_wsp(
    swhcor_array,
    wspcor_array,
    N_array,
    grid_lon,
    grid_lat,
    time_array,
    output="../data/ifremer_swh/IFREMER_binned_alt_wsp_swh_1993_2015.nc",
    summary=summary,
)
