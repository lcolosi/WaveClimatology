# WW3 Deresolution SWH (hs) and WSP program

# Path to access python functions
import sys

sys.path.append("../tools/")

# Path to access raw data (CHANGE PATH IF DATA IS NOT in /data DIRECTORY)
data_path = "../data/"

# Import Libraries
import numpy as np
from netCDF4 import Dataset, num2date
import glob
import datetime

# my functions
from averaging_stats import running_mean
from data_processing import shift_grid
from save_netcdf_fields import (
    add_global_atrributes,
    save_netcdf_deresolved_hs,
    save_netcdf_deresolved_wsp,
)

# Set dimensions
nt, nlon, nlat = 8400, 360, 133
initial_year = 1993
final_year = 2015

# Initialize variables
filenames_hs = []
filenames_wnd = []

# create year loop:
for iyear in range(initial_year, final_year + 1, 1):

    # obtain filenames:
    # Hs:
    files = sorted(glob.glob(data_path + "WW3-GLOB-30M_" + "%s" % iyear + "*_hs.nc"))
    filenames_hs += files
    # Wnd:
    files = sorted(glob.glob(data_path + "WW3-GLOB-30M_" + "%s" % iyear + "*_wnd.nc"))
    filenames_wnd += files


# Set longitude and latitude vectors
######## hs ########
nc_hs = Dataset(filenames_hs[0], "r")
lon_hs = nc_hs.variables["longitude"][:]
lat_hs = nc_hs.variables["latitude"][:]

######## wnd ########
nc_wnd = Dataset(filenames_wnd[0], "r")
lon_wnd = nc_wnd.variables["longitude"][:]
lat_wnd = nc_wnd.variables["latitude"][:]

# Deresolve longitude and latitude
######## hs ########
lon_hs_dr = np.floor(running_mean(lon_hs, k_dim=[2, 1], task="deresolve"))
lat_hs_dr = np.floor(running_mean(lat_hs, k_dim=[2, 1], task="deresolve"))

######## wnd ########
lon_wnd_dr = np.floor(running_mean(lon_wnd, k_dim=[2, 1], task="deresolve"))
lat_wnd_dr = np.floor(running_mean(lat_wnd, k_dim=[2, 1], task="deresolve"))


# Set resolution of longitude and latitude. Also set the desired indices for where latitude will be truncated:
dlon = abs(lon_hs_dr[1] - lon_hs_dr[0])
dlat = abs(lat_hs_dr[1] - lat_hs_dr[0])
lat_min, lat_max = -66, 66

######## Deresolve Hs ########
# initialize masked array:
hs_ww3_cfsr_d = np.ma.masked_all([nt, nlat, nlon])
time_d_hs = []

# initialize counter:
i = 0

# Loop through filenames:
for f in filenames_hs:

    # Set nc variable:
    nc_hs = Dataset(f, "r")

    # call swh and time:
    hs = nc_hs.variables["hs"][:]
    itime = num2date(nc_hs.variables["time"][:], nc_hs.variables["time"].units)

    # Initialize loop variable:
    days = np.array([d.day for d in itime])

    # Loop through days:
    for iday in np.unique(days):

        # Initialize index for ith year
        ind_day = days == iday

        # call data:
        hs_h = hs[ind_day, :, :]

        # Average temporally
        hs_d = np.ma.mean(hs_h, axis=0)

        # Decrease the resolution of the the wsp matrix via convolution:
        hs_conv = running_mean(data=hs_d, k_dim=[2, 2], task="deresolve")

        # Extract hs data from -66 to 66 degrees:
        trunc = abs(int((min(lat_hs_dr) - lat_min) / dlat))
        hs_c = hs_conv[trunc : len(lat_hs_dr) - trunc - 2, :]

        # Center the data over the Pacific:
        hs_shift, lon_shift = shift_grid(hs_c, lon_hs_dr, dlon=180)

        # Save wsp and time data:
        hs_ww3_cfsr_d[i, :, :] = hs_shift
        time_d_hs.append(itime[0])

        # counter sum
        i += 1
        print(i)


######## Deresolve WSP ########
# initialize masked array:
wsp_ww3_cfsr_d = np.ma.masked_all([nt, nlat, nlon])
time_d_wnd = []

# initialize counter:
i = 0

# Loop through filenames
for f in filenames_wnd:

    # set nc variable :
    nc_wnd = Dataset(f, "r")

    # call wind velocity components and time:
    uwnd = nc_wnd.variables["uwnd"][:]
    vwnd = nc_wnd.variables["vwnd"][:]
    itime = num2date(nc_wnd.variables["time"][:], nc_wnd.variables["time"].units)

    # compute wind speed:
    wsp = np.sqrt((uwnd ** 2) + (vwnd ** 2))

    # Initialize loop variable:
    days = np.array([d.day for d in itime])

    # Loop through days:
    for iday in np.unique(days):

        # Initialize index for ith year
        ind_day = days == iday

        # call data:
        wsp_h = wsp[ind_day, :, :]

        # Average temporally
        wsp_d = np.ma.mean(wsp_h, axis=0)

        # Decrease the resolution of the the wsp matrix via convolution:
        wsp_conv = running_mean(data=wsp_d, k_dim=[2, 2], task="deresolve")

        # Extract hs data from -66 to 66 degrees:
        trunc = abs(int((min(lat_wnd_dr) - lat_min) / dlat))
        wsp_c = wsp_conv[trunc : len(lat_wnd_dr) - trunc - 2, :]

        # Center the data over the Pacific:
        wsp_shift, lon_shift = shift_grid(wsp_c, lon_wnd_dr, dlon=180)

        # save the average daily array in a 3D array:
        wsp_ww3_cfsr_d[i, :, :] = wsp_shift
        time_d_wnd.append(itime[0])

        # counter sum
        i += 1
        print(i)

# Extract latitudes from -66 to 66 degrees:
lat_hs_dr = lat_hs_dr[trunc : len(lat_hs_dr) - trunc - 2]
lat_wnd_dr = lat_wnd_dr[trunc : len(lat_wnd_dr) - trunc - 2]

# Loop through years to save deresolved Hs and Wsp data

######## Hs ########
# set year time arrays that correspond to the year at which wsp data was collected
years = np.array([y.year for y in time_d_hs])

# Convert time_d to an array:
time_array = np.array(time_d_hs)

# loop through years and intialize yearly variables for binned SWH:
for iyear in np.unique(years):

    # Initialize index for ith year
    ind_year = years == iyear

    # Call ith year for deresolved wsp and time
    hs_iyear = hs_ww3_cfsr_d[ind_year, :, :]
    time_iyear = time_array[ind_year]

    # Set summary variable:
    summary = (
        "Data contained in this netCDF file is derived from the wave hindcast product produced by French Research Institute for Exploitation of the Sea (IFREMER) using the WAVE-height, WATer depth and Current Hindcasting III (WW3) wave model forced by Climate Forecast System Reanalysis (CFSR) winds (ftp://ftp.ifremer.fr/ifremer/ww3/HINDCAST). Thus, this data is an intermediate product. Here, WW3 wind speed (WSP) is averaged in 1 degree longitude by 1 degree latitude bins from 66N to 66S over the time period of 1 January %s"
        % iyear
        + " to 31 December %s" % iyear
        + " at daily intervals. Data points are averaged using a running mean. Deresolved WSP data is stored in a 3-dimensional (time, latitude, longitude) masked array."
    )

    # Save ith year binned swh data
    save_netcdf_deresolved_hs(
        hs=hs_iyear,
        lon=lon_hs_dr,
        lat=lat_hs_dr,
        time=time_iyear,
        output="../data/ww3_swh/ww3_deresolved_hs_%s" % iyear + ".nc",
        summary=summary,
    )


######## WSP ########
# set year time arrays that correspond to the year at which wsp data was collected
years = np.array([y.year for y in time_d_wnd])

# Convert time_d to an array:
time_array = np.array(time_d_wnd)

# loop through years and intialize yearly variables for binned SWH:
for iyear in np.unique(years):

    # Initialize index for ith year
    ind_year = years == iyear

    # Call ith year for deresolved wsp and time
    wsp_iyear = wsp_ww3_cfsr_d[ind_year, :, :]
    time_iyear = time_array[ind_year]

    # Set summary variable:
    summary = (
        "Data contained in this netCDF file is derived from the wave hindcast product produced by French Research Institute for Exploitation of the Sea (IFREMER) using the WAVE-height, WATer depth and Current Hindcasting III (WW3) wave model forced by Climate Forecast System Reanalysis (CFSR) winds (ftp://ftp.ifremer.fr/ifremer/ww3/HINDCAST). Thus, this data is an intermediate product. Here, WW3 significant wave height (Hs) is averaged in 1 degree longitude by 1 degree latitude bins from 66N to 66S over the time period of 1 January %s"
        % iyear
        + " to 31 December %s" % iyear
        + " at daily intervals. Data points are averaged using a running mean. Deresolved Hs data is stored in a 3-dimensional (time, latitude, longitude) masked array."
    )

    # Save ith year binned swh data
    save_netcdf_deresolved_wsp(
        wsp=wsp_iyear,
        lon=lon_wnd_dr,
        lat=lat_wnd_dr,
        time=time_iyear,
        output="../data/ww3_wsp/ww3_deresolved_wsp_%s" % iyear + ".nc",
        summary=summary,
    )
