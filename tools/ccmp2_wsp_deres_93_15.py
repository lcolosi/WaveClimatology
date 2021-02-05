# CCMP2 Deresolution WSP program

# Path to access python functions
import sys

sys.path.append("../tools/")

# Path to access raw data (CHANGE PATH IF DATA IS NOT in ../data DIRECTORY)
data_path = "../data/"

# Import Libraries
import numpy as np
from netCDF4 import Dataset, num2date
import glob

# my functions
from averaging_stats import running_mean
from save_netcdf_fields import add_global_atrributes, save_netcdf_deresolved_wsp

# Set dimensions
nt, nlon, nlat = 8400, 360, 133
initial_year = 1993
final_year = 2015

# Set filenames
filenames = []

# create year loop:
for iyear in range(initial_year, final_year + 1, 1):

    # obtain filenames from all days within the directory:
    files = sorted(
        glob.glob(
            data_path + "CCMP_Wind_Analysis_" + "%s" % iyear + "*_V02.0_L3.0_RSS.nc"
        )
    )

    # concatinate all filenames into a list
    filenames += files

# Set longitude and latitude vectors
nc = Dataset(filenames[0], "r")
lon = nc.variables["longitude"][:]
lat = nc.variables["latitude"][:]

# Deresolve longitude and latitude
lon_dr = running_mean(lon, k_dim=[4, 1], task="deresolve")
lat_dr = running_mean(lat, k_dim=[4, 1], task="deresolve")

# Set resolution of longitude and latitude. Also set the desired indices for where latitude will be truncated:
dlon = abs(lon_dr[1] - lon_dr[0])
dlat = abs(lat_dr[1] - lat_dr[0])
lat_min, lat_max = -66, 66

# Deresolve WSP

# initialize variables
wsp_ccmp_d = np.ma.masked_all([nt, nlat, nlon])
time_d = []

# initialize counter:
i = 0

# Loop through filenames:
for f in filenames:

    # Set nc variable:
    nc = Dataset(f, "r")

    # Call wind components adn time:
    uwnd = nc.variables["uwnd"][:]
    vwnd = nc.variables["vwnd"][:]
    itime = num2date(nc.variables["time"][:], nc.variables["time"].units)

    # Compute wind speed:
    wsp_h = np.sqrt((uwnd ** 2) + (vwnd ** 2))

    # Average temporally
    wsp_d = np.mean(wsp_h, axis=0)

    # Decrease the resolution of the the wsp matrix via convolution:
    wsp_conv = running_mean(data=wsp_d, k_dim=[4, 4], task="deresolve")

    # Extract wsp data from -66 to 66 degrees:
    trunc = abs(int((min(lat_dr) - lat_min) / dlat))
    wsp_c = wsp_conv[trunc : len(lat_dr) - trunc, :]

    # Save wsp and time data:
    wsp_ccmp_d[i, :, :] = wsp_c
    time_d.append(itime[0])

    # counter sum
    i += 1

# Extract latitudes from -66 to 66 degrees:
lat_dr = lat_dr[trunc : len(lat_dr) - trunc]

# Loop through years to save deresolved WSP data
# set year time arrays that correspond to the year at which wsp data was collected
years = np.array([y.year for y in time_d])

# Convert time_d to an array:
time_array = np.array(time_d)

# loop through years and intialize yearly variables for binned SWH:
for iyear in np.unique(years):

    # Initialize index for ith year
    ind_year = years == iyear

    # Call ith year for deresolved wsp and time
    wsp_iyear = wsp_ccmp_d[ind_year, :, :]
    time_iyear = time_array[ind_year]

    # Set summary variable:
    summary = (
        "Data contained in this netCDF file is derived from the Cross Calibrated Multi-Platform version 2 (CCMP2) wind vector analysis produced by Remote Sensing Systems product (available at www.remss.com). Thus, this data is an intermediate product. Here, CCMP2 wind speed (WSP) is averaged in 1 degree longitude by 1 degree latitude bins from 66N to 66S over the time period of 1 January %s"
        % iyear
        + " to 31 December %s" % iyear
        + " at daily intervals. Data points are averaged using a running mean. Deresolved WSP data is stored in a 3-dimensional (time, latitude, longitude) masked array."
    )

    # Save ith year binned swh data
    save_netcdf_deresolved_wsp(
        wsp=wsp_iyear,
        lon=lon_dr,
        lat=lat_dr,
        time=time_iyear,
        output="../data/ccmp2_wsp/CCMP2_deresolved_wsp_%s" % iyear + ".nc",
        summary=summary,
    )
