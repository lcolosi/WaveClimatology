# WAVEWATCH III Probability of Swell program

# Path to access python functions
import sys

sys.path.append("../tools/")

# Path to access raw data (CHANGE PATH IF DATA IS NOT in /data DIRECTORY)
data_path = "../data/"

# Import Libraries
import numpy as np
from netCDF4 import Dataset, num2date
import glob
from collections import OrderedDict

# My functions
import cartopy_figs as cart
from save_netcdf_fields import add_global_atrributes, save_netcdf_prob_swell

# Set time extents for data set
initial_year = 1993
final_year = 2015

# Initialize variables
filenames_fp = []
filenames_wnd = []

# create year loop:
for iyear in range(initial_year, final_year + 1, 1):

    # obtain filenames:
    files_fp = sorted(
        glob.glob(data_path + "fp/WW3-GLOB-30M_" + "%s" % iyear + "*_fp.nc")
    )
    files_wnd = sorted(
        glob.glob(data_path + "Wnd/WW3-GLOB-30M_" + "%s" % iyear + "*_wnd.nc")
    )

    # concatinate all filenames into a list
    filenames_fp += files_fp
    filenames_wnd += files_wnd

# Set longitude and latitude vectors
######## wnd ########
nc_wnd = Dataset(filenames_wnd[0], "r")
lon_wnd = nc_wnd.variables["longitude"][:]
lat_wnd = nc_wnd.variables["latitude"][:]

######## fp ########
nc_fp = Dataset(filenames_fp[0], "r")
lon_fp = nc_fp.variables["longitude"][:]
lat_fp = nc_fp.variables["latitude"][:]

# Set dimensions
nt, nlat, nlon = len(filenames_wnd), len(lat_wnd), len(lon_wnd)

# Initialize variables and constants
N_swell = np.ma.masked_all([nt, nlat, nlon])
N_total = np.ma.masked_all([nt, nlat, nlon])
time = []
g = 9.81

# Loop through filenames:
for ifile in range(nt):

    # Set nc variables:
    nc_wnd = Dataset(filenames_wnd[ifile], "r")
    nc_fp = Dataset(filenames_fp[ifile], "r")

    # call wind velocity components, peak frequency, and time:
    uwnd = nc_wnd.variables["uwnd"][:]
    vwnd = nc_wnd.variables["vwnd"][:]
    fp = nc_fp.variables["fp"][:]
    itime = num2date(nc_fp.variables["time"][:], nc_fp.variables["time"].units)

    # close netCDF file
    nc_wnd.close()
    nc_fp.close()

    # compute wind speed:
    wsp = np.sqrt((uwnd ** 2) + (vwnd ** 2))

    # save time:
    time.append(itime[0])

    # Compute phase speed
    cp = g / (2 * np.pi * fp)

    # compute wave age:
    wave_age = cp / wsp

    # Separate Swell and Wind-Sea Events
    swell_ind = wave_age > 1.2

    # count the number of swell events and wind_sea events
    N_s = np.sum(swell_ind, axis=0)
    N_t = np.ma.ones([nlat, nlon]) * np.shape(swell_ind)[0]

    # save number of swell events and total events
    N_swell[ifile, :, :] = N_s
    N_total[ifile, :, :] = N_t

# Compute probability of swell:

# intialize ordered dictionary
prob_s = OrderedDict()
prob_s["p_month"], prob_s["p_season"] = [], []

######## Seasonally ########

# initialize time indices
months = np.array([m.month for m in time])

# loop through seasons DJF, MAM, JJA, SON in time series:
# initialize season variables
seasons = [[12, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]]

# initialize seasonal loop
for s in range(0, 4):

    # call season:
    season = seasons[s]

    # initialize season index
    bool_1 = months == season[0]
    bool_2 = months == season[1]
    bool_3 = months == season[2]
    ind_1 = bool_1 | bool_2
    ind = ind_1 | bool_3

    # call wave age from season
    N_swell_s = N_swell[ind, :, :]
    N_total_s = N_total[ind, :, :]

    # count the total number of events
    N_swell_t = np.sum(N_swell_s, axis=0)
    N_total_t = np.sum(N_total_s, axis=0)

    # Compute probability of swell
    p_season = N_swell_t / N_total_t

    # Save in dictionary:
    prob_s["p_season"].append(p_season)

# Convert list of 2D arrays to a single 3D array for probability of swell:
prob_swell_s = np.transpose(np.ma.dstack(prob_s["p_season"]), (2, 0, 1))

######## Monthly ########

# initialize time indices
months = np.array([m.month for m in time])

# loop through months in time series:
# initialize variables
month = np.arange(1, 13, 1)

# initialize month loop
for m in range(0, 12):

    # call month:
    imonth = month[m]

    # initialize month indices
    ind = months == imonth

    # call wave age from season
    N_swell_m = N_swell[ind, :, :]
    N_total_m = N_total[ind, :, :]

    # count the total number of events
    N_swell_t = np.sum(N_swell_m, axis=0)
    N_total_t = np.sum(N_total_m, axis=0)

    # Compute probability of swell
    p_month = N_swell_t / N_total_t

    # Save in dictionary:
    prob_s["p_month"].append(p_month)

# Convert list of 2D arrays to a single 3D array for probability of swell:
prob_swell_m = np.transpose(np.ma.dstack(prob_s["p_month"]), (2, 0, 1))

# Set summary variable:
output = "../data/prob_swell/WW3_probability_swell.nc"
summary = "Data contained in this netCDF file is derived from the wave hindcast peak wave frequency (fp) and wind speed (WSP) product produced by French Research Institute for Exploitation of the Sea (IFREMER) using the WAVE-height, WATer depth and Current Hindcasting III (WW3) wave model forced by Climate Forecast System Reanalysis (CFSR) winds (ftp://ftp.ifremer.fr/ifremer/ww3/HINDCAST). Thus, this data is an intermediate product. Here, the wave age, defined as the ratio between peak phase speed and wind speed 10 meters above the ocean surface, is calculated across the globe from 66N to 66S over the time period of 1 January 1993 to 31 December 2015 and is used to categorize wave fields as either dominated by swell or wind-sea. For wave age less than or equal to 1.2, the wave field is considered to be dominated by wind-sea. For wave age greater than 1.2, the wave field is considered to be dominated by swell. Probability of swell illustrate the fraction of time the wave field is swell-dominated relative to the total number of wave measurements. Probability of swell is computed seasonally and monthly and are stored in a 3-dimensional (time, latitude, longitude) masked arrays."

# Save seasonal and monthly Probability of Swell data
save_netcdf_prob_swell(
    monthly_prob_swell=prob_swell_m,
    seasonal_prob_swell=prob_swell_s,
    lon=lon_wnd,
    lat=lat_wnd,
    monthly_time=np.arange(1, 13, 1),
    seasonal_time=np.arange(1, 5, 1),
    output=output,
    summary=summary,
)
