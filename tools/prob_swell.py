# WAVEWATCH III Probability of Swell program

# Path to access python functions
import sys

sys.path.append("../tools/")

# Path to access raw data
# data_path = '../data/'
data_path = "/zdata/downloads/ww3_CFSR/"

# Import Libraries
import numpy as np
from netCDF4 import Dataset, num2date
import glob
import datetime

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

# Call Fp and WSP data from WW3 netcdf files
######## fp ########
# initialize masked array:
time_fp = []

# initialize counter:
cn = 0

# Loop through filenames:
for f in filenames_fp:

    # Set nc variable:
    nc_fp = Dataset(f, "r")

    # call peak frequency and time:
    fp_h = nc_fp.variables["fp"][:]
    itime = num2date(nc_fp.variables["time"][:], nc_fp.variables["time"].units)

    # save the hourly wsp array and time:
    if cn == 0:
        fp = np.ma.copy(fp_h)
    elif cn > 0:
        fp = np.ma.concatenate((fp, fp_h), axis=0)
    time_fp.append(itime[0])

    # counter sum
    cn += 1
    print(cn)

######## WSP ########
# initialize time array:
time_wnd = []

# initialize counter:
cn = 0

# Loop through filenames
for f in filenames_wnd:

    # set nc variable :
    nc_wnd = Dataset(f, "r")

    # call wind velocity components and time:
    uwnd = nc_wnd.variables["uwnd"][:]
    vwnd = nc_wnd.variables["vwnd"][:]
    itime = num2date(nc_wnd.variables["time"][:], nc_wnd.variables["time"].units)

    # compute wind speed:
    wsp_h = np.sqrt((uwnd ** 2) + (vwnd ** 2))

    # save the hourly wsp array and time:
    if cn == 0:
        wsp = np.ma.copy(wsp_h)
    elif cn > 0:
        wsp = np.ma.concatenate((wsp, wsp_h), axis=0)
    time_wnd.append(itime[0])

    # counter sum
    cn += 1
    print(cn)

# initialize variables:
g = 9.81

# compute phase speed:
cp = g / (2 * np.pi * fp)

# compute wave age:
wave_age = cp / wsp

# Compute probability of swell:

######## Seasonally ########
# intialize ordered dictionary
prob_s = OrderedDict()
prob_s["p_swell"], prob_s["p_wind"] = [], []

# initialize time indices
months = np.array([m.month for m in time_c])

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
    wave_age_s = wave_age[ind, :, :]

    # Separate Swell and Wind-Sea Events
    swell_ind = wave_age_s > 1.2
    wind_ind = wave_age_s <= 1.2

    # count the number of swell events and wind_sea events
    N_swell = np.sum(swell_ind, axis=0)
    N_wind = np.sum(wind_ind, axis=0)
    N_total = N_swell + N_wind

    # Compute probability of swell and Wind-Sea
    p_swell = N_swell / N_total
    p_wind = N_wind / N_total

    # Save in dictionary:
    prob_s["p_swell"].append(p_swell)
    prob_s["p_wind"].append(p_wind)

# Convert list of 2D arrays to a single 3D array for probability of swell:
prob_swell_s = np.ma.dstack(prob_s["p_swell"])

######## Monthly ########
# intialize ordered dictionary
prob_m = OrderedDict()
prob_m["p_swell"], prob_m["p_wind"] = [], []

# initialize time indices
months = np.array([m.month for m in time_c])

# loop through months in time series:
# initialize variables
month = np.arange(1, 13, 1)

# initialize month loop
for m in range(0, 12):

    # call month:
    imonth = month[m]

    # initialize month indices
    ind = months == imonth

    # call wave age from month
    wave_age_s = wave_age[ind, :, :]

    # Separate Swell and Wind-Sea Events
    swell_ind = wave_age_s > 1.2
    wind_ind = wave_age_s <= 1.2

    # count the number of swell events and wind_sea events
    N_swell = np.sum(swell_ind, axis=0)
    N_wind = np.sum(wind_ind, axis=0)
    N_total = N_swell + N_wind

    # Compute probability of swell and Wind-Sea
    p_swell = N_swell / N_total
    p_wind = N_wind / N_total

    # Save in dictionary:
    prob_m["p_swell"].append(p_swell)
    prob_m["p_wind"].append(p_wind)

# Convert list of 2D arrays to a single 3D array for probability of swell:
prob_swell_m = np.ma.dstack(prob_m["p_swell"])

# Set summary variable:
summary = "Data contained in this netCDF file is derived from the wave hindcast peak wave frequency (fp) and wind speed (WSP) product produced by French Research Institute for Exploitation of the Sea (IFREMER) using the WAVE-height, WATer depth and Current Hindcasting III (WW3) wave model forced by Climate Forecast System Reanalysis (CFSR) winds (ftp://ftp.ifremer.fr/ifremer/ww3/HINDCAST). Thus, this data is an intermediate product. Here, the wave age, defined as the ratio between peak phase speed and wind speed 10 meters above the ocean surface, is calculated across the globe from 66N to 66S over the time period of 1 January 1993 to 31 December 2015 and is used to categorize wave fields as either dominated by swell or wind-sea. For wave age less than or equal to 1.2, the wave field is considered to be dominated by wind-sea. For wave age greater than 1.2, the wave field is considered to be dominated by swell. Probability of swell illustrate the fraction of time the wave field is swell-dominated relative to the total number of wave measurements. Probability of swell is computed seasonally and monthly and are stored in a 3-dimensional (time, latitude, longitude) masked arrays."

# Save seasonal and monthly Probability of Swell data
save_netcdf_prob_swell(
    monthly_prob_swell=prob_swell_m,
    seasonal_prob_swell=prob_swell_s,
    lon=lon_wnd,
    lat=lat_wnd,
    monthly_time=np.arange(1, 13, 1),
    seasonal_time=np.arange(1, 4, 1),
    output="../data/WW3_probability_swell.nc",
    summary=summary,
)
