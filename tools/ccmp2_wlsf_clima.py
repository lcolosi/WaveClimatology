# CCMP2 Weighted Least Squares fit program for monthly climatological data

# Path to access python functions
import sys

sys.path.append("../tools/")

# Path to access intermediate data
data_path = "../data/"

# Import Libraries
import numpy as np
from netCDF4 import Dataset, num2date

# Import functions
from data_processing import import_data
from averaging_stats import clima_mean
from lsf import weighted_least_square_fit, LSF_parameters, uncertainty_phase_amp
from save_netcdf_fields import save_netcdf_lsf_parameters

# Set dimensions for data of space and time
nt, nlon, nlat = 12, 360, 133

# Call wsp data:
wsp, time, lat, lon = import_data("CCMP2_wsp", data_path)

# Calculate monthly climatologies
wsp_clima_dict = clima_mean(date_time=np.ma.array(time), data=wsp)
wsp_clima_mean = np.ma.array(wsp_clima_dict["mean"])
wsp_clima_std = np.ma.array(wsp_clima_dict["std"])
wsp_clima_n = np.ma.array(wsp_clima_dict["N"])

# call monthly decorrelation scale
nc = Dataset(data_path + "CCMP2_wsp_decor_time_scale.nc", "r")
decor = nc.variables["decor_scale"][:]
time_decor = num2date(nc.variables["time"][:], nc.variables["time"].units)

# Calculate monthly climatologies of monthly decorrelation scales
decor_clima_dict = clima_mean(date_time=time_decor, data=decor)
decor_clima_mean = np.ma.array(decor_clima_dict["mean"])

# Compute standard error of the mean
wsp_clima_n_eff = wsp_clima_n / decor_clima_mean
wsp_clima_stdm = wsp_clima_std / np.sqrt(wsp_clima_n_eff)

# Apply stdm mask to mean:
wsp_clima_mean = np.ma.masked_where(np.ma.getmask(wsp_clima_stdm), wsp_clima_mean)

# initialze variables
parameters = 5  # Fit mean, annual cycle, and semi-annual cycle to data

# initialize masked arrays:
wsp_amp1_m = np.ma.masked_all([nlat, nlon])
wsp_phase1_m = np.ma.masked_all([nlat, nlon])
wsp_amp2_m = np.ma.masked_all([nlat, nlon])
wsp_phase2_m = np.ma.masked_all([nlat, nlon])
wsp_fve_m = np.ma.masked_all([nlat, nlon])
wsp_amp1_sig_p = np.ma.masked_all([nlat, nlon])
wsp_phase1_sig_p = np.ma.masked_all([nlat, nlon])
wsp_amp2_sig_p = np.ma.masked_all([nlat, nlon])
wsp_phase2_sig_p = np.ma.masked_all([nlat, nlon])
wsp_count_m = np.ma.masked_all([nlat, nlon])

# Loop through latitude and longitude
for ilat in range(0, nlat):
    for ilon in range(0, nlon):

        # call time series from monthly climatological mean and standard error:
        wsp_mean_ts = wsp_clima_mean[:, ilat, ilon]
        wsp_stdm_ts = wsp_clima_stdm[:, ilat, ilon]

        # Find number of non-masked values for each season of the year:
        ival_wsp_w = np.count_nonzero(
            ~np.ma.getmask(np.ma.hstack((wsp_mean_ts[11], wsp_mean_ts[0:2])))
        )
        ival_wsp_spr = np.count_nonzero(~np.ma.getmask(wsp_mean_ts[2:5]))
        ival_wsp_sum = np.count_nonzero(~np.ma.getmask(wsp_mean_ts[5:8]))
        ival_wsp_f = np.count_nonzero(~np.ma.getmask(wsp_mean_ts[8:11]))

        # Count number of observations of mean and stdm swh:
        wsp_count_m[ilat, ilon] = len(wsp_mean_ts[np.ma.nonzero(wsp_mean_ts)])

        # place lsf calculation condition:
        if (
            ival_wsp_w > 0
            and ival_wsp_spr > 0
            and ival_wsp_sum > 0
            and ival_wsp_f > 0
            and wsp_count_m[ilat, ilon] >= 5
        ):

            # compute weighted least square fit:
            wsp_hfit_w, x_wsp_w, x_wsp_sigma = weighted_least_square_fit(
                data=np.ma.copy(wsp_mean_ts),
                sigma=np.ma.copy(wsp_stdm_ts),
                trend="sinusoidal",
                parameters=parameters,
                period=12,
            )

            # compute parameters of least square fit
            (
                wsp_res,
                wsp_rms,
                wsp_amp1,
                wsp_phase1,
                wsp_amp2,
                wsp_phase2,
                wsp_fve,
            ) = LSF_parameters(
                data=np.ma.copy(wsp_mean_ts),
                model=wsp_hfit_w,
                x_solution=x_wsp_w,
                trend="sinusoidal",
                parameters=parameters,
                lsf="weighted",
                sigma=np.ma.copy(wsp_stdm_ts),
            )

            # Compute uncertainties of least squares fit parameters:
            (
                sigma_phase_1,
                sigma_amp_1,
                sigma_phase_2,
                sigma_amp_2,
            ) = uncertainty_phase_amp(
                parameters=x_wsp_w,
                parameters_sigma=x_wsp_sigma,
                cycles="two",
                trend=None,
            )

            # save data into 2D arrays:
            wsp_amp1_m[ilat, ilon] = wsp_amp1
            wsp_phase1_m[ilat, ilon] = wsp_phase1
            wsp_amp2_m[ilat, ilon] = wsp_amp2
            wsp_phase2_m[ilat, ilon] = wsp_phase2
            wsp_fve_m[ilat, ilon] = wsp_fve
            wsp_amp1_sig_p[ilat, ilon] = sigma_amp_1
            wsp_phase1_sig_p[ilat, ilon] = sigma_phase_1
            wsp_amp2_sig_p[ilat, ilon] = sigma_amp_2
            wsp_phase2_sig_p[ilat, ilon] = sigma_phase_2

# Save lsf parameters in netCDF files:
# Initialize variables
output = "../data/CCMP2_wsp_lsf_parameters.nc"
summary = "Data contained in this netCDF file is derived from the Cross Calibrated Multi-Platform version 2 (CCMP2) wind vector analysis produced by Remote Sensing Systems product (available at www.remss.com). Thus, this data is an intermediate product. Here, the parameters of the weighted least-squares fit to the mean, annual, and semi-annual cycles of wind speed (WSP) monthly climatologies, computed from CCMP2 WSP daily and 1 degree by 1 degree averaged data spanning from 1 January 1993 to 31 December 2015, along their uncertainties are stored in a 2-dimensional (latitude, longitude) masked array. Parameters include annual cycle amplitude and phase, semi-annual cycle amplitude and phase, and fraction of variance explained by model. Models were only fitted to monthly climatologies that had at least one data point per season and a total of 5 or more data points. Grid points over land or with insufficient data are masked. Uncertainties, computed using error propagation for amplitude and phase, represent the standard error of the mean amplitude and phase. Uncertainty is not computed for fraction of variance explained. Phase is in units of radians. Convert to units of months by multiplying annual phase by 6/pi and semi-annual phase by 3/pi. CCMP2 measures WSP 10 meters above the ocean surface."

# Save in NetCDF
save_netcdf_lsf_parameters(
    wsp_amp1_m,
    wsp_phase1_m,
    wsp_amp2_m,
    wsp_phase2_m,
    wsp_fve_m,
    wsp_amp1_sig_p,
    wsp_phase1_sig_p,
    wsp_amp2_sig_p,
    wsp_phase2_sig_p,
    lon,
    lat,
    output,
    summary,
)
