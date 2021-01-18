# WW3 Weighted Least Squares fit program for monthly climatology data

# Path to access python functions
import sys

sys.path.append("../tools/")

# Path to access intermediate data
data_path_ws = "../data/ww3_swh/"
data_path_ww = "../data/ww3_wsp/"
data_path_decor = "../data/decor_scales/"

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

# Call data:
wsp, time, lat, lon = import_data("WW3_wsp", data_path_ww)
swh, time, lat, lon = import_data("WW3_swh", data_path_ws)

# Calculate the monthly averaged from 1993 to 2015
#### SWH ####
swh_clima_dict = clima_mean(date_time=np.ma.array(time), data=swh)
swh_clima_mean = np.ma.array(swh_clima_dict["mean"])
swh_clima_std = np.ma.array(swh_clima_dict["std"])
swh_clima_n = np.ma.array(swh_clima_dict["N"])
#### WSP ####
wsp_clima_dict = clima_mean(date_time=np.ma.array(time), data=wsp)
wsp_clima_mean = np.ma.array(wsp_clima_dict["mean"])
wsp_clima_std = np.ma.array(wsp_clima_dict["std"])
wsp_clima_n = np.ma.array(wsp_clima_dict["N"])

# call monthly decorrelation scale
#### SWH ####
nc_swh = Dataset(data_path_decor + "WW3_swh_decor_time_scale.nc", "r")
decor_swh = nc_swh.variables["decor_scale"][:]
time_decor_swh = num2date(nc_swh.variables["time"][:], nc_swh.variables["time"].units)
#### WSP ####
nc_wsp = Dataset(data_path_decor + "WW3_wsp_decor_time_scale.nc", "r")
decor_wsp = nc_wsp.variables["decor_scale"][:]
time_decor_wsp = num2date(nc_wsp.variables["time"][:], nc_wsp.variables["time"].units)

# Calculate monthly climatologies of monthly decorrelation scales
#### SWH ####
decor_clima_dict = clima_mean(date_time=time_decor_swh, data=decor_swh)
swh_decor_clima_mean = np.ma.array(decor_clima_dict["mean"])
#### WSP ####
decor_clima_dict = clima_mean(date_time=time_decor_wsp, data=decor_wsp)
wsp_decor_clima_mean = np.ma.array(decor_clima_dict["mean"])

# Compute standard error of the mean
swh_clima_n_eff = swh_clima_n / swh_decor_clima_mean
swh_clima_stdm = swh_clima_std / np.sqrt(swh_clima_n_eff)
wsp_clima_n_eff = wsp_clima_n / wsp_decor_clima_mean
wsp_clima_stdm = wsp_clima_std / np.sqrt(wsp_clima_n_eff)

# Apply stdm mask to mean:
swh_clima_mean = np.ma.masked_where(np.ma.getmask(swh_clima_stdm), swh_clima_mean)
wsp_clima_mean = np.ma.masked_where(np.ma.getmask(wsp_clima_stdm), wsp_clima_mean)

# initialze variables
parameters = 5  # Fit mean, annual cycle, and semi-annual cycle to data

# initialize masked arrays:
########### WSP ###########
wsp_whfit_m = np.ma.masked_all([nt, nlat, nlon])
wsp_amp1_m = np.ma.masked_all([nlat, nlon])
wsp_phase1_m = np.ma.masked_all([nlat, nlon])
wsp_amp2_m = np.ma.masked_all([nlat, nlon])
wsp_phase2_m = np.ma.masked_all([nlat, nlon])
wsp_fve_m = np.ma.masked_all([nlat, nlon])
wsp_amp1_sig_p = np.ma.masked_all([nlat, nlon])
wsp_phase1_sig_p = np.ma.masked_all([nlat, nlon])
wsp_amp2_sig_p = np.ma.masked_all([nlat, nlon])
wsp_phase2_sig_p = np.ma.masked_all([nlat, nlon])

########### SWH ###########
swh_whfit_m = np.ma.masked_all([nt, nlat, nlon])
swh_amp1_m = np.ma.masked_all([nlat, nlon])
swh_phase1_m = np.ma.masked_all([nlat, nlon])
swh_amp2_m = np.ma.masked_all([nlat, nlon])
swh_phase2_m = np.ma.masked_all([nlat, nlon])
swh_fve_m = np.ma.masked_all([nlat, nlon])
swh_amp1_sig_p = np.ma.masked_all([nlat, nlon])
swh_phase1_sig_p = np.ma.masked_all([nlat, nlon])
swh_amp2_sig_p = np.ma.masked_all([nlat, nlon])
swh_phase2_sig_p = np.ma.masked_all([nlat, nlon])

# Loop through latitude and longitude
for ilat in range(0, nlat):
    for ilon in range(0, nlon):

        ########### WSP ###########

        # call time series from monthly mean and standard error:
        wsp_mean_ts = wsp_clima_mean[:, ilat, ilon]
        wsp_stdm_ts = wsp_clima_stdm[:, ilat, ilon]

        # Find number of non-masked values for each season of the year:
        ival_wsp_w = np.count_nonzero(
            ~np.ma.getmask(np.ma.hstack((wsp_mean_ts[11], wsp_mean_ts[0:2])))
        )
        ival_wsp_spr = np.count_nonzero(~np.ma.getmask(wsp_mean_ts[2:5]))
        ival_wsp_sum = np.count_nonzero(~np.ma.getmask(wsp_mean_ts[5:8]))
        ival_wsp_f = np.count_nonzero(~np.ma.getmask(wsp_mean_ts[8:11]))

        # Count number of observations of mean and stdm wsp:
        wsp_count = len(wsp_mean_ts[np.ma.nonzero(wsp_mean_ts)])

        # place condition:
        if (
            ival_wsp_w > 0
            and ival_wsp_spr > 0
            and ival_wsp_sum > 0
            and ival_wsp_f > 0
            and wsp_count >= 5
        ):

            # compute least square fit:
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

            # Compute uncertainties of parameters least squares fit:
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

        ########### SWH ###########

        # call time series from monthly mean and standard error:
        swh_mean_ts = swh_clima_mean[:, ilat, ilon]
        swh_stdm_ts = swh_clima_stdm[:, ilat, ilon]

        # Find number of non-masked values for each season of the year:
        ival_swh_w = np.count_nonzero(
            ~np.ma.getmask(np.ma.hstack((swh_mean_ts[11], swh_mean_ts[0:2])))
        )
        ival_swh_spr = np.count_nonzero(~np.ma.getmask(swh_mean_ts[2:5]))
        ival_swh_sum = np.count_nonzero(~np.ma.getmask(swh_mean_ts[5:8]))
        ival_swh_f = np.count_nonzero(~np.ma.getmask(swh_mean_ts[8:11]))

        # Count number of observations of mean swh:
        swh_count = len(swh_mean_ts[np.ma.nonzero(swh_mean_ts)])

        # place condition:
        if (
            ival_swh_w > 0
            and ival_swh_spr > 0
            and ival_swh_sum > 0
            and ival_swh_f > 0
            and swh_count >= 5
        ):

            # compute least square fit:
            swh_hfit_w, x_swh_w, x_swh_sigma = weighted_least_square_fit(
                data=np.ma.copy(swh_mean_ts),
                sigma=np.ma.copy(swh_stdm_ts),
                trend="sinusoidal",
                parameters=parameters,
                period=12,
            )

            # compute parameters of least square fit
            (
                swh_res,
                swh_rms,
                swh_amp1,
                swh_phase1,
                swh_amp2,
                swh_phase2,
                swh_fve,
            ) = LSF_parameters(
                data=np.ma.copy(swh_mean_ts),
                model=swh_hfit_w,
                x_solution=x_swh_w,
                trend="sinusoidal",
                parameters=parameters,
                lsf="weighted",
                sigma=np.ma.copy(swh_stdm_ts),
            )

            # Compute uncertainties of parameters least squares fit:
            (
                sigma_phase_1,
                sigma_amp_1,
                sigma_phase_2,
                sigma_amp_2,
            ) = uncertainty_phase_amp(
                parameters=x_swh_w,
                parameters_sigma=x_swh_sigma,
                cycles="two",
                trend=None,
            )

            # save data into 2D arrays:
            swh_amp1_m[ilat, ilon] = swh_amp1
            swh_phase1_m[ilat, ilon] = swh_phase1
            swh_amp2_m[ilat, ilon] = swh_amp2
            swh_phase2_m[ilat, ilon] = swh_phase2
            swh_fve_m[ilat, ilon] = swh_fve
            swh_amp1_sig_p[ilat, ilon] = sigma_amp_1
            swh_phase1_sig_p[ilat, ilon] = sigma_phase_1
            swh_amp2_sig_p[ilat, ilon] = sigma_amp_2
            swh_phase2_sig_p[ilat, ilon] = sigma_phase_2

# Save WW3 swh and wsp data in netCDF files:
########### SWH ###########
output_ws = "../data/lsf_parameters/WW3_swh_lsf_parameters.nc"
summary_ws = "Data contained in this netCDF file is derived from the wave hindcast significant wave height (SWH) product produced by French Research Institute for Exploitation of the Sea (IFREMER) using the WAVE-height, WATer depth and Current Hindcasting III (WW3) wave model forced by Climate Forecast System Reanalysis (CFSR) winds (ftp://ftp.ifremer.fr/ifremer/ww3/HINDCAST). Thus, this data is an intermediate product. Here, the parameters of the weighted least-squares fit to the mean, annual, and semi-annual cycles of SWH monthly climatologies, computed from the WW3 SWH daily and 1 degree by 1 degree averaged data spanning from 1 January 1993 to 31 December 2015, along their uncertainties are stored in a 2-dimensional (latitude, longitude) masked array. Parameters include annual cycle amplitude and phase, semi-annual cycle amplitude and phase, and fraction of variance explained by model. Models were only fitted to monthly climatologies that had at least one data point per season and a total of 5 or more data points. Grid points over land or with insufficient data are masked. Uncertainties, computed using error propagation for amplitude and phase, represent the standard error of the mean amplitude and phase. Uncertainty is not computed for fraction of variance explained. Phase is in units of radians. Convert to units of months by multiplying annual phase by 6/pi and semi-annual phase by 3/pi."

# Save to netCDF file:
save_netcdf_lsf_parameters(
    swh_amp1_m,
    swh_phase1_m,
    swh_amp2_m,
    swh_phase2_m,
    swh_fve_m,
    swh_amp1_sig_p,
    swh_phase1_sig_p,
    swh_amp2_sig_p,
    swh_phase2_sig_p,
    lon,
    lat,
    output_ws,
    summary_ws,
)


########### WSP ###########
output_ww = "../data/lsf_parameters/WW3_wsp_lsf_parameters.nc"
summary_ww = "Data contained in this netCDF file is derived from the wave hindcast wind speed (WSP) product produced by French Research Institute for Exploitation of the Sea (IFREMER) using the WAVE-height, WATer depth and Current Hindcasting III (WW3) wave model forced by Climate Forecast System Reanalysis (CFSR) winds (ftp://ftp.ifremer.fr/ifremer/ww3/HINDCAST). Thus, this data is an intermediate product. Here, the parameters of the weighted least-squares fit to the mean, annual, and semi-annual cycles of SWH monthly climatologies, computed from the WW3 WSP daily and 1 degree by 1 degree averaged data spanning from 1 January 1993 to 31 December 2015, along their uncertainties are stored in a 2-dimensional (latitude, longitude) masked array. Parameters include annual cycle amplitude and phase, semi-annual cycle amplitude and phase, and fraction of variance explained by model. Models were only fitted to monthly climatologies that had at least one data point per season and a total of 5 or more data points. Grid points over land or with insufficient data are masked. Uncertainties, computed using error propagation for amplitude and phase, represent the standard error of the mean amplitude and phase. Uncertainty is not computed for fraction of variance explained. Phase is in units of radians. Convert to units of months by multiplying annual phase by 6/pi and semi-annual phase by 3/pi. WW3 models WSP 10 meters above the ocean surface."

# Save to netCDF file:
save_netcdf_lsf_parameters(
    wsp_amp1_m,
    wsp_phase1_m,
    wsp_amp2_m,
    wsp_phase2_m,
    swh_fve_m,
    swh_amp1_sig_p,
    swh_phase1_sig_p,
    swh_amp2_sig_p,
    swh_phase2_sig_p,
    lon,
    lat,
    output_ww,
    summary_ww,
)
