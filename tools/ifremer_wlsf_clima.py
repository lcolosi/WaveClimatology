# Ifremer Weighted Least Squares fit program for monthly climatology data

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

# call data:
swh, time, lat, lon = import_data("IFREMER_swh", data_path)

# Calculate the monthly climatologies of swh from 1993 to 2015
swh_clima_dict = clima_mean(date_time=np.ma.array(time), data=swh)
swh_clima_mean = np.ma.array(swh_clima_dict["mean"])
swh_clima_std = np.ma.array(swh_clima_dict["std"])
swh_clima_n = np.ma.array(swh_clima_dict["N"])

# Mask all standard deviations equal to zero:
swh_clima_std_m = np.ma.masked_equal(swh_clima_std, 0)

# call monthly decorrelation scale
nc = Dataset(data_path + "IFREMER_swh_decor_time_scale.nc", "r")
decor = nc.variables["decor_scale"][:]
time_decor = num2date(nc.variables["time"][:], nc.variables["time"].units)

# Calculate monthly climatologies of monthly decorrelation scales
decor_clima_dict = clima_mean(date_time=time_decor, data=decor)
decor_clima_mean = np.ma.array(decor_clima_dict["mean"])

# Compute standard error of the mean
swh_n_eff = swh_clima_n / decor_clima_mean
swh_clima_stdm = swh_clima_std_m / np.sqrt(swh_n_eff)

# Apply stdm mask to mean:
swh_clima_mean = np.ma.masked_where(np.ma.getmask(swh_clima_stdm), swh_clima_mean)

# initialze variables
parameters = 5  # Fit mean, annual cycle, and semi-annual cycle to data

# initialize masked arrays:
swh_amp1_m = np.ma.masked_all([nlat, nlon])
swh_phase1_m = np.ma.masked_all([nlat, nlon])
swh_amp2_m = np.ma.masked_all([nlat, nlon])
swh_phase2_m = np.ma.masked_all([nlat, nlon])
swh_fve_m = np.ma.masked_all([nlat, nlon])
swh_amp1_sig_p = np.ma.masked_all([nlat, nlon])
swh_phase1_sig_p = np.ma.masked_all([nlat, nlon])
swh_amp2_sig_p = np.ma.masked_all([nlat, nlon])
swh_phase2_sig_p = np.ma.masked_all([nlat, nlon])
swh_count_m = np.ma.masked_all([nlat, nlon])

# Loop through latitude and longitude
for ilat in range(0, nlat):
    for ilon in range(0, nlon):

        # call time series from monthly climatological mean and standard error:
        swh_mean_ts = swh_clima_mean[:, ilat, ilon]
        swh_stdm_ts = swh_clima_stdm[:, ilat, ilon]

        # Find number of non-masked values for each season of the year:
        ival_swh_w = np.count_nonzero(
            ~np.ma.getmask(np.ma.hstack((swh_mean_ts[11], swh_mean_ts[0:2])))
        )
        ival_swh_spr = np.count_nonzero(~np.ma.getmask(swh_mean_ts[2:5]))
        ival_swh_sum = np.count_nonzero(~np.ma.getmask(swh_mean_ts[5:8]))
        ival_swh_f = np.count_nonzero(~np.ma.getmask(swh_mean_ts[8:11]))

        # Count number of observations of mean:
        swh_count = len(swh_mean_ts[np.ma.nonzero(swh_mean_ts)])

        # place lsf calculation condition:
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

# Save lsf parameters in netCDF files:
# Initialize variables
output = "../data/IFREMER_swh_lsf_parameters.nc"
summary = "Data contained in this netCDF file is derived from the French Research Institute for Exploitation of the Sea (IFREMER) cross-calibrated along-track satellite altimetry significant wave height (SWH) product (ftp://ftp.ifremer.fr/ifremer/cersat/products/swath/altimeters/waves). Thus, this data is an intermediate product. Here, the parameters of the weighted least-squares fit to the mean, annual, and semi-annual cycles of SWH monthly climatologies, computed from the IFREMER binned SWH daily data collected from 1 January 1993 to 31 December 2015, along their uncertainties are stored in a 2-dimensional (latitude, longitude) masked array. Parameters include annual cycle amplitude and phase, semi-annual cycle amplitude and phase, and fraction of variance explained by model. Models were only fitted to monthly climatologies that had at least one data point per season and a total of 5 or more data points. Grid points over land or with insufficient data are masked. Uncertainties, computed using error propagation for amplitude and phase, represent the standard error of the mean amplitude and phase. Uncertainty is not computed for fraction of variance explained. Phase is in units of radians. Convert to units of months by multiplying annual phase by 6/pi and semi-annual phase by 3/pi."

# Save to netCDF files:
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
    output,
    summary,
)
