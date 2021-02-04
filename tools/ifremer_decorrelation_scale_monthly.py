# Ifremer SWH Monthly Decorrelation time scales from 1993 to 2015

# Objectives of Notebook: Compute the decorrelation scale on a monthly basis

# Path to access python functions
import sys

sys.path.append("../tools/")

# Path to access intermediate data
data_path = "../data/ifremer_swh/"

# libraries
import numpy as np
from netCDF4 import Dataset, num2date

# my functions
from data_processing import import_data
from averaging_stats import monthly_average
from lsf import least_square_fit
from decorrelation_scale import decor_scale
from save_netcdf_fields import add_global_atrributes, save_netcdf_decor_scale

# set time and space variables
nt, nlon, nlat = 8400, 360, 133

# call data:
swh, time, lat, lon = import_data("IFREMER_swh", data_path)

# Use monthly average function to partition data and time into monthly segments:
swh_month_dict = monthly_average(np.array(time), swh)

# Initialize monthly partitioned swh and time:
swh_monthly_time = np.ma.array(swh_month_dict["time"])
swh_monthly_data = np.ma.array(swh_month_dict["data"])

# Compute decorrelation time scales
# set variables
ntime = swh_monthly_data.shape[0]
decor = np.zeros((ntime, nlat, nlon))

# Loop over time
for itime in range(0, ntime):
    print(itime)

    # call data
    ts_month = swh_monthly_data[itime]

    # Loop over horizontal space
    for ilat in range(0, nlat):
        for ilon in range(0, nlon):

            # Call data for each grid point:
            ts_grid = ts_month[:, ilat, ilon]

            # Do not preform decorrelation analysis if time series has all masked values.
            if ts_grid.mask.all() == False:

                # Count the number of data point in the time series
                ndata = np.count_nonzero(~np.ma.getmask(ts_grid))

                # Proceed with decorrelation time scale calculation if time series has more than one data point.
                if ndata > 1:

                    # detrend grid point:
                    ts_trend, x_trend = least_square_fit(
                        data=ts_grid, trend="linear", parameters=2, period=None
                    )

                    # remove linear trend:
                    ts_detrend = ts_grid - ts_trend

                    # Compute decorrelation time scale:
                    (
                        ds,
                        ds_N,
                        coef_pos,
                        coef_neg,
                        upper_95_CI,
                        lower_95_CI,
                    ) = decor_scale(
                        data=ts_grid,
                        window=len(ts_grid),
                        lag=len(ts_grid),
                        method="integral_unbiased_coef",
                        dt=1,
                        units=["days", "days"],
                    )

                    # Save decorrelation time scale and autocorrelation function
                    decor[itime, ilat, ilon] = ds

                # Set decorrelation time scale to 1 if only one time step exists in the time series
                elif ndata == 1:

                    # Compute decorrelation time scale:
                    ds = 1

                    # Save decorrelation time scale and autocorrelation function
                    decor[itime, ilat, ilon] = ds

# Save data in a NetCDF file:
# Initialize variables
output = "../data/decor_scales/IFREMER_swh_decor_time_scale.nc"
summary = "Data contained in this netCDF file is derived from the French Research Institute for Exploitation of the Sea (IFREMER) cross-calibrated along-track satellite altimetry significant wave height (SWH) product (ftp://ftp.ifremer.fr/ifremer/cersat/products/swath/altimeters/waves). Thus, this data is an intermediate product. Here, the decorrelation time scales are computed from integrals of the lagged covariance for each month from January 1993 to December 2015 across the globe from 66N to 66S. Decorrelation time scales for Hs are stored in a 3-dimensional (time, latitude, longitude) masked array."

# Save in NetCDF
save_netcdf_decor_scale(decor, lon, lat, swh_monthly_time, output, summary)
