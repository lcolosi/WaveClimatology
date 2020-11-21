# WW3 Hs Monthly Decorrelation time scales from 1993 to 2015

# Objectives of Notebook: Compute the decorrelation scale on a monthly basis

import sys
sys.path.append('/zdata/home/lcolosi/python_functions/')

#libraries
import numpy as np
from netCDF4 import Dataset, num2date

#my functions 
from monthly_mean import monthly_average
from decorrelation_ts import decor_scale
from unweighted_least_square_fit import least_square_fit 

#set time and space variables 
nt, nlon, nlat = 8400, 360, 133

#set filename
filename = '/zdata/downloads/colosi_data_bk/binned_data/WW3/CFSR/lc_binned_data/ww3_hs_daily_deresolved_data_93_15.nc'

#set nc variable in order to read attributes and obtained data: 
nc = Dataset(filename, 'r')

#call data
lon = nc.variables['lon'][:]
lat = nc.variables['lat'][:]
time = num2date(nc.variables['time'][:], nc.variables['time'].units)
swh = nc.variables['swh'][:]

#compute monthly averages 
swh_month_dict = monthly_average(date_time = time, data = swh, size = '3d')

swh_monthly_time = np.ma.array(swh_month_dict['time'])
swh_monthly_data = np.ma.array(swh_month_dict['data'])
swh_monthly_mean = np.ma.array(swh_month_dict['mean'])
swh_monthly_median = np.ma.array(swh_month_dict['median'])
swh_monthly_std = np.ma.array(swh_month_dict['std'])
swh_monthly_n = np.ma.array(swh_month_dict['N'])

#Compute decorrelation time scales 
# set variables 
ntime = swh_monthly_data.shape[0]
nmonth = [len(swh_monthly_data[imonth]) for imonth in range(ntime)]
decor = np.zeros((ntime,nlat,nlon))
autocor = np.ma.masked_all([np.max(nmonth)*2,nlat, nlon])

#Loop over time 
for itime in range(0,ntime): 
    print(itime)
            
    #call data 
    ts_month = swh_monthly_data[itime]
    
    #Loop over horizontal space
    for ilat in range(0,nlat):
        for ilon in range(0,nlon):
            
            #Call data for each grid point: 
            ts_grid = ts_month[:,ilat,ilon]
            
            #Do not preform decorrelation analysis if time series has all masked values: 
            if ts_grid.mask.all() == False:  
                
                #Count the number of data point in the time series 
                ndata = np.count_nonzero(~np.ma.getmask(ts_grid))
                        
                #Proceed with decorrelation time scale calculation if time series has more than one data point.
                if ndata > 1:
            
                    #Preform any necessary detrending:
                    #detrend grid point: 
                    ts_trend, x_trend = least_square_fit(data = ts_grid, trend = 'linear', parameters = 2, period = 12)

                    #remove linear trend: 
                    ts_detrend = ts_grid - ts_trend 

                    #Compute decorrelation time scale: 
                    ds, ds_N, coef_pos, coef_neg, upper_95_CI, lower_95_CI = decor_scale(data = ts_detrend, window = len(ts_detrend), lag = len(ts_detrend), method = 'integral_unbiased_coef', dt = 1, units = ['days', 'days'])

                    #Save decorrelation time scale and autocorrelation function
                    decor[itime,ilat,ilon] = ds
                    autocor_l = len(np.hstack((coef_neg, coef_pos)))
                    autocor[:autocor_l,ilat,ilon] = np.hstack((coef_neg, coef_pos))
                
                #Set decorrelation time scale to 1 if only one time step exists in the time series
                elif ndata == 1: 
                    
                    #Compute decorrelation time scale: 
                    ds = 1
                    
                    #Save decorrelation time scale and autocorrelation function
                    decor[itime,ilat,ilon] = ds
                    autocor[:2,ilat,ilon] = np.ma.array((1,0))
                    
#Save data in a NetCDF file: 
#Initialize variables
output = '/zdata/downloads/colosi_data_bk/ucsd_lib_data_repo/WW3_swh_decor_time_scale.nc'
summary = 'Data contained in this netCDF file is derived from the wave hindcast significant wave height (SWH) product produced by French Research Institute for Exploitation of the Sea (IFREMER) using the WAVE-height, WATer depth and Current Hindcasting III (WW3) wave model forced by Climate Forecast System Reanalysis (CFSR) winds (ftp://ftp.ifremer.fr/ifremer/ww3/HINDCAST). Thus, this data is an intermediate product. Here, the decorrelation time scales are computed from integrals of the lagged covariance for each month from January 1993 to December 2015 across the globe from 66N to 66S. Decorrelation time scales are stored in a 3-dimensional (time, latitude, longitude) masked array.'

#Save in NetCDF
save_netcdf.save_netcdf_decor_scale(decor, lon, lat, swh_monthly_time, output, summary)

# Save data in a npz file: 
np.savez('/zdata/downloads/colosi_data_bk/npz_data/decor_scale/monthly_global_ds_ww3_hs_int', swh_time_monthly = swh_monthly_time, swh_decor_monthly = decor, swh_autocor = autocor, swh_autocor_mask = np.ma.getmask(autocor))
