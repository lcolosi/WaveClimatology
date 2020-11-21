# CCMP Version 2 WSP Monthly Decorrelation time scales from 1993 to 2015

# Objectives of Notebook: Compute the decorrelation scale on a monthly basis

import sys
sys.path.append('../tools/')

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
filename = '../data/ccmp_v2_wsp_daily_data_93_15.nc'

#set nc variable in order to read attributes and obtained data: 
nc = Dataset(filename, 'r')

#call data
lon = nc.variables['lon'][:]
lat = nc.variables['lat'][:]
time = num2date(nc.variables['time'][:], nc.variables['time'].units)
wsp = nc.variables['wsp'][:]

#compute monthly averages 
wsp_month_data = monthly_average(date_time = time, data = wsp, size = '3d')

#For wsp:
wsp_monthly_time = np.ma.array(wsp_month_data['time'])
wsp_monthly_data = np.ma.array(wsp_month_data['data'])
wsp_monthly_mean = np.ma.array(wsp_month_data['mean'])
wsp_monthly_median = np.ma.array(wsp_month_data['median'])
wsp_monthly_std = np.ma.array(wsp_month_data['std'])
wsp_monthly_n = np.ma.array(wsp_month_data['N'])

#Compute decorrelation time scales 
# set variables 
ntime = wsp_monthly_data.shape[0]
nmonth = [len(wsp_monthly_data[imonth]) for imonth in range(ntime)]
decor = np.zeros((ntime,nlat,nlon))
autocor = np.ma.masked_all([np.max(nmonth)*2,nlat, nlon])

#Loop over time 
for itime in range(0,ntime): 
    print(itime)
            
    #call data 
    ts_month = wsp_monthly_data[itime]
    
    #Loop over horizontal space
    for ilat in range(0,nlat):
        for ilon in range(0,nlon):
            
            #Call data for each grid point: 
            ts_grid = ts_month[:,ilat,ilon]
            
            #Do not preform decorrelation analysis if time series has all masked values. 
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
output = '/zdata/downloads/colosi_data_bk/ucsd_lib_data_repo/CCMP2_wsp_decor_time_scale.nc'
summary = 'Data contained in this netCDF file is derived from the Cross Calibrated Multi-Platform version 2 (CCMP2) wind vector analysis produced by Remote Sensing Systems product (available at www.remss.com). Thus, this data is an intermediate product. Here, the decorrelation time scales are computed from integrals of the lagged covariance for each month from January 1993 to December 2015 across the globe from 66N to 66S. Decorrelation time scales are stored in a 3-dimensional (time, latitude, longitude) masked array.'

#Save in NetCDF
save_netcdf.save_netcdf_decor_scale(decor, lon, lat, wsp_monthly_time, output, summary)

# Save data in a npz file: 
np.savez('/zdata/downloads/colosi_data_bk/npz_data/decor_scale/monthly_global_ds_ccmp2_wsp_int', wsp_time_monthly = wsp_monthly_time, wsp_decor_monthly = decor, wsp_autocor = autocor, wsp_autocor_mask = np.ma.getmask(autocor))
