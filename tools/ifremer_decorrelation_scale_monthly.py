# Ifremer SWH Monthly Decorrelation time scales from 1993 to 2015

# Objectives of Notebook: Compute the decorrelation scale on a monthly basis

import sys
sys.path.append('../tools/')

#libraries
import numpy as np  
from netCDF4 import Dataset, num2date

#my functions 
from monthly_mean import monthly_average
from decorrelation_ts import decor_scale

#set time and space variables
nt, nlon, nlat = 8400, 360, 133

#Set filename
filename = '../data/ifremer_swh_daily_binned_data_93_16_bia.nc'

#set nc variable: 
nc =  Dataset(filename, 'r')

#call data
swh = nc.variables['swh'][:]
lon = nc.variables['lon'][:]
lat = nc.variables['lat'][:]
time_i = num2date(nc.variables['time'][:], nc.variables['time'].units) 

#resrict data to -66 to 66:
ii = np.where(abs(lat) <= 66)[0]
swh = swh[:,ii,:]
lat = lat[ii]

#find initial and final indices: 
#create year vector: 
years = np.array([y.year for y in time_i])

#create boolean arrays and combine them:
ind_92 = years != 1992
ind_16 = years != 2016
ind_time = ind_92*ind_16

#use the compress function to find all indices that do not lie in 2016 and extract slices of matirx along the time axis from swh
swh_c = np.compress(ind_time, swh, axis = 0)

#extract the time steps: 
time_c = time_i[ind_time]

#Calculate monthly average
swh_month_data = monthly_average(date_time = time_c, data = swh_c, size = '3d')

#For swh:
swh_monthly_time = np.ma.array(swh_month_data['time'])
swh_monthly_data = np.ma.array(swh_month_data['data'])
swh_monthly_mean = np.ma.array(swh_month_data['mean'])
swh_monthly_median = np.ma.array(swh_month_data['median'])
swh_monthly_std = np.ma.array(swh_month_data['std'])
swh_monthly_n = np.ma.array(swh_month_data['N'])

#Compute decorrelation time scales 
#set variables 
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
            
            #Do not preform decorrelation analysis if time series has all masked values. 
            if ts_grid.mask.all() == False:  
            
                
                #Count the number of data point in the time series 
                ndata = np.count_nonzero(~np.ma.getmask(ts_grid))
                        
                #Proceed with decorrelation time scale calculation if time series has more than one data point.
                if ndata > 1: 

                    #Compute decorrelation time scale: 
                    ds, ds_N, coef_pos, coef_neg, upper_95_CI, lower_95_CI = decor_scale(data = ts_grid, window = len(ts_grid), lag = len(ts_grid), method = 'integral_unbiased_coef', dt = 1, units = ['days', 'days'])

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
output = '/zdata/downloads/colosi_data_bk/ucsd_lib_data_repo/IFREMER_swh_decor_time_scale.nc'
summary = 'Data contained in this netCDF file is derived from the French Research Institute for Exploitation of the Sea (IFREMER) cross-calibrated along-track satellite altimetry significant wave height (SWH) product (ftp://ftp.ifremer.fr/ifremer/cersat/products/swath/altimeters/waves). Thus, this data is an intermediate product. Here, the decorrelation time scales are computed from integrals of the lagged covariance for each month from January 1993 to December 2015 across the globe from 66N to 66S. Decorrelation time scales are stored in a 3-dimensional (time, latitude, longitude) masked array.'

#Save in NetCDF
save_netcdf.save_netcdf_decor_scale(decor, lon, lat, swh_monthly_time, output, summary)
                    
#Save data in a npz file: 
np.savez('/zdata/downloads/colosi_data_bk/npz_data/decor_scale/monthly_global_ds_ifremer_swh_int', swh_time_monthly = swh_monthly_time, swh_decor_monthly = decor, swh_autocor = autocor, swh_autocor_mask = np.ma.getmask(autocor))