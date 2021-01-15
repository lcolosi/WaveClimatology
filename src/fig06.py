"""
Figure 6 Caption: 
(left column) Northern Hemisphere wind speed in SWARs, averaged over June, July,
and August. (right column) IFREMER SWH (solid blue) and CCMP2 WSP (solid red) 
climatologies extracted from the outlined 4$^{\circ}$ by 4$^{\circ}$ boxes within
SWARs. Blue shading represents the standard error of the mean, dotted blue is 
the annual cycle weighted least-squares fitted to monthly climatology for mean 
SWH of the hemisphere ocean basin the SWAR is located in, and black solid is the 
residual between SWH regional climatology and annual cycle. SWARs include 
Northern California (A and B), Southern Caribbean Sea (C and D), and North Africa
near the coast of Morocco and western Sahara (E and F). (Comparison plots showing
equivalent quantities for WW3 and IFREMER SWH and CFSR and CCMP2 WSP can be found
in Figure~S9A,C,E of the supplementary material).
"""

# Path to access python functions
import sys

sys.path.append("../tools/")

# Path to access intermediate data
data_path_i = "../data/ifremer_swh/"
data_path_c = "../data/ccmp2_wsp/"
data_path_decor = "../data/decor_scales/"

# libraries
import numpy as np
from netCDF4 import Dataset, num2date
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cmocean.cm as cmo
import matplotlib.patches as mpatches

# my functions
from data_processing import import_data
from averaging_stats import clima_mean, stat_moments_temporal
from lsf import weighted_least_square_fit, LSF_parameters
from regional_clima_figs import regional_clima, regional_clima_plot
import cartopy_figs as cart

# call IFREMER SWH and CCMP2 WSP processed data:
swh, time_s, lat_s, lon_s = import_data("IFREMER_swh", data_path_i)
wsp, time_w, lat_w, lon_w = import_data("CCMP2_wsp", data_path_c)

# Call decorrelation time scales
###### SWH ######
nc_swh = Dataset(data_path_decor + "IFREMER_swh_decor_time_scale.nc", "r")
decor_swh = nc_swh.variables["decor_scale"][:]
time_decor_swh = num2date(nc_swh.variables["time"][:], nc_swh.variables["time"].units)
###### WSP ######
nc_wsp = Dataset(data_path_decor + "CCMP2_wsp_decor_time_scale.nc", "r")
decor_wsp = nc_wsp.variables["decor_scale"][:]
time_decor_wsp = num2date(nc_wsp.variables["time"][:], nc_wsp.variables["time"].units)

# Compute WSP statistical moments seasonally
wsp_stats_s = stat_moments_temporal(wsp, time_w, "seasonally", "sample")
wsp_moments_mean = np.ma.array(wsp_stats_s["mean"])

# Calculate monthly climatologies for SWH and WSP data and decorrelation scales:
###### SWH ######
swh_clima_dict = clima_mean(date_time=np.ma.array(time_s), data=swh)
swh_mean_i = np.ma.array(swh_clima_dict["mean"])
swh_var_i = np.ma.array(swh_clima_dict["var"])
swh_n_i = np.ma.array(swh_clima_dict["N"])
decor_clima_dict_s = clima_mean(date_time=time_decor_swh, data=decor_swh)
decor_mean_swh_i = np.ma.array(decor_clima_dict_s["mean"])

###### WSP ######
wsp_clima_dict = clima_mean(date_time=np.ma.array(time_w), data=wsp)
wsp_mean_c = np.ma.array(wsp_clima_dict["mean"])
wsp_var_c = np.ma.array(wsp_clima_dict["var"])
wsp_n_c = np.ma.array(wsp_clima_dict["N"])
decor_clima_dict_w = clima_mean(date_time=time_decor_wsp, data=decor_wsp)
decor_mean_wsp_c = np.ma.array(decor_clima_dict_w["mean"])

## Compute basin-scale annual cycle fit ##
# Set partition indices for the world oceans excluding the equatorial regions and marginal seas:
ind_nph_lat, ind_nph_lon = [116, 133], [106, 266]
ind_npl_lat, ind_npl_lon = [75, 117], [106, 275]
ind_nawh_lat, ind_nawh_lon = [116, 133], [295, 361]
ind_nawl_lat, ind_nawl_lon = [77, 117], [274, 361]
ind_nae_lat, ind_nae_lon = [77, 133], [0, 24]
ind_sp_lat, ind_sp_lon = [0, 60], [131, 292]
ind_saw_lat, ind_saw_lon = [0, 68], [291, 360]
ind_sae_lat, ind_sae_lon = [0, 68], [0, 24]
ind_nioh_lat, ind_nioh_lon = [89, 95], [60, 105]
ind_niol_lat, ind_niol_lon = [77, 90], [51, 105]
ind_mioh_lat, ind_mioh_lon = [69, 78], [40, 98]
ind_miol_lat, ind_miol_lon = [59, 70], [59, 70]
ind_sio_lat, ind_sio_lon = [0, 60], [23, 132]

# Partition mean, variance, number of observations, and decorrelation scales
####### Mean #######
swh_nph = swh_mean_i[
    :, ind_nph_lat[0] : ind_nph_lat[1], ind_nph_lon[0] : ind_nph_lon[1]
]
swh_npl = swh_mean_i[
    :, ind_npl_lat[0] : ind_npl_lat[1], ind_npl_lon[0] : ind_nph_lon[1]
]
swh_nawh = swh_mean_i[
    :, ind_nawh_lat[0] : ind_nawh_lat[1], ind_nawh_lon[0] : ind_nawh_lon[1]
]
swh_nawl = swh_mean_i[
    :, ind_nawl_lat[0] : ind_nawl_lat[1], ind_nawl_lon[0] : ind_nawl_lon[1]
]
swh_nae = swh_mean_i[
    :, ind_nae_lat[0] : ind_nae_lat[1], ind_nae_lon[0] : ind_nae_lon[1]
]
swh_sp = swh_mean_i[:, ind_sp_lat[0] : ind_sp_lat[1], ind_sp_lon[0] : ind_sp_lon[1]]
swh_saw = swh_mean_i[
    :, ind_saw_lat[0] : ind_saw_lat[1], ind_saw_lon[0] : ind_saw_lon[1]
]
swh_sae = swh_mean_i[
    :, ind_sae_lat[0] : ind_sae_lat[1], ind_sae_lon[0] : ind_sae_lon[1]
]
swh_nioh = swh_mean_i[
    :, ind_nioh_lat[0] : ind_nioh_lat[1], ind_nioh_lon[0] : ind_nioh_lon[1]
]
swh_niol = swh_mean_i[
    :, ind_niol_lat[0] : ind_niol_lat[1], ind_niol_lon[0] : ind_niol_lon[1]
]
swh_mioh = swh_mean_i[
    :, ind_mioh_lat[0] : ind_mioh_lat[1], ind_mioh_lon[0] : ind_mioh_lon[1]
]
swh_miol = swh_mean_i[
    :, ind_miol_lat[0] : ind_miol_lat[1], ind_miol_lon[0] : ind_miol_lon[1]
]
swh_sio = swh_mean_i[
    :, ind_sio_lat[0] : ind_sio_lat[1], ind_sio_lon[0] : ind_sio_lon[1]
]
swh_reg = [
    swh_nph,
    swh_npl,
    swh_nawh,
    swh_nawl,
    swh_nae,
    swh_sp,
    swh_saw,
    swh_sae,
    swh_nioh,
    swh_niol,
    swh_mioh,
    swh_miol,
    swh_sio,
]

####### Variance #######
swh_sig_nph = swh_var_i[
    :, ind_nph_lat[0] : ind_nph_lat[1], ind_nph_lon[0] : ind_nph_lon[1]
]
swh_sig_npl = swh_var_i[
    :, ind_npl_lat[0] : ind_npl_lat[1], ind_npl_lon[0] : ind_nph_lon[1]
]
swh_sig_nawh = swh_var_i[
    :, ind_nawh_lat[0] : ind_nawh_lat[1], ind_nawh_lon[0] : ind_nawh_lon[1]
]
swh_sig_nawl = swh_var_i[
    :, ind_nawl_lat[0] : ind_nawl_lat[1], ind_nawl_lon[0] : ind_nawl_lon[1]
]
swh_sig_nae = swh_var_i[
    :, ind_nae_lat[0] : ind_nae_lat[1], ind_nae_lon[0] : ind_nae_lon[1]
]
swh_sig_sp = swh_var_i[:, ind_sp_lat[0] : ind_sp_lat[1], ind_sp_lon[0] : ind_sp_lon[1]]
swh_sig_saw = swh_var_i[
    :, ind_saw_lat[0] : ind_saw_lat[1], ind_saw_lon[0] : ind_saw_lon[1]
]
swh_sig_sae = swh_var_i[
    :, ind_sae_lat[0] : ind_sae_lat[1], ind_sae_lon[0] : ind_sae_lon[1]
]
swh_sig_nioh = swh_var_i[
    :, ind_nioh_lat[0] : ind_nioh_lat[1], ind_nioh_lon[0] : ind_nioh_lon[1]
]
swh_sig_niol = swh_var_i[
    :, ind_niol_lat[0] : ind_niol_lat[1], ind_niol_lon[0] : ind_niol_lon[1]
]
swh_sig_mioh = swh_var_i[
    :, ind_mioh_lat[0] : ind_mioh_lat[1], ind_mioh_lon[0] : ind_mioh_lon[1]
]
swh_sig_miol = swh_var_i[
    :, ind_miol_lat[0] : ind_miol_lat[1], ind_miol_lon[0] : ind_miol_lon[1]
]
swh_sig_sio = swh_var_i[
    :, ind_sio_lat[0] : ind_sio_lat[1], ind_sio_lon[0] : ind_sio_lon[1]
]
swh_reg_sig = [
    swh_sig_nph,
    swh_sig_npl,
    swh_sig_nawh,
    swh_sig_nawl,
    swh_sig_nae,
    swh_sig_sp,
    swh_sig_saw,
    swh_sig_sae,
    swh_sig_nioh,
    swh_sig_niol,
    swh_sig_mioh,
    swh_sig_miol,
    swh_sig_sio,
]

####### Number of data points #######
swh_n_nph = swh_n_i[:, ind_nph_lat[0] : ind_nph_lat[1], ind_nph_lon[0] : ind_nph_lon[1]]
swh_n_npl = swh_n_i[:, ind_npl_lat[0] : ind_npl_lat[1], ind_npl_lon[0] : ind_nph_lon[1]]
swh_n_nawh = swh_n_i[
    :, ind_nawh_lat[0] : ind_nawh_lat[1], ind_nawh_lon[0] : ind_nawh_lon[1]
]
swh_n_nawl = swh_n_i[
    :, ind_nawl_lat[0] : ind_nawl_lat[1], ind_nawl_lon[0] : ind_nawl_lon[1]
]
swh_n_nae = swh_n_i[:, ind_nae_lat[0] : ind_nae_lat[1], ind_nae_lon[0] : ind_nae_lon[1]]
swh_n_sp = swh_n_i[:, ind_sp_lat[0] : ind_sp_lat[1], ind_sp_lon[0] : ind_sp_lon[1]]
swh_n_saw = swh_n_i[:, ind_saw_lat[0] : ind_saw_lat[1], ind_saw_lon[0] : ind_saw_lon[1]]
swh_n_sae = swh_n_i[:, ind_sae_lat[0] : ind_sae_lat[1], ind_sae_lon[0] : ind_sae_lon[1]]
swh_n_nioh = swh_n_i[
    :, ind_nioh_lat[0] : ind_nioh_lat[1], ind_nioh_lon[0] : ind_nioh_lon[1]
]
swh_n_niol = swh_n_i[
    :, ind_niol_lat[0] : ind_niol_lat[1], ind_niol_lon[0] : ind_niol_lon[1]
]
swh_n_mioh = swh_n_i[
    :, ind_mioh_lat[0] : ind_mioh_lat[1], ind_mioh_lon[0] : ind_mioh_lon[1]
]
swh_n_miol = swh_n_i[
    :, ind_miol_lat[0] : ind_miol_lat[1], ind_miol_lon[0] : ind_miol_lon[1]
]
swh_n_sio = swh_n_i[:, ind_sio_lat[0] : ind_sio_lat[1], ind_sio_lon[0] : ind_sio_lon[1]]
swh_reg_n = [
    swh_nph,
    swh_n_npl,
    swh_n_nawh,
    swh_n_nawl,
    swh_n_nae,
    swh_n_sp,
    swh_n_saw,
    swh_n_sae,
    swh_n_nioh,
    swh_n_niol,
    swh_n_mioh,
    swh_n_miol,
    swh_n_sio,
]

####### Decorrelation time scale #######
decor_mean_swh_nph = decor_mean_swh_i[
    :, ind_nph_lat[0] : ind_nph_lat[1], ind_nph_lon[0] : ind_nph_lon[1]
]
decor_mean_swh_npl = decor_mean_swh_i[
    :, ind_npl_lat[0] : ind_npl_lat[1], ind_npl_lon[0] : ind_nph_lon[1]
]
decor_mean_swh_nawh = decor_mean_swh_i[
    :, ind_nawh_lat[0] : ind_nawh_lat[1], ind_nawh_lon[0] : ind_nawh_lon[1]
]
decor_mean_swh_nawl = decor_mean_swh_i[
    :, ind_nawl_lat[0] : ind_nawl_lat[1], ind_nawl_lon[0] : ind_nawl_lon[1]
]
decor_mean_swh_nae = decor_mean_swh_i[
    :, ind_nae_lat[0] : ind_nae_lat[1], ind_nae_lon[0] : ind_nae_lon[1]
]
decor_mean_swh_sp = decor_mean_swh_i[
    :, ind_sp_lat[0] : ind_sp_lat[1], ind_sp_lon[0] : ind_sp_lon[1]
]
decor_mean_swh_saw = decor_mean_swh_i[
    :, ind_saw_lat[0] : ind_saw_lat[1], ind_saw_lon[0] : ind_saw_lon[1]
]
decor_mean_swh_sae = decor_mean_swh_i[
    :, ind_sae_lat[0] : ind_sae_lat[1], ind_sae_lon[0] : ind_sae_lon[1]
]
decor_mean_swh_nioh = decor_mean_swh_i[
    :, ind_nioh_lat[0] : ind_nioh_lat[1], ind_nioh_lon[0] : ind_nioh_lon[1]
]
decor_mean_swh_niol = decor_mean_swh_i[
    :, ind_niol_lat[0] : ind_niol_lat[1], ind_niol_lon[0] : ind_niol_lon[1]
]
decor_mean_swh_mioh = decor_mean_swh_i[
    :, ind_mioh_lat[0] : ind_mioh_lat[1], ind_mioh_lon[0] : ind_mioh_lon[1]
]
decor_mean_swh_miol = decor_mean_swh_i[
    :, ind_miol_lat[0] : ind_miol_lat[1], ind_miol_lon[0] : ind_miol_lon[1]
]
decor_mean_swh_sio = decor_mean_swh_i[
    :, ind_sio_lat[0] : ind_sio_lat[1], ind_sio_lon[0] : ind_sio_lon[1]
]
swh_reg_dcor = [
    decor_mean_swh_nph,
    decor_mean_swh_npl,
    decor_mean_swh_nawh,
    decor_mean_swh_nawl,
    decor_mean_swh_nae,
    decor_mean_swh_sp,
    decor_mean_swh_saw,
    decor_mean_swh_sae,
    decor_mean_swh_nioh,
    decor_mean_swh_niol,
    decor_mean_swh_mioh,
    decor_mean_swh_miol,
    decor_mean_swh_sio,
]

# For each partition:
# Compute Mean
swh_reg_m = [np.ma.mean(np.ma.mean(ireg, axis=1), axis=1) for ireg in swh_reg]
swh_reg_m_n = [
    np.ma.mean(np.ma.array([swh_reg_m[0], swh_reg_m[1]]), axis=0),
    np.ma.mean(np.ma.array([swh_reg_m[2], swh_reg_m[3], swh_reg_m[4]]), axis=0),
    swh_reg_m[5],
    np.ma.mean(np.ma.array([swh_reg_m[6], swh_reg_m[7]]), axis=0),
    np.ma.mean(
        np.ma.array(
            [swh_reg_m[8], swh_reg_m[9], swh_reg_m[10], swh_reg_m[11], swh_reg_m[12]]
        ),
        axis=0,
    ),
]

# Compute average Variance
swh_reg_var_m = [np.ma.mean(np.ma.mean(ireg, axis=1), axis=1) for ireg in swh_reg_sig]
swh_reg_var_m_n = [
    np.ma.mean(np.ma.array([swh_reg_var_m[0], swh_reg_var_m[1]]), axis=0),
    np.ma.mean(
        np.ma.array([swh_reg_var_m[2], swh_reg_var_m[3], swh_reg_var_m[4]]), axis=0
    ),
    swh_reg_var_m[5],
    np.ma.mean(np.ma.array([swh_reg_var_m[6], swh_reg_var_m[7]]), axis=0),
    np.ma.mean(
        np.ma.array(
            [
                swh_reg_var_m[8],
                swh_reg_var_m[9],
                swh_reg_var_m[10],
                swh_reg_var_m[11],
                swh_reg_var_m[12],
            ]
        ),
        axis=0,
    ),
]

# compute the average of the decorrelation time scale and number of data points for each region:
swh_reg_dcor_m = [np.ma.mean(np.ma.mean(ireg, axis=1), axis=1) for ireg in swh_reg_dcor]
swh_reg_dcor_m_n = [
    np.ma.mean(np.ma.array([swh_reg_dcor_m[0], swh_reg_dcor_m[1]]), axis=0),
    np.ma.mean(
        np.ma.array([swh_reg_dcor_m[2], swh_reg_dcor_m[3], swh_reg_dcor_m[4]]), axis=0
    ),
    swh_reg_dcor_m[5],
    np.ma.mean(np.ma.array([swh_reg_dcor_m[6], swh_reg_dcor_m[7]]), axis=0),
    np.ma.mean(
        np.ma.array(
            [
                swh_reg_dcor_m[8],
                swh_reg_dcor_m[9],
                swh_reg_dcor_m[10],
                swh_reg_dcor_m[11],
                swh_reg_dcor_m[12],
            ]
        ),
        axis=0,
    ),
]

# compute the average of number of data points for each region:
swh_reg_n_m = [np.ma.mean(np.ma.mean(ireg, axis=1), axis=1) for ireg in swh_reg_n]
swh_reg_n_m_n = [
    np.ma.mean(np.ma.array([swh_reg_n_m[0], swh_reg_n_m[1]]), axis=0),
    np.ma.mean(np.ma.array([swh_reg_n_m[2], swh_reg_n_m[3], swh_reg_n_m[4]]), axis=0),
    swh_reg_n_m[5],
    np.ma.mean(np.ma.array([swh_reg_n_m[6], swh_reg_n_m[7]]), axis=0),
    np.ma.mean(
        np.ma.array(
            [
                swh_reg_n_m[8],
                swh_reg_n_m[9],
                swh_reg_n_m[10],
                swh_reg_n_m[11],
                swh_reg_n_m[12],
            ]
        ),
        axis=0,
    ),
]

# Compute degrees of freedom or N_effective for stdm calculation:
n_eff_reg = []
for ireg in range(5):
    n_eff_reg.append(swh_reg_n_m_n[ireg] / swh_reg_dcor_m_n[0])

# compute the standard error of the mean
swh_reg_stdm = []
for ireg in range(5):
    swh_reg_stdm.append(np.sqrt(swh_reg_var_m_n[ireg]) / np.sqrt(n_eff_reg[ireg]))

# Perform a weighted least squares annual cycle fit to the mean regional data
hfit_basin_swh, x_data_basin, x_data_sigma_basin = [], [], []
for ireg in range(5):

    # Perform weigthed least squares fit
    hfit, x_data, x_data_sigma = weighted_least_square_fit(
        data=np.ma.copy(swh_reg_m_n[ireg]),
        sigma=np.ma.copy(swh_reg_stdm[ireg]),
        trend="sinusoidal",
        parameters=3,
        period=12,
    )

    # save data
    hfit_basin_swh.append(hfit)
    x_data_basin.append(x_data)
    x_data_sigma_basin.append(x_data_sigma)

# Compute regional climatologies:

####### Northern California #######
# intialize variable
ngrid = 4
lon_grid = 232
lat_grid = 103
location = ["NH", "west"]

# compute SWH regional climatology
(
    swh_mean_nc,
    swh_stdm_nc,
    swh_hfit,
    swh_x_data,
    swh_residual,
    grid_cor_nc,
) = regional_clima(
    swh_mean_i,
    swh_var_i,
    swh_n_i,
    lon_s,
    lat_s,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_swh_i,
    location,
    lsf="weighted",
    parameters=5,
)

# compute WSP regional climatology
(
    wsp_mean_nc,
    wsp_stdm_nc,
    wsp_hfit,
    wsp_x_data,
    wsp_residual,
    grid_cor_nc,
) = regional_clima(
    wsp_mean_c,
    wsp_var_c,
    wsp_n_c,
    lon_w,
    lat_w,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_wsp_c,
    location,
    lsf="weighted",
    parameters=5,
)

# Perform weigthed least squares fit for basin-scale climatology projected onto regional climatology
hfit_reg_nc, x_data_reg, x_data_sigma_reg = weighted_least_square_fit(
    swh_mean_nc,
    swh_stdm_nc,
    trend="sinusoidal",
    parameters=2,
    period=12,
    phase=np.arctan2(x_data_basin[0][1], x_data_basin[0][2]),
)

####### Southern Caribbean #######
# intialize variable
ngrid = 4
lon_grid = 284
lat_grid = 80
location = ["NH", "west"]

# compute SWH regional climatology
(
    swh_mean_sc,
    swh_stdm_sc,
    swh_hfit,
    swh_x_data,
    swh_residual,
    grid_cor_sc,
) = regional_clima(
    swh_mean_i,
    swh_var_i,
    swh_n_i,
    lon_s,
    lat_s,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_swh_i,
    location,
    lsf="weighted",
    parameters=5,
)

# compute WSP regional climatology
(
    wsp_mean_sc,
    wsp_stdm_sc,
    wsp_hfit,
    wsp_x_data,
    wsp_residual,
    grid_cor_sc,
) = regional_clima(
    wsp_mean_c,
    wsp_var_c,
    wsp_n_c,
    lon_w,
    lat_w,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_wsp_c,
    location,
    lsf="weighted",
    parameters=5,
)

# Perform weigthed least squares fit for basin-scale climatology projected onto regional climatology
hfit_reg_sc, x_data_reg, x_data_sigma_reg = weighted_least_square_fit(
    swh_mean_sc,
    swh_stdm_sc,
    trend="sinusoidal",
    parameters=2,
    period=12,
    phase=np.arctan2(x_data_basin[1][1], x_data_basin[1][2]),
)

####### North Africa (Morocco) #######
# intialize variable
ngrid = 4
lon_grid = 345
lat_grid = 96
location = ["NH", "west"]

# compute SWH regional climatology
(
    swh_mean_na,
    swh_stdm_na,
    swh_hfit,
    swh_x_data,
    swh_residual,
    grid_cor_na,
) = regional_clima(
    swh_mean_i,
    swh_var_i,
    swh_n_i,
    lon_s,
    lat_s,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_swh_i,
    location,
    lsf="weighted",
    parameters=5,
)

# compute WSP regional climatology
(
    wsp_mean_na,
    wsp_stdm_na,
    wsp_hfit,
    wsp_x_data,
    wsp_residual,
    grid_cor_na,
) = regional_clima(
    wsp_mean_c,
    wsp_var_c,
    wsp_n_c,
    lon_w,
    lat_w,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_wsp_c,
    location,
    lsf="weighted",
    parameters=5,
)

# Perform weigthed least squares fit for basin-scale climatology projected onto regional climatology
hfit_reg_na, x_data_reg, x_data_sigma_reg = weighted_least_square_fit(
    swh_mean_na,
    swh_stdm_na,
    trend="sinusoidal",
    parameters=2,
    period=12,
    phase=np.arctan2(x_data_basin[1][1], x_data_basin[1][2]),
)

# Set variables for plotting
projection = ccrs.PlateCarree(central_longitude=180.0)
time = np.arange(1, 13, 1)
time_ticks = [
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
    "Jul",
    "Aug",
    "Sep",
    "Oct",
    "Nov",
    "Dec",
]  # x-axis tick labels for NH
xlim = [0, 13]
ylim = [[-1, 4], [5, 10]]
levels = np.arange(2, 8.1, 0.1)
resolution = "10m"

# initialize subplot axes:
fig = plt.figure(figsize=(22, 28))

############## Subplot 1  ##############
# California coast
ax1 = fig.add_subplot(321, projection=ccrs.PlateCarree(central_longitude=180.0))
cart.set_subplots(
    ax1, projection, resolution, lon_min=44, lon_max=66, lat_min=30, lat_max=44
)
cs1 = ax1.contourf(
    lon_w,
    lat_w,
    wsp_moments_mean[2, :, :],
    levels,
    cmap=cmo.thermal,
    extend="both",
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
cart.grid_lines_rc(
    ax=ax1,
    xticks=[-133, -129, -125, -121, -117],
    yticks=[28, 32, 36, 40, 44],
    fontsize=20,
    linewidth=1,
    color="gray",
    alpha=0.3,
    linestyle="--",
    grid=True,
)
cax1 = plt.axes([0.1, 0.06, 0.4, 0.03])
cart.set_cbar(
    cs1,
    cax1,
    fig,
    orientation="horizontal",
    extend="both",
    cbar_label="$m\,s^{-1}$",
    nbins=7,
    fontsize=20,
    cbar_ticks=[],
    task="regular",
)
cart.subplot_label(
    ax1,
    xdist_label=0.1,
    ydist_label=0.88,
    subplot_label="A",
    form="box",
    fs_shade=28,
    fs_main=20,
    color="black",
)

# Set regional climatology grid box
ax1.add_patch(
    mpatches.Rectangle(
        xy=[-129, 36],
        width=4,
        height=4,
        edgecolor="black",
        facecolor=None,
        linestyle="-",
        linewidth=2,
        fill=False,
        alpha=1,
        transform=ccrs.PlateCarree(),
    )
)

############## Subplot 2 ##############
# California Coast
ax2 = fig.add_subplot(322)
regional_clima_plot(
    ax2,
    swh_mean_nc,
    swh_stdm_nc,
    hfit_reg_nc,
    wsp_mean_nc,
    wsp_stdm_nc,
    None,
    None,
    None,
    None,
    None,
    time,
    time_ticks,
    xlim,
    ylim,
    subplot_label="B",
    fontsize=20,
    linewidth=1.5,
    task="C_I_res",
    grid=False,
)

############## Subplot 3  ##############
# South Caribbean
ax3 = fig.add_subplot(323, projection=ccrs.PlateCarree(central_longitude=180.0))
cart.set_subplots(
    ax3, projection, resolution, lon_min=97, lon_max=119, lat_min=7, lat_max=24
)
cs3 = ax3.contourf(
    lon_w,
    lat_w,
    wsp_moments_mean[2, :, :],
    levels,
    cmap=cmo.thermal,
    extend="both",
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
cart.grid_lines_rc(
    ax=ax3,
    xticks=[-84, -80, -76, -72, -68, -64, -60],
    yticks=[9, 13, 17, 21, 25],
    fontsize=20,
    linewidth=1,
    color="gray",
    alpha=0.3,
    linestyle="--",
    grid=True,
)
cart.subplot_label(
    ax3,
    xdist_label=0.1,
    ydist_label=0.88,
    subplot_label="C",
    form="box",
    fs_shade=28,
    fs_main=20,
    color="black",
)

# Set regional climatology grid box
ax3.add_patch(
    mpatches.Rectangle(
        xy=[-76, 13],
        width=4,
        height=4,
        edgecolor="black",
        facecolor=None,
        linestyle="-",
        linewidth=2,
        fill=False,
        alpha=1,
        transform=ccrs.PlateCarree(),
    )
)

############## Subplot 4 ##############
# South Caribbean Sea
ax4 = fig.add_subplot(324)
regional_clima_plot(
    ax4,
    swh_mean_sc,
    swh_stdm_sc,
    hfit_reg_sc,
    wsp_mean_sc,
    wsp_stdm_sc,
    None,
    None,
    None,
    None,
    None,
    time,
    time_ticks,
    xlim,
    ylim,
    subplot_label="D",
    fontsize=20,
    linewidth=1.5,
    task="C_I_res",
    grid=False,
)

############## Subplot 5  ##############
# North Africa
ax5 = fig.add_subplot(325, projection=ccrs.PlateCarree(central_longitude=180.0))
cart.set_subplots(
    ax5, projection, resolution, lon_min=154, lon_max=174, lat_min=19, lat_max=34
)
cs5 = ax5.contourf(
    lon_w,
    lat_w,
    wsp_moments_mean[2, :, :],
    levels,
    cmap=cmo.thermal,
    extend="both",
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
cart.grid_lines_rc(
    ax5,
    xticks=[-27, -23, -19, -15, -11, -7],
    yticks=[22, 26, 30, 34, 38, 42],
    fontsize=20,
    linewidth=1,
    color="gray",
    alpha=0.3,
    linestyle="--",
    grid=True,
)
cart.subplot_label(
    ax5,
    xdist_label=0.1,
    ydist_label=0.88,
    subplot_label="E",
    form="box",
    fs_shade=28,
    fs_main=20,
    color="black",
)

# Set regional climatology grid box
ax5.add_patch(
    mpatches.Rectangle(
        xy=[-15, 30],
        width=4,
        height=4,
        edgecolor="black",
        facecolor=None,
        linestyle="-",
        linewidth=2,
        fill=False,
        alpha=1,
        transform=ccrs.PlateCarree(),
    )
)

############## Subplot 6 ##############
# North Africa (Morocco)
ax6 = fig.add_subplot(326)
regional_clima_plot(
    ax6,
    swh_mean_na,
    swh_stdm_na,
    hfit_reg_na,
    wsp_mean_na,
    wsp_stdm_na,
    None,
    None,
    None,
    None,
    None,
    time,
    time_ticks,
    xlim,
    ylim,
    subplot_label="F",
    fontsize=20,
    linewidth=1.5,
    task="C_I_res",
    grid=False,
)

# adjust spacing for the entire figure (not the subplot)
fig.subplots_adjust(wspace=0.3, hspace=0.15)

# save figure
fig.savefig(fname="../figs/fig06", bbox_inches="tight", dpi=300)
