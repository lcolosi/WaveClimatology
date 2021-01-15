"""
Figure S9 Caption:
Regional climatologies with Ifremer SWH (Solid blue), CCMP2 (Solid red), WW3 SWH
(dashed blue), and WW3 WSP (dashed red). Same regions as in Figure~7.
"""

# Path to access python functions
import sys

sys.path.append("../tools/")

# Path to access intermediate data
data_path_i = "../data/ifremer_swh/"
data_path_c = "../data/ccmp2_wsp/"
data_path_ws = "../data/ww3_swh/"
data_path_ww = "../data/ww3_wsp/"
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

# call IFREMER SWH, CCMP2 WSP, WW3 SWH, and WW3 WSP processed data:
swh_i, time_i, lat_i, lon_i = import_data("IFREMER_swh", data_path_i)
wsp_c, time_c, lat_c, lon_c = import_data("CCMP2_wsp", data_path_c)
swh_w, time_ws, lat_ws, lon_ws = import_data("WW3_swh", data_path_ws)
wsp_w, time_ww, lat_ww, lon_ww = import_data("WW3_wsp", data_path_ww)

# Call decorrelation time scales
###### IFREMER SWH ######
nc_swh_i = Dataset(data_path_decor + "IFREMER_swh_decor_time_scale.nc", "r")
decor_swh_i = nc_swh_i.variables["decor_scale"][:]
time_decor_swh_i = num2date(
    nc_swh_i.variables["time"][:], nc_swh_i.variables["time"].units
)
###### CCMP2 WSP ######
nc_wsp_c = Dataset(data_path_decor + "CCMP2_wsp_decor_time_scale.nc", "r")
decor_wsp_c = nc_wsp_c.variables["decor_scale"][:]
time_decor_wsp_c = num2date(
    nc_wsp_c.variables["time"][:], nc_wsp_c.variables["time"].units
)
#### WW3 SWH ####
nc_swh_w = Dataset(data_path_decor + "WW3_swh_decor_time_scale.nc", "r")
decor_swh_w = nc_swh_w.variables["decor_scale"][:]
time_decor_swh_w = num2date(
    nc_swh_w.variables["time"][:], nc_swh_w.variables["time"].units
)
#### WW3 WSP ####
nc_wsp_w = Dataset(data_path_decor + "WW3_wsp_decor_time_scale.nc", "r")
decor_wsp_w = nc_wsp_w.variables["decor_scale"][:]
time_decor_wsp_w = num2date(
    nc_wsp_w.variables["time"][:], nc_wsp_w.variables["time"].units
)

# Calculate monthly climatologies for SWH and WSP data and decorrelation scales:
###### IFREMER SWH ######
swh_clima_dict = clima_mean(date_time=np.ma.array(time_i), data=swh_i)
swh_mean_i = np.ma.array(swh_clima_dict["mean"])
swh_var_i = np.ma.array(swh_clima_dict["var"])
swh_n_i = np.ma.array(swh_clima_dict["N"])
decor_clima_dict_s = clima_mean(date_time=time_decor_swh_i, data=decor_swh_i)
decor_mean_swh_i = np.ma.array(decor_clima_dict_s["mean"])

###### CCMP2 WSP ######
wsp_clima_dict = clima_mean(date_time=np.ma.array(time_c), data=wsp_c)
wsp_mean_c = np.ma.array(wsp_clima_dict["mean"])
wsp_var_c = np.ma.array(wsp_clima_dict["var"])
wsp_n_c = np.ma.array(wsp_clima_dict["N"])
decor_clima_dict_w = clima_mean(date_time=time_decor_wsp_c, data=decor_wsp_c)
decor_mean_wsp_c = np.ma.array(decor_clima_dict_w["mean"])

#### WW3 SWH ####
swh_clima_dict = clima_mean(date_time=np.ma.array(time_ws), data=swh_w)
swh_mean_w = np.ma.array(swh_clima_dict["mean"])
swh_var_w = np.ma.array(swh_clima_dict["var"])
swh_n_w = np.ma.array(swh_clima_dict["N"])
decor_clima_dict_s = clima_mean(date_time=time_decor_swh_w, data=decor_swh_w)
decor_mean_swh_w = np.ma.array(decor_clima_dict_s["mean"])

#### WW3 WSP ####
wsp_clima_dict = clima_mean(date_time=np.ma.array(time_ww), data=wsp_w)
wsp_mean_w = np.ma.array(wsp_clima_dict["mean"])
wsp_var_w = np.ma.array(wsp_clima_dict["var"])
wsp_n_w = np.ma.array(wsp_clima_dict["N"])
decor_clima_dict_w = clima_mean(date_time=time_decor_wsp_w, data=decor_wsp_w)
decor_mean_wsp_w = np.ma.array(decor_clima_dict_w["mean"])

# Compute regional climatologies:

####### Northern California #######
# intialize variable
ngrid = 4
lon_grid = 232
lat_grid = 103
location = ["NH", "west"]

# compute ifremer swh regional climatology
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
    lon_i,
    lat_i,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_swh_i,
    location,
    lsf="weighted",
    parameters=5,
)

# compute ccmp2 wsp regional climatology
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
    lon_c,
    lat_c,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_wsp_c,
    location,
    lsf="weighted",
    parameters=5,
)

# compute ww3 swh regional climatology
(
    swh_m_mean_nc,
    swh_m_stdm_nc,
    swh_m_hfit,
    swh_x_data,
    swh_residual,
    grid_cor_nc,
) = regional_clima(
    swh_mean_w,
    swh_var_w,
    swh_n_w,
    lon_ws,
    lat_ws,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_swh_w,
    location,
    lsf="weighted",
    parameters=5,
)

# compute ww3 wsp regional climatology
(
    wsp_m_mean_nc,
    wsp_m_stdm_nc,
    wsp_m_hfit,
    wsp_x_data,
    wsp_residual,
    grid_cor_nc,
) = regional_clima(
    wsp_mean_w,
    wsp_var_w,
    wsp_n_w,
    lon_ws,
    lat_ws,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_wsp_w,
    location,
    lsf="weighted",
    parameters=5,
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
    lon_i,
    lat_i,
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
    lon_c,
    lat_c,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_wsp_c,
    location,
    lsf="weighted",
    parameters=5,
)

# compute ww3 swh regional climatology
(
    swh_m_mean_sc,
    swh_m_stdm_sc,
    swh_m_hfit,
    swh_x_data,
    swh_residual,
    grid_cor_sc,
) = regional_clima(
    swh_mean_w,
    swh_var_w,
    swh_n_w,
    lon_ws,
    lat_ws,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_swh_w,
    location,
    lsf="weighted",
    parameters=5,
)

# compute ww3 wsp regional climatology
(
    wsp_m_mean_sc,
    wsp_m_stdm_sc,
    wsp_m_hfit,
    wsp_x_data,
    wsp_residual,
    grid_cor_sc,
) = regional_clima(
    wsp_mean_w,
    wsp_var_w,
    wsp_n_w,
    lon_ws,
    lat_ws,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_wsp_w,
    location,
    lsf="weighted",
    parameters=5,
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
    lon_i,
    lat_i,
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
    lon_c,
    lat_c,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_wsp_c,
    location,
    lsf="weighted",
    parameters=5,
)

# compute ww3 swh regional climatology
(
    swh_m_mean_na,
    swh_m_stdm_na,
    swh_m_hfit,
    swh_x_data,
    swh_residual,
    grid_cor_na,
) = regional_clima(
    swh_mean_w,
    swh_var_w,
    swh_n_w,
    lon_ws,
    lat_ws,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_swh_w,
    location,
    lsf="weighted",
    parameters=5,
)

# compute ww3 wsp regional climatology
(
    wsp_m_mean_na,
    wsp_m_stdm_na,
    wsp_m_hfit,
    wsp_x_data,
    wsp_residual,
    grid_cor_na,
) = regional_clima(
    wsp_mean_w,
    wsp_var_w,
    wsp_n_w,
    lon_ws,
    lat_ws,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_wsp_w,
    location,
    lsf="weighted",
    parameters=5,
)

####### Western Australia #######
# intialize variable
ngrid = 4
lon_grid = 108
lat_grid = 35
location = ["SH", "east"]

# compute SWH regional climatology
(
    swh_mean_wa,
    swh_stdm_wa,
    swh_hfit,
    swh_x_data,
    swh_residual,
    grid_cor_wa,
) = regional_clima(
    swh_mean_i,
    swh_var_i,
    swh_n_i,
    lon_i,
    lat_i,
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
    wsp_mean_wa,
    wsp_stdm_wa,
    wsp_hfit,
    wsp_x_data,
    wsp_residual,
    grid_cor_wa,
) = regional_clima(
    wsp_mean_c,
    wsp_var_c,
    wsp_n_c,
    lon_c,
    lat_c,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_wsp_c,
    location,
    lsf="weighted",
    parameters=5,
)

# compute ww3 swh regional climatology
(
    swh_m_mean_wa,
    swh_m_stdm_wa,
    swh_m_hfit,
    swh_x_data,
    swh_residual,
    grid_cor_wa,
) = regional_clima(
    swh_mean_w,
    swh_var_w,
    swh_n_w,
    lon_ws,
    lat_ws,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_swh_w,
    location,
    lsf="weighted",
    parameters=5,
)

# compute ww3 wsp regional climatology
(
    wsp_m_mean_wa,
    wsp_m_stdm_wa,
    wsp_m_hfit,
    wsp_x_data,
    wsp_residual,
    grid_cor_wa,
) = regional_clima(
    wsp_mean_w,
    wsp_var_w,
    wsp_n_w,
    lon_ws,
    lat_ws,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_wsp_w,
    location,
    lsf="weighted",
    parameters=5,
)

####### Peru-Chile Coast #######
# intialize variable
ngrid = 4
lon_grid = 283
lat_grid = 30
location = ["SH", "west"]

# compute SWH regional climatology
(
    swh_mean_pc,
    swh_stdm_pc,
    swh_hfit,
    swh_x_data,
    swh_residual,
    grid_cor_pc,
) = regional_clima(
    swh_mean_i,
    swh_var_i,
    swh_n_i,
    lon_i,
    lat_i,
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
    wsp_mean_pc,
    wsp_stdm_pc,
    wsp_hfit,
    wsp_x_data,
    wsp_residual,
    grid_cor_pc,
) = regional_clima(
    wsp_mean_c,
    wsp_var_c,
    wsp_n_c,
    lon_c,
    lat_c,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_wsp_c,
    location,
    lsf="weighted",
    parameters=5,
)

# compute ww3 swh regional climatology
(
    swh_m_mean_pc,
    swh_m_stdm_pc,
    swh_m_hfit,
    swh_x_data,
    swh_residual,
    grid_cor_pc,
) = regional_clima(
    swh_mean_w,
    swh_var_w,
    swh_n_w,
    lon_ws,
    lat_ws,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_swh_w,
    location,
    lsf="weighted",
    parameters=5,
)

# compute ww3 wsp regional climatology
(
    wsp_m_mean_pc,
    wsp_m_stdm_pc,
    wsp_m_hfit,
    wsp_x_data,
    wsp_residual,
    grid_cor_pc,
) = regional_clima(
    wsp_mean_w,
    wsp_var_w,
    wsp_n_w,
    lon_ws,
    lat_ws,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_wsp_w,
    location,
    lsf="weighted",
    parameters=5,
)

####### South Africa (Namibian) #######
# intialize variable
ngrid = 4
lon_grid = 11
lat_grid = 36
location = ["SH", "east"]

# compute SWH regional climatology
(
    swh_mean_sa,
    swh_stdm_sa,
    swh_hfit,
    swh_x_data,
    swh_residual,
    grid_cor_sa,
) = regional_clima(
    swh_mean_i,
    swh_var_i,
    swh_n_i,
    lon_i,
    lat_i,
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
    wsp_mean_sa,
    wsp_stdm_sa,
    wsp_hfit,
    wsp_x_data,
    wsp_residual,
    grid_cor_sa,
) = regional_clima(
    wsp_mean_c,
    wsp_var_c,
    wsp_n_c,
    lon_c,
    lat_c,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_wsp_c,
    location,
    lsf="weighted",
    parameters=5,
)

# compute ww3 swh regional climatology
(
    swh_m_mean_sa,
    swh_m_stdm_sa,
    swh_m_hfit,
    swh_x_data,
    swh_residual,
    grid_cor_sa,
) = regional_clima(
    swh_mean_w,
    swh_var_w,
    swh_n_w,
    lon_ws,
    lat_ws,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_swh_w,
    location,
    lsf="weighted",
    parameters=5,
)

# compute ww3 wsp regional climatology
(
    wsp_m_mean_sa,
    wsp_m_stdm_sa,
    wsp_m_hfit,
    wsp_x_data,
    wsp_residual,
    grid_cor_sa,
) = regional_clima(
    wsp_mean_w,
    wsp_var_w,
    wsp_n_w,
    lon_ws,
    lat_ws,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_wsp_w,
    location,
    lsf="weighted",
    parameters=5,
)

####### Arabian Sea #######
# intialize variable
ngrid = 4
lon_grid = 51
lat_grid = 70
location = ["SH", "east"]

# compute SWH regional climatology
(
    swh_mean_as,
    swh_stdm_as,
    swh_hfit,
    swh_x_data,
    swh_residual,
    grid_cor_as,
) = regional_clima(
    swh_mean_i,
    swh_var_i,
    swh_n_i,
    lon_i,
    lat_i,
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
    wsp_mean_as,
    wsp_stdm_as,
    wsp_hfit,
    wsp_x_data,
    wsp_residual,
    grid_cor_as,
) = regional_clima(
    wsp_mean_c,
    wsp_var_c,
    wsp_n_c,
    lon_c,
    lat_c,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_wsp_c,
    location,
    lsf="weighted",
    parameters=5,
)

# compute ww3 swh regional climatology
(
    swh_m_mean_as,
    swh_m_stdm_as,
    swh_m_hfit,
    swh_x_data,
    swh_residual,
    grid_cor_as,
) = regional_clima(
    swh_mean_w,
    swh_var_w,
    swh_n_w,
    lon_ws,
    lat_ws,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_swh_w,
    location,
    lsf="weighted",
    parameters=5,
)

# compute ww3 wsp regional climatology
(
    wsp_m_mean_as,
    wsp_m_stdm_as,
    wsp_m_hfit,
    wsp_x_data,
    wsp_residual,
    grid_cor_as,
) = regional_clima(
    wsp_mean_w,
    wsp_var_w,
    wsp_n_w,
    lon_ws,
    lat_ws,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_wsp_w,
    location,
    lsf="weighted",
    parameters=5,
)

# Set variables for plotting
time = np.arange(1, 13, 1)
time_ticks_nh = [
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
]
time_ticks_sh = [
    "Jul",
    "Aug",
    "Sep",
    "Oct",
    "Nov",
    "Dec",
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
]
xlim = [0, 13]
ylim = [[1, 4], [4, 10]]

# Create figure and axes
fig, axes = plt.subplots(4, 2, figsize=(16, 20))
ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8 = axes.flatten()

########### Subplot 1 ###########
# California Coast
regional_clima_plot(
    ax1,
    swh_mean_nc,
    swh_stdm_nc,
    None,
    wsp_mean_nc,
    wsp_stdm_nc,
    None,
    swh_m_mean_nc,
    swh_m_stdm_nc,
    wsp_m_mean_nc,
    wsp_m_stdm_nc,
    time,
    time_ticks_nh,
    xlim,
    ylim,
    subplot_label="A",
    fontsize=20,
    linewidth=1.5,
    task="ww3",
    grid=False,
)

########### Subplot 2 ###########
# West Coast of Australia
regional_clima_plot(
    ax2,
    swh_mean_wa,
    swh_stdm_wa,
    None,
    wsp_mean_wa,
    wsp_stdm_wa,
    None,
    swh_m_mean_wa,
    swh_m_stdm_wa,
    wsp_m_mean_wa,
    wsp_m_stdm_wa,
    time,
    time_ticks_sh,
    xlim,
    ylim,
    subplot_label="B",
    fontsize=20,
    linewidth=1.5,
    task="ww3",
    grid=False,
)

########### Subplot 3 ###########
# South Caribbean Sea
regional_clima_plot(
    ax3,
    swh_mean_sc,
    swh_stdm_sc,
    None,
    wsp_mean_sc,
    wsp_stdm_sc,
    None,
    swh_m_mean_sc,
    swh_m_stdm_sc,
    wsp_m_mean_sc,
    wsp_m_stdm_sc,
    time,
    time_ticks_nh,
    xlim,
    ylim,
    subplot_label="C",
    fontsize=20,
    linewidth=1.5,
    task="ww3",
    grid=False,
)

########### Subplot 4 ###########
# Peru-Chile Coast
regional_clima_plot(
    ax4,
    swh_mean_pc,
    swh_stdm_pc,
    None,
    wsp_mean_pc,
    wsp_stdm_pc,
    None,
    swh_m_mean_pc,
    swh_m_stdm_pc,
    wsp_m_mean_pc,
    wsp_m_stdm_pc,
    time,
    time_ticks_sh,
    xlim,
    ylim,
    subplot_label="D",
    fontsize=20,
    linewidth=1.5,
    task="ww3",
    grid=False,
)

########### Subplot 5 ###########
# North Africa (Morocco)
regional_clima_plot(
    ax5,
    swh_mean_na,
    swh_stdm_na,
    None,
    wsp_mean_na,
    wsp_stdm_na,
    None,
    swh_m_mean_na,
    swh_m_stdm_na,
    wsp_m_mean_na,
    wsp_m_stdm_na,
    time,
    time_ticks_nh,
    xlim,
    ylim,
    subplot_label="E",
    fontsize=20,
    linewidth=1.5,
    task="ww3",
    grid=False,
)

########### Subplot 6 ###########
# South Africa Coast (Namibia)
regional_clima_plot(
    ax6,
    swh_mean_sa,
    swh_stdm_sa,
    None,
    wsp_mean_sa,
    wsp_stdm_sa,
    None,
    swh_m_mean_sa,
    swh_m_stdm_sa,
    wsp_m_mean_sa,
    wsp_m_stdm_sa,
    time,
    time_ticks_sh,
    xlim,
    ylim,
    subplot_label="F",
    fontsize=20,
    linewidth=1.5,
    task="ww3",
    grid=False,
)

########### Subplot 7 ###########
# suppress output for axis 7
ax7.axis("off")

########### Subplot 8 ###########
# Arabian Sea
regional_clima_plot(
    ax8,
    swh_mean_as,
    swh_stdm_as,
    None,
    wsp_mean_as,
    wsp_stdm_as,
    None,
    swh_m_mean_as,
    swh_m_stdm_as,
    wsp_m_mean_as,
    wsp_m_stdm_as,
    time,
    time_ticks_sh,
    xlim,
    ylim=[[0.5, 4], [3, 12]],
    subplot_label="G",
    fontsize=20,
    linewidth=1.5,
    task="ww3",
    grid=False,
)

# make sure the figure looks good:
fig.tight_layout()

# adjust spacing for the entire figure
fig.subplots_adjust(wspace=0.35, hspace=0.2)

# save figure:
fig.savefig(fname="../figs/figS09", bbox_inches="tight", dpi=300)
