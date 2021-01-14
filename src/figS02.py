"""
Figure S2 Caption: 
Fractional uncertainty of amplitude of annual cycle for (A) IFREMER SWH and (B)
CCMP2 WSP; amplitude of semi-annual cycle for (C) IFREMER SWH and (D) CCMP2 WSP.
"""

# Path to access python functions
import sys

sys.path.append("../tools/")

# Path to access intermediate data
data_path = "../data/"

# libraries
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cmocean.cm as cmo
import matplotlib.patches as mpatches

# my functions
import cartopy_figs as cart

# Set filenames
filename_ifremer = data_path + "IFREMER_swh_lsf_parameters.nc"
filename_ccmp2 = data_path + "CCMP2_wsp_lsf_parameters.nc"

# Call data
######### Ifremer #########
nc_i = Dataset(filename_ifremer, "r")
lon = nc_i.variables["lon"][:]
lat = nc_i.variables["lat"][:]
a_amp_i = nc_i.variables["a_amp"][:]
s_amp_i = nc_i.variables["s_amp"][:]
a_amp_unc_i = nc_i.variables["a_amp_unc"][:]
s_amp_unc_i = nc_i.variables["s_amp_unc"][:]

######### CCMP2 #########
nc_c = Dataset(filename_ccmp2, "r")
a_amp_c = nc_c.variables["a_amp"][:]
s_amp_c = nc_c.variables["s_amp"][:]
a_amp_unc_c = nc_c.variables["a_amp_unc"][:]
s_amp_unc_c = nc_c.variables["s_amp_unc"][:]

# Compute the relative uncertainty
a_ratio_i = a_amp_unc_i / a_amp_i
s_ratio_i = s_amp_unc_i / s_amp_i
a_ratio_c = a_amp_unc_c / a_amp_c
s_ratio_c = s_amp_unc_c / s_amp_c

# Initialize variables for plotting
projection = ccrs.PlateCarree(central_longitude=180.0)
xticks = [0, 45, 90, 135, -360, -180, -135, -90, -45]
yticks = [-60, -40, -20, 0, 20, 40, 60]
levels = np.arange(0.0, 2 + 0.1, 0.1)
resolution = "110m"

# Create figure and axes
fig, axes = plt.subplots(2, 2, figsize=(14, 6), subplot_kw={"projection": projection})
ax1, ax2, ax3, ax4 = axes.flatten()

############## Subplot 1  #################
cart.set_subplots(
    ax1, projection, resolution, lon_min=-180, lon_max=179, lat_min=-66, lat_max=66
)
cs1 = ax1.contourf(
    lon,
    lat,
    a_ratio_i,
    levels,
    cmap=cmo.diff,
    extend="both",
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
cart.set_grid_ticks(
    ax1,
    projection=ccrs.PlateCarree(),
    xticks=xticks,
    yticks=yticks,
    xlabels=False,
    ylabels=True,
    grid=False,
    fontsize=12,
    color="black",
)
cax1 = plt.axes([0.48, 0.35, 0.012, 0.32])
cart.set_cbar(
    cs1,
    cax1,
    fig,
    orientation="vertical",
    extend="both",
    cbar_label="",
    nbins=5,
    fontsize=12,
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
    fs_main=18,
    color="black",
)

############## Subplot 2  #################
cart.set_subplots(
    ax2, projection, resolution, lon_min=-180, lon_max=179, lat_min=-66, lat_max=66
)
cs2 = ax2.contourf(
    lon,
    lat,
    a_ratio_c,
    levels,
    cmap=cmo.balance,
    extend="both",
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
cart.set_grid_ticks(
    ax2,
    projection=ccrs.PlateCarree(),
    xticks=xticks,
    yticks=yticks,
    xlabels=False,
    ylabels=False,
    grid=False,
    fontsize=12,
    color="black",
)
cax2 = plt.axes([0.91, 0.35, 0.012, 0.32])
cart.set_cbar(
    cs2,
    cax2,
    fig,
    orientation="vertical",
    extend="both",
    cbar_label="",
    nbins=5,
    fontsize=12,
    cbar_ticks=[],
    task="regular",
)
cart.subplot_label(
    ax2,
    xdist_label=0.1,
    ydist_label=0.88,
    subplot_label="B",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="black",
)

############## Subplot 3  #################
cart.set_subplots(
    ax3, projection, resolution, lon_min=-180, lon_max=179, lat_min=-66, lat_max=66
)
cs3 = ax3.contourf(
    lon,
    lat,
    s_ratio_i,
    levels,
    cmap=cmo.diff,
    extend="both",
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
cart.set_grid_ticks(
    ax3,
    projection=ccrs.PlateCarree(),
    xticks=xticks,
    yticks=yticks,
    xlabels=True,
    ylabels=True,
    grid=False,
    fontsize=12,
    color="black",
)
cart.subplot_label(
    ax3,
    xdist_label=0.1,
    ydist_label=0.88,
    subplot_label="C",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="black",
)

############## Subplot 4  #################
cart.set_subplots(
    ax4, projection, resolution, lon_min=-180, lon_max=179, lat_min=-66, lat_max=66
)
cs4 = ax4.contourf(
    lon,
    lat,
    s_ratio_c,
    levels,
    cmap=cmo.balance,
    extend="both",
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
cart.set_grid_ticks(
    ax4,
    projection=ccrs.PlateCarree(),
    xticks=xticks,
    yticks=yticks,
    xlabels=True,
    ylabels=False,
    grid=False,
    fontsize=12,
    color="black",
)
cart.subplot_label(
    ax4,
    xdist_label=0.1,
    ydist_label=0.88,
    subplot_label="D",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="black",
)

# adjust spacing for the entire figure
fig.subplots_adjust(wspace=0.25, hspace=0.02)

# save figure
fig.savefig(fname="../figs/figS02", bbox_inches="tight", dpi=300)
