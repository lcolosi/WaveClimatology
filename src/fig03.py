"""
Figure 3 Caption: 
Fraction of variance explained by weighted annual and semi-annual least squares
fit for IFREMER SWH (A) and CCMP2 WSP (B) from January 1st, 1993 to December 
31st, 2015.
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
from matplotlib import cm

# my functions
import cartopy_figs as cart

# Call data:
filename_ifremer = data_path + "IFREMER_swh_lsf_parameters.nc"
filename_ccmp2 = data_path + "CCMP2_wsp_lsf_parameters.nc"

######### Ifremer #########
nc_i = Dataset(filename_ifremer, "r")
lon = nc_i.variables["lon"][:]
lat = nc_i.variables["lat"][:]
fve_i = nc_i.variables["fve"][:]

######### CCMP2 #########
nc_c = Dataset(filename_ccmp2, "r")
fve_c = nc_c.variables["fve"][:]

# Initialize variables for plotting
projection = ccrs.PlateCarree(central_longitude=180.0)
levels = np.arange(75, 101, 1)
xticks = [0, 45, 90, 135, -360, -180, -135, -90, -45]
yticks = [-60, -40, -20, 0, 20, 40, 60]
resolution = "110m"

# Create figure and axes
fig, axes = plt.subplots(1, 2, figsize=(24, 20), subplot_kw={"projection": projection})
ax1, ax2 = axes.flatten()

############## Subplot 1  #################
cart.set_subplots(
    ax1, projection, resolution, lon_min=-180, lon_max=179, lat_min=-66, lat_max=66
)
cs1 = ax1.contourf(
    lon,
    lat,
    fve_i * 100,
    levels=levels,
    cmap=cm.YlOrRd_r,
    extend="min",
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
cart.set_grid_ticks(
    ax1,
    projection=ccrs.PlateCarree(),
    xticks=xticks,
    yticks=yticks,
    xlabels=True,
    ylabels=True,
    grid=False,
    fontsize=15,
    color="black",
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
    fve_c * 100,
    levels=levels,
    cmap=cm.YlOrRd_r,
    extend="min",
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
cart.set_grid_ticks(
    ax2,
    projection=ccrs.PlateCarree(),
    xticks=xticks,
    yticks=yticks,
    xlabels=True,
    ylabels=False,
    grid=False,
    fontsize=15,
    color="black",
)
cax2 = plt.axes([0.36, 0.36, 0.3, 0.02])
cart.set_cbar(
    cs2,
    cax2,
    fig,
    orientation="horizontal",
    extend="both",
    cbar_label="$\%$",
    nbins=6,
    fontsize=20,
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

# adjust spacing for the entire figure (not the subplot)
fig.subplots_adjust(wspace=0.05, hspace=0.2)

# save figure
fig.savefig(fname="../figs/fig03", bbox_inches="tight", dpi=300)
