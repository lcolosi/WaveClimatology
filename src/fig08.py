"""
Figure 8 Caption: 
Fraction of variance explained by weighted annual and semi-annual least squares
Seasonal progression of probability of swell using wave age criterion (3) and 
WW3 peak frequency and WSP from January 1st, 1993 to December 31st, 2015 where 
(A) DJF, (B) MAM, (C) JJA, and (D) SON.
"""

# Path to access python functions
import sys

sys.path.append("../tools/")

# Path to access intermediate data
data_path = "../data/prob_swell/"

# libraries
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
import cmocean.cm as cmo

# my functions
import cartopy_figs as cart

# Set filename
filename = data_path + "WW3_probability_swell.nc"

# Call data:
nc = Dataset(filename, "r")
lon = nc.variables["lon"][:]
lat = nc.variables["lat"][:]
prob_swell = nc.variables["seasonal_prob_swell"][:]

# Adjust prob_swell and longitude for plotting.
prob_swell, lon = add_cyclic_point(prob_swell, coord=lon)

# Initialize variables for plotting
projection = ccrs.PlateCarree(central_longitude=180.0)
levels = np.arange(75, 101, 1)
xticks = [0, 45, 90, 135, -360, -180, -135, -90, -45]
yticks = [-60, -40, -20, 0, 20, 40, 60]
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
    prob_swell[0] * 100,
    levels=levels,
    cmap=cmo.thermal,
    extend="min",
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
    fontsize=14,
    color="black",
)
cax1 = plt.axes([0.40, 0.07, 0.2, 0.02])
cart.set_cbar(
    cs1,
    cax1,
    fig,
    orientation="horizontal",
    extend="min",
    cbar_label="$\%$",
    nbins=6,
    fontsize=14,
    cbar_ticks=[],
    task="regular",
)
cart.subplot_label(
    ax1,
    xdist_label=0.1,
    ydist_label=0.86,
    subplot_label="A",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="black",
)
cart.subplot_label(
    ax1,
    xdist_label=0.24,
    ydist_label=0.86,
    subplot_label="DJF",
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
    prob_swell[1] * 100,
    levels=levels,
    cmap=cmo.thermal,
    extend="min",
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
    fontsize=14,
    color="black",
)
cart.subplot_label(
    ax2,
    xdist_label=0.1,
    ydist_label=0.86,
    subplot_label="B",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="black",
)
cart.subplot_label(
    ax2,
    xdist_label=0.24,
    ydist_label=0.86,
    subplot_label="MAM",
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
    prob_swell[2] * 100,
    levels=levels,
    cmap=cmo.thermal,
    extend="min",
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
    fontsize=14,
    color="black",
)
cart.subplot_label(
    ax3,
    xdist_label=0.1,
    ydist_label=0.86,
    subplot_label="C",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="black",
)
cart.subplot_label(
    ax3,
    xdist_label=0.24,
    ydist_label=0.86,
    subplot_label="JJA",
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
    prob_swell[3] * 100,
    levels=levels,
    cmap=cmo.thermal,
    extend="min",
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
    fontsize=14,
    color="black",
)
cart.subplot_label(
    ax4,
    xdist_label=0.1,
    ydist_label=0.86,
    subplot_label="D",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="black",
)
cart.subplot_label(
    ax4,
    xdist_label=0.24,
    ydist_label=0.86,
    subplot_label="SON",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="black",
)

# adjust spacing for the entire figure (not the subplot)
fig.subplots_adjust(wspace=0.08, hspace=0.001)

# save figure
fig.savefig(fname="../figs/fig08", bbox_inches="tight", dpi=300)
