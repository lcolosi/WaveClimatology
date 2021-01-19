"""
Figure S1 Caption: 
(A) Ifremer SWH annual cycle phase map and (B) June through August seasonal 
probability of swell in Polynesian island region illustrating island shadowing.
"""

# Path to access python functions
import sys

sys.path.append("../tools/")

# Path to access intermediate data
data_path_lsf = "../data/lsf_parameters/"
data_path_ps = "../data/prob_swell/"

# libraries
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
import cmocean.cm as cmo
import matplotlib.patches as mpatches

# my functions
import cartopy_figs as cart

# Set filenames
filename_lsf = data_path_lsf + "IFREMER_swh_lsf_parameters.nc"
filename_ps = data_path_ps + "WW3_probability_swell.nc"

# Call wlsf phase parameter data:
nc_i = Dataset(filename_lsf, "r")
lon_i = nc_i.variables["lon"][:]
lat_i = nc_i.variables["lat"][:]
a_phase = nc_i.variables["a_phase"][:]

# Call probability of swell data:
nc_ps = Dataset(filename_ps, "r")
lon_ps = nc_ps.variables["lon"][:]
lat_ps = nc_ps.variables["lat"][:]
prob_swell = nc_ps.variables["seasonal_prob_swell"][:]

#Adjust prob_swell and longitude for plotting.
prob_swell, lon_ps = add_cyclic_point(prob_swell, coord=lon_ps)

# Initialize variables for plotting
projection = ccrs.PlateCarree(central_longitude=180.0)
xticks = [160, 180, -160, -140]
yticks = [-20, -10, 0, 10]
resolution = "10m"

# Create figure and axes
fig, axes = plt.subplots(2, 1, figsize=(16, 12), subplot_kw={"projection": projection})
ax1, ax2 = axes.flatten()

############## Subplot 1  #################
cart.set_subplots(
    ax1, projection, resolution, lon_min=-30, lon_max=50, lat_min=-25, lat_max=11
)
cs1 = ax1.pcolormesh(
    lon_i,
    lat_i,
    a_phase,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
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
    fontsize=15,
    color="black",
)
cax1 = plt.axes([0.82, 0.56, 0.015, 0.28])
cart.set_cbar(
    cs1,
    cax1,
    fig,
    orientation="vertical",
    extend="both",
    cbar_label="$months$",
    nbins=7,
    fontsize=15,
    cbar_ticks=[
        np.arange(-np.pi, np.pi + 0.5, (np.pi + np.pi) / 6).tolist(),
        ["", "Jun", "Aug", "Oct", "Dec", "Feb", "Apr", "Jun"],
    ],
    task="custom ticks",
)
cart.subplot_label(
    ax1,
    xdist_label=0.05,
    ydist_label=0.9,
    subplot_label="A",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="black",
)

# Place minor islands on map. This requires a working NaturalEarthFeature installation from cartopy. If you have that,
# you may uncomment this line:
# ax1.add_feature(cfeature.NaturalEarthFeature('physical', 'minor_islands', resolution))

# Set grid boxes for wave shadowing
ax1.add_patch(
    mpatches.Rectangle(
        xy=[176, -18],
        width=6,
        height=5,
        edgecolor="black",
        facecolor=None,
        linestyle="-",
        linewidth=2,
        fill=False,
        alpha=1,
        transform=ccrs.PlateCarree(),
    )
)
ax1.add_patch(
    mpatches.Rectangle(
        xy=[212, -17],
        width=6,
        height=5,
        edgecolor="black",
        facecolor=None,
        linestyle="-",
        linewidth=2,
        fill=False,
        alpha=1,
        transform=ccrs.PlateCarree(),
    )
)


############## Subplot 2  #################
cart.set_subplots(
    ax2, projection, resolution, lon_min=-30, lon_max=50, lat_min=-25, lat_max=11
)
cs2 = ax2.pcolormesh(
    lon_ps,
    lat_ps,
    prob_swell[2, :, :] * 100,
    vmin=75,
    vmax=100,
    cmap=cmo.thermal,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
cart.set_grid_ticks(
    ax2,
    projection=ccrs.PlateCarree(),
    xticks=xticks,
    yticks=yticks,
    xlabels=True,
    ylabels=True,
    grid=False,
    fontsize=15,
    color="black",
)
cax2 = plt.axes([0.82, 0.14, 0.015, 0.28])
cart.set_cbar(
    cs2,
    cax2,
    fig,
    orientation="vertical",
    extend="min",
    cbar_label="$\%$",
    nbins=6,
    fontsize=15,
    cbar_ticks=[],
    task="regular",
)
cart.subplot_label(
    ax2,
    xdist_label=0.05,
    ydist_label=0.9,
    subplot_label="B",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="black",
)

# Place minor islands on map. This requires a working NaturalEarthFeature installation from cartopy. If you have that,
# you may uncomment this line:
# ax2.add_feature(cfeature.NaturalEarthFeature('physical', 'minor_islands', resolution))

# adjust spacing for the entire figure
fig.subplots_adjust(wspace=0.2, hspace=0.2)

# save figure
fig.savefig(fname="../figs/figS01", bbox_inches="tight", dpi=300)
