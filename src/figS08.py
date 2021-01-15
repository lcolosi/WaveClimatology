"""
Figure S8 Caption: 
Phase of annual cycle for (A) WW3 SWH and (B) CFSR WSP; phase of semi-annual 
cycle for (C) WW3 SWH and (D) CFSR WSP; amplitude of annual cycle for (E) WW3 
SWH and (F) CFSR WSP; amplitude of semi-annual cycle for (G) WW3 SWH and (H) 
CFSR WSP. Grid points with a amplitude less than or equal to 2 standard deviations
are considered not statistically significant and masked white; the same pixels are
also masked for phase. See section 2.3 for details of computation.
"""

# Path to access python functions
import sys

sys.path.append("../tools/")

# Path to access intermediate data
data_path = "../data/lsf_parameters/"

# libraries
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cmocean.cm as cmo

# my functions
import cartopy_figs as cart

# Call data:
filename_swh = data_path + "WW3_swh_lsf_parameters.nc"
filename_wsp = data_path + "WW3_wsp_lsf_parameters.nc"

######### SWH #########
nc_s = Dataset(filename_swh, "r")
lon = nc_s.variables["lon"][:]
lat = nc_s.variables["lat"][:]
a_amp_s = nc_s.variables["a_amp"][:]
a_phase_s = nc_s.variables["a_phase"][:]
s_amp_s = nc_s.variables["s_amp"][:]
s_phase_s = nc_s.variables["s_phase"][:]
a_amp_unc_s = nc_s.variables["a_amp_unc"][:]
a_phase_unc_s = nc_s.variables["a_phase_unc"][:]
s_amp_unc_s = nc_s.variables["s_amp_unc"][:]
s_phase_unc_s = nc_s.variables["s_phase_unc"][:]

######### WSP #########
nc_w = Dataset(filename_wsp, "r")
a_amp_w = nc_w.variables["a_amp"][:]
a_phase_w = nc_w.variables["a_phase"][:]
s_amp_w = nc_w.variables["s_amp"][:]
s_phase_w = nc_w.variables["s_phase"][:]
a_amp_unc_w = nc_w.variables["a_amp_unc"][:]
a_phase_unc_w = nc_w.variables["a_phase_unc"][:]
s_amp_unc_w = nc_w.variables["s_amp_unc"][:]
s_phase_unc_w = nc_w.variables["s_phase_unc"][:]

# Set noise to signal ratio criteria
ns = 5 / 10

# Mask not statistical significance swh and wsp data
######### SWH #########
# Compute the relative uncertainty
a_ratio_s = a_amp_unc_s / a_amp_s
s_ratio_s = s_amp_unc_s / s_amp_s
# Mask not statistically significant grid points
a_mask_s = np.ma.getmask(np.ma.masked_greater_equal(a_ratio_s, (ns)))
s_mask_s = np.ma.getmask(np.ma.masked_greater_equal(s_ratio_s, (ns)))
# Apply statistical significance masks
a_phase_s_m = np.ma.masked_where(a_mask_s, a_phase_s)
a_amp_s_m = np.ma.masked_where(a_mask_s, a_amp_s)
s_phase_s_m = np.ma.masked_where(s_mask_s, s_phase_s)
s_amp_s_m = np.ma.masked_where(s_mask_s, s_amp_s)

######### WSP #########
# Compute the relative uncertainty
a_ratio_w = a_amp_unc_w / a_amp_w
s_ratio_w = s_amp_unc_w / s_amp_w
# Mask not statistically significant grid points
a_mask_w = np.ma.getmask(np.ma.masked_greater_equal(a_ratio_w, (ns)))
s_mask_w = np.ma.getmask(np.ma.masked_greater_equal(s_ratio_w, (ns)))
# Apply statistical significance masks
a_phase_w_m = np.ma.masked_where(a_mask_w, a_phase_w)
a_amp_w_m = np.ma.masked_where(a_mask_w, a_amp_w)
s_phase_w_m = np.ma.masked_where(s_mask_w, s_phase_w)
s_amp_w_m = np.ma.masked_where(s_mask_w, s_amp_w)

# Initialize variables for plotting
projection = ccrs.PlateCarree(central_longitude=180.0)
levels = np.arange(0, 1.6, 0.1)
resolution = "110m"

# Create figure and axes
fig, axes = plt.subplots(4, 2, figsize=(16, 12), subplot_kw={"projection": projection})
ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8 = axes.flatten()

############## Subplot 1  #################
cart.set_subplots(
    ax1, projection, resolution, lon_min=-180, lon_max=179, lat_min=-66, lat_max=66
)
cs1 = ax1.pcolormesh(
    lon,
    lat,
    a_phase_s_m,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=180.0),
)
cart.set_grid_ticks(
    ax1,
    projection=ccrs.PlateCarree(),
    xticks=[0, 45, 90, 135, -360, -180, -135, -90, -45 - 0],
    yticks=[-60, -40, -20, 0, 20, 40, 60],
    xlabels=False,
    ylabels=True,
    grid=False,
    fontsize=12,
    color="black",
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

############## Subplot 2  #################
cart.set_subplots(
    ax2, projection, resolution, lon_min=-180, lon_max=179, lat_min=-66, lat_max=66
)
cs2 = ax2.pcolormesh(
    lon,
    lat,
    a_phase_w_m,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=180.0),
)
cart.set_grid_ticks(
    ax2,
    projection=ccrs.PlateCarree(),
    xticks=[0, 45, 90, 135, -360, -180, -135, -90, -45],
    yticks=[-60, -40, -20, 0, 20, 40, 60],
    xlabels=False,
    ylabels=False,
    grid=False,
    fontsize=12,
    color="black",
)
cax2 = plt.axes([0.91, 0.712, 0.013, 0.16])
cart.set_cbar(
    cs2,
    cax2,
    fig,
    orientation="vertical",
    extend="both",
    cbar_label="$months$",
    nbins=7,
    fontsize=12,
    cbar_ticks=[
        np.arange(-np.pi, np.pi + 0.5, (np.pi + np.pi) / 6).tolist(),
        ["Jun", "Aug", "Oct", "Dec", "Feb", "Apr", "Jun"],
    ],
    task="custom ticks",
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
cs3 = ax3.pcolormesh(
    lon,
    lat,
    s_phase_s_m,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=180.0),
)
cart.set_grid_ticks(
    ax3,
    projection=ccrs.PlateCarree(),
    xticks=[0, 45, 90, 135, -360, -180, -135, -90, -45],
    yticks=[-60, -40, -20, 0, 20, 40, 60],
    xlabels=False,
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
cs4 = ax4.pcolormesh(
    lon,
    lat,
    s_phase_w_m,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=180.0),
)
cart.set_grid_ticks(
    ax4,
    projection=ccrs.PlateCarree(),
    xticks=[0, 45, 90, 135, -360, -180, -135, -90, -45],
    yticks=[-60, -40, -20, 0, 20, 40, 60],
    xlabels=False,
    ylabels=False,
    grid=False,
    fontsize=12,
    color="black",
)
cax4 = plt.axes([0.91, 0.517, 0.013, 0.16])
cart.set_cbar(
    cs4,
    cax4,
    fig,
    orientation="vertical",
    extend="both",
    cbar_label="$months$",
    nbins=6,
    fontsize=12,
    cbar_ticks=[
        np.arange(-np.pi, np.pi + 0.5, (np.pi + np.pi) / 5).tolist(),
        ["Jan", "Feb", "Mar", "Apr", "May", "Jun"],
    ],
    task="custom ticks",
)
cart.subplot_label(
    ax4,
    xdist_label=0.1,
    ydist_label=0.88,
    subplot_label="D",
    fs_shade=28,
    form="box",
    fs_main=18,
    color="black",
)

############## Subplot 5  #################
cart.set_subplots(
    ax5, projection, resolution, lon_min=-180, lon_max=179, lat_min=-66, lat_max=66
)
cs5 = ax5.contourf(
    lon,
    lat,
    a_amp_s_m,
    levels=levels,
    cmap=cmo.thermal,
    extend="both",
    transform=ccrs.PlateCarree(central_longitude=180.0),
)
cart.set_grid_ticks(
    ax5,
    projection=ccrs.PlateCarree(),
    xticks=[0, 45, 90, 135, -360, -180, -135, -90, -45],
    yticks=[-60, -40, -20, 0, 20, 40, 60],
    xlabels=False,
    ylabels=True,
    grid=False,
    fontsize=12,
    color="black",
)
cax5 = plt.axes([0.49, 0.18, 0.013, 0.24])
cart.set_cbar(
    cs5,
    cax5,
    fig,
    orientation="vertical",
    extend="both",
    cbar_label="$m$",
    nbins=5,
    fontsize=12,
    cbar_ticks=[],
    task="regular",
)
cart.subplot_label(
    ax5,
    xdist_label=0.1,
    ydist_label=0.88,
    subplot_label="E",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="black",
)

############## Subplot 6  #################
cart.set_subplots(
    ax6, projection, resolution, lon_min=-180, lon_max=179, lat_min=-66, lat_max=66
)
cs6 = ax6.contourf(
    lon,
    lat,
    a_amp_w_m,
    levels=levels,
    cmap=cmo.thermal,
    extend="both",
    transform=ccrs.PlateCarree(central_longitude=180.0),
)
cart.set_grid_ticks(
    ax6,
    projection=ccrs.PlateCarree(),
    xticks=[0, 45, 90, 135, -360, -180, -135, -90, -45],
    yticks=[-60, -40, -20, 0, 20, 40, 60],
    xlabels=False,
    ylabels=False,
    grid=False,
    fontsize=12,
    color="black",
)
cax6 = plt.axes([0.91, 0.18, 0.013, 0.24])
cart.set_cbar(
    cs6,
    cax6,
    fig,
    orientation="vertical",
    extend="both",
    cbar_label="$m\,s^{-1}$",
    nbins=5,
    fontsize=12,
    cbar_ticks=[],
    task="regular",
)
cart.subplot_label(
    ax6,
    xdist_label=0.1,
    ydist_label=0.88,
    subplot_label="F",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="black",
)

############## Subplot 7  #################
cart.set_subplots(
    ax7, projection, resolution, lon_min=-180, lon_max=179, lat_min=-66, lat_max=66
)
cs7 = ax7.contourf(
    lon,
    lat,
    s_amp_s_m,
    levels=levels,
    cmap=cmo.thermal,
    extend="both",
    transform=ccrs.PlateCarree(central_longitude=180.0),
)
cart.set_grid_ticks(
    ax7,
    projection=ccrs.PlateCarree(),
    xticks=[0, 45, 90, 135, -360, -180, -135, -90, -45],
    yticks=[-60, -40, -20, 0, 20, 40, 60],
    xlabels=True,
    ylabels=True,
    grid=False,
    fontsize=12,
    color="black",
)
cart.subplot_label(
    ax7,
    xdist_label=0.1,
    ydist_label=0.88,
    subplot_label="G",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="black",
)

############## Subplot 8  #################
cart.set_subplots(
    ax8, projection, resolution, lon_min=-180, lon_max=179, lat_min=-66, lat_max=66
)
cs8 = ax8.contourf(
    lon,
    lat,
    s_amp_w_m,
    levels=levels,
    cmap=cmo.thermal,
    extend="both",
    transform=ccrs.PlateCarree(central_longitude=180.0),
)
cart.set_grid_ticks(
    ax8,
    projection=ccrs.PlateCarree(),
    xticks=[0, 45, 90, 135, -360, -180, -135, -90, -45 - 0],
    yticks=[-60, -40, -20, 0, 20, 40, 60],
    xlabels=True,
    ylabels=False,
    grid=False,
    fontsize=12,
    color="black",
)
cart.subplot_label(
    ax8,
    xdist_label=0.1,
    ydist_label=0.88,
    subplot_label="H",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="black",
)

# adjust spacing for the entire figure (not the subplot)
fig.subplots_adjust(wspace=0.20, hspace=0.1)

# save figure
fig.savefig(fname="../figs/figS08", bbox_inches="tight", dpi=300)
