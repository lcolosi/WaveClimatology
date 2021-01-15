"""
Figure 1 caption: 
DJF seasonal mean for (A) IFREMER SWH and (B) CCMP2 WSP; JJA seasonal
mean for (C) IFREMER SWH and (D) CCMP2 WSP; DJF standard deviation of daily
data for (E) IFREMER SWH and (F) CCMP2 WSP; JJA standard deviation of daily
data for (G) IFREMER SWH and (H) CCMP2 WSP.
"""

# Path to access python functions
import sys

sys.path.append("../tools/")

# Path to access intermediate data
data_path_i = "../data/ifremer_swh/"
data_path_c = "../data/ccmp2_wsp/"

# libraries
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cmocean.cm as cmo
import matplotlib.patches as mpatches

# my functions
from intermediate_data_import import import_data
from statistical_moments_temporal import stat_moments_temporal
import cartopy_figs as cart

# call IFREMER SWH and CCMP2 WSP processed data:
swh, time_s, lat_s, lon_s = import_data("IFREMER_swh", data_path_i)
wsp, time_w, lat_w, lon_w = import_data("CCMP2_wsp", data_path_c)

# Compute statistical moments seasonally
swh_stats_s = stat_moments_temporal(swh, time_s, "seasonally", "sample")
wsp_stats_s = stat_moments_temporal(wsp, time_w, "seasonally", "sample")

# Set mean and varaince statistics:
swh_mean = np.ma.array(swh_stats_s["mean"])
swh_variance = np.ma.array(swh_stats_s["var"])
wsp_mean = np.ma.array(wsp_stats_s["mean"])
wsp_variance = np.ma.array(wsp_stats_s["var"])

# Compute standard deviation
swh_std = np.ma.sqrt(swh_variance)
wsp_std = np.ma.sqrt(wsp_variance)

# Initialize variables for plotting
projection = ccrs.PlateCarree(central_longitude=180.0)
levels_sm = np.arange(0, 4.6, 0.1)
levels_wm = np.arange(0, 12.1, 0.1)
levels_sv = np.arange(0.0, 2.1, 0.1)
levels_wv = np.arange(0, 4.1, 0.1)
resolution = "50m"

# Create figure and axes
fig, axes = plt.subplots(4, 2, figsize=(16, 12), subplot_kw={"projection": projection})
ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8 = axes.flatten()

############## Subplot 1  #################
cart.set_subplots(
    ax1, projection, resolution, lon_min=-180, lon_max=179, lat_min=-66, lat_max=66
)
cs1 = ax1.contourf(
    lon_s,
    lat_s,
    swh_mean[0, :, :],
    levels_sm,
    cmap=cmo.haline,
    extend="both",
    transform=ccrs.PlateCarree(central_longitude=0.0),
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
cax1 = plt.axes([0.49, 0.57, 0.013, 0.24])
cart.set_cbar(
    cs1,
    cax1,
    fig,
    orientation="vertical",
    extend="both",
    cbar_label="$m$",
    nbins=9,
    fontsize=12,
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
    xdist_label=0.22,
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
    lon_w,
    lat_w,
    wsp_mean[0, :, :],
    levels_wm,
    cmap=cmo.thermal,
    extend="both",
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
cart.set_grid_ticks(
    ax2,
    projection=ccrs.PlateCarree(),
    xticks=[0, 45, 90, 135, -360, -180, -135, -90, -45 - 0],
    yticks=[-60, -40, -20, 0, 20, 40, 60],
    xlabels=False,
    ylabels=False,
    grid=False,
    fontsize=12,
    color="black",
)
cax2 = plt.axes([0.91, 0.57, 0.013, 0.24])
cart.set_cbar(
    cs2,
    cax2,
    fig,
    orientation="vertical",
    extend="both",
    cbar_label="$m\,{s}^{-1}$",
    nbins=7,
    fontsize=12,
    cbar_ticks=[],
    task="regular",
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
    xdist_label=0.22,
    ydist_label=0.86,
    subplot_label="DJF",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="black",
)

# Set regional climatology grid boxes for SH
ax2.add_patch(
    mpatches.Rectangle(
        xy=[108, -31],
        width=4,
        height=4,
        edgecolor="black",
        facecolor=None,
        linestyle="-",
        linewidth=1,
        fill=False,
        alpha=1,
        transform=ccrs.PlateCarree(),
    )
)
ax2.add_patch(
    mpatches.Rectangle(
        xy=[-78, -36],
        width=4,
        height=4,
        edgecolor="black",
        facecolor=None,
        linestyle="-",
        linewidth=1,
        fill=False,
        alpha=1,
        transform=ccrs.PlateCarree(),
    )
)
ax2.add_patch(
    mpatches.Rectangle(
        xy=[10, -30],
        width=4,
        height=4,
        edgecolor="black",
        facecolor=None,
        linestyle="-",
        linewidth=1,
        fill=False,
        alpha=1,
        transform=ccrs.PlateCarree(),
    )
)
ax2.add_patch(
    mpatches.Rectangle(
        xy=[51, 4],
        width=4,
        height=4,
        edgecolor="black",
        facecolor=None,
        linestyle="-",
        linewidth=1,
        fill=False,
        alpha=1,
        transform=ccrs.PlateCarree(),
    )
)

# set labels for grid boxes
cart.subplot_label(
    ax2,
    xdist_label=0.37,
    ydist_label=0.31,
    subplot_label="7A",
    form="box",
    fs_shade=12,
    fs_main=12,
    color="tab:blue",
)
cart.subplot_label(
    ax2,
    xdist_label=0.84,
    ydist_label=0.31,
    subplot_label="7C",
    form="box",
    fs_shade=12,
    fs_main=12,
    color="tab:blue",
)
cart.subplot_label(
    ax2,
    xdist_label=0.07,
    ydist_label=0.37,
    subplot_label="7E",
    form="box",
    fs_shade=12,
    fs_main=12,
    color="tab:blue",
)
cart.subplot_label(
    ax2,
    xdist_label=0.1,
    ydist_label=0.61,
    subplot_label="7G",
    form="box",
    fs_shade=12,
    fs_main=12,
    color="tab:blue",
)

############## Subplot 3  #################
cart.set_subplots(
    ax3, projection, resolution, lon_min=-180, lon_max=179, lat_min=-66, lat_max=66
)
cs3 = ax3.contourf(
    lon_s,
    lat_s,
    swh_mean[2, :, :],
    levels_sm,
    cmap=cmo.haline,
    extend="both",
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
cart.set_grid_ticks(
    ax3,
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
    xdist_label=0.22,
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
    lon_w,
    lat_w,
    wsp_mean[2, :, :],
    levels_wm,
    cmap=cmo.thermal,
    extend="both",
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
cart.set_grid_ticks(
    ax4,
    projection=ccrs.PlateCarree(),
    xticks=[0, 45, 90, 135, -360, -180, -135, -90, -45 - 0],
    yticks=[-60, -40, -20, 0, 20, 40, 60],
    xlabels=False,
    ylabels=False,
    grid=False,
    fontsize=12,
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
    xdist_label=0.22,
    ydist_label=0.86,
    subplot_label="JJA",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="black",
)

# Set regional climatology grid boxes for NH
ax4.add_patch(
    mpatches.Rectangle(
        xy=[-129, 36],
        width=4,
        height=4,
        edgecolor="black",
        facecolor=None,
        linestyle="-",
        linewidth=1,
        fill=False,
        alpha=1,
        transform=ccrs.PlateCarree(),
    )
)
ax4.add_patch(
    mpatches.Rectangle(
        xy=[-76, 13],
        width=4,
        height=4,
        edgecolor="black",
        facecolor=None,
        linestyle="-",
        linewidth=1,
        fill=False,
        alpha=1,
        transform=ccrs.PlateCarree(),
    )
)
ax4.add_patch(
    mpatches.Rectangle(
        xy=[-15, 30],
        width=4,
        height=4,
        edgecolor="black",
        facecolor=None,
        linestyle="-",
        linewidth=1,
        fill=False,
        alpha=1,
        transform=ccrs.PlateCarree(),
    )
)

# set labels for grid boxes
cart.subplot_label(
    ax4,
    xdist_label=0.70,
    ydist_label=0.83,
    subplot_label="6A",
    form="box",
    fs_shade=12,
    fs_main=12,
    color="tab:blue",
)
cart.subplot_label(
    ax4,
    xdist_label=0.82,
    ydist_label=0.50,
    subplot_label="6C",
    form="box",
    fs_shade=12,
    fs_main=12,
    color="tab:blue",
)
cart.subplot_label(
    ax4,
    xdist_label=0.97,
    ydist_label=0.62,
    subplot_label="6E",
    form="box",
    fs_shade=12,
    fs_main=12,
    color="tab:blue",
)

############## Subplot 5  #################
cart.set_subplots(
    ax5, projection, resolution, lon_min=-180, lon_max=179, lat_min=-66, lat_max=66
)
cs5 = ax5.contourf(
    lon_s,
    lat_s,
    swh_std[0, :, :],
    levels_sv,
    cmap=cmo.haline,
    extend="both",
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
cart.set_grid_ticks(
    ax5,
    projection=ccrs.PlateCarree(),
    xticks=[0, 45, 90, 135, -360, -180, -135, -90, -45 - 0],
    yticks=[-60, -40, -20, 0, 20, 40, 60],
    xlabels=False,
    ylabels=True,
    grid=False,
    fontsize=12,
    color="black",
)
cax3 = plt.axes([0.49, 0.19, 0.013, 0.24])
cart.set_cbar(
    cs5,
    cax3,
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
    ydist_label=0.86,
    subplot_label="E",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="black",
)
cart.subplot_label(
    ax5,
    xdist_label=0.22,
    ydist_label=0.86,
    subplot_label="DJF",
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
    lon_w,
    lat_w,
    wsp_std[0, :, :],
    levels_wv,
    cmap=cmo.thermal,
    extend="both",
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
cart.set_grid_ticks(
    ax6,
    projection=ccrs.PlateCarree(),
    xticks=[0, 45, 90, 135, -360, -180, -135, -90, -45 - 0],
    yticks=[-60, -40, -20, 0, 20, 40, 60],
    xlabels=False,
    ylabels=False,
    grid=False,
    fontsize=12,
    color="black",
)
cax4 = plt.axes([0.91, 0.19, 0.013, 0.24])
cart.set_cbar(
    cs6,
    cax4,
    fig,
    orientation="vertical",
    extend="both",
    cbar_label="$m\,{s}^{-1}$",
    nbins=6,
    fontsize=12,
    cbar_ticks=[],
    task="regular",
)
cart.subplot_label(
    ax6,
    xdist_label=0.1,
    ydist_label=0.86,
    subplot_label="F",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="black",
)
cart.subplot_label(
    ax6,
    xdist_label=0.22,
    ydist_label=0.86,
    subplot_label="DJF",
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
    lon_s,
    lat_s,
    swh_std[2, :, :],
    levels_sv,
    cmap=cmo.haline,
    extend="both",
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
cart.set_grid_ticks(
    ax7,
    projection=ccrs.PlateCarree(),
    xticks=[0, 45, 90, 135, -360, -180, -135, -90, -45 - 0],
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
    ydist_label=0.86,
    subplot_label="G",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="black",
)
cart.subplot_label(
    ax7,
    xdist_label=0.2,
    ydist_label=0.86,
    subplot_label="JJA",
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
    lon_w,
    lat_w,
    wsp_std[2, :, :],
    levels_wv,
    cmap=cmo.thermal,
    extend="both",
    transform=ccrs.PlateCarree(central_longitude=0.0),
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
    ydist_label=0.86,
    subplot_label="H",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="black",
)
cart.subplot_label(
    ax8,
    xdist_label=0.2,
    ydist_label=0.86,
    subplot_label="JJA",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="black",
)

# adjust spacing for figure
fig.subplots_adjust(wspace=0.19, hspace=0.1)

# save figure
fig.savefig(fname="../figs/fig01", bbox_inches="tight", dpi=300)
