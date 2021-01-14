"""
Figure S5 Caption: 
Same as Figure~4 from main text with geographic locations of regional 
climatologies for SWARs in Figures~6,~7.
"""

# Path to access python functions
import sys

sys.path.append("../tools/")

# Path to access intermediate data
data_path = "../data/"

# libraries
import numpy as np
from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import cartopy.crs as ccrs
import cmocean.cm as cmo

# my functions
import cartopy_figs as cart

# Call data:
filename_ccmp2 = data_path + "CCMP2_wsp_lsf_parameters.nc"

######### CCMP2 #########
nc_c = Dataset(filename_ccmp2, "r")
lon = nc_c.variables["lon"][:]
lat = nc_c.variables["lat"][:]
a_amp_c = nc_c.variables["a_amp"][:]
a_phase_c = nc_c.variables["a_phase"][:]
a_amp_unc_c = nc_c.variables["a_amp_unc"][:]

# Set noise to signal ratio criteria
ns = 5 / 10

# Mask not statistical significance ccmp2 data
# Compute the relative uncertainty
a_ratio_c = a_amp_unc_c / a_amp_c
# Mask not statistically significant grid points
a_mask_c = np.ma.getmask(np.ma.masked_greater_equal(a_ratio_c, (ns)))
# Apply statistical significance masks
a_phase_c_m = np.ma.masked_where(a_mask_c, a_phase_c)

# Partition the world oceans excluding the equatorial regions and marginal seas:
n_indian_h = a_phase_c_m[89:95, 60:105]
n_indian_l = a_phase_c_m[77:90, 51:105]
m_indian_h = a_phase_c_m[69:78, 40:98]
m_indian_l = a_phase_c_m[59:70, 40:105]
s_indian = a_phase_c_m[0:60, 23:132]
n_pacific_h = a_phase_c_m[116:133, 106:266]
n_pacific_l = a_phase_c_m[75:117, 106:275]
s_pacific = a_phase_c_m[0:60, 131:292]
n_atlantic_w_h = a_phase_c_m[116:133, 295:361]
n_atlantic_w_l = a_phase_c_m[77:117, 274:361]
n_atlantic_e = a_phase_c_m[77:133, 0:24]
s_atlantic_w = a_phase_c_m[0:68, 291:360]
s_atlantic_e = a_phase_c_m[0:68, 0:24]

# Initialize longitude and latitude variables for each partition:
# North indian ocean high
lat_reg_nih = lat[89:95]
lon_reg_nih = lon[60:105]

# North indian ocean low
lat_reg_nil = lat[77:90]
lon_reg_nil = lon[51:105]

# Mid indian ocean high
lat_reg_mih = lat[69:78]
lon_reg_mih = lon[40:98]

# Mid indian Ocean low
lat_reg_mil = lat[59:70]
lon_reg_mil = lon[40:105]

# South Indian Ocean
lat_reg_si = lat[0:60]
lon_reg_si = lon[23:132]

# north pacific high
lat_reg_nph = lat[116:133]
lon_reg_nph = lon[106:266]

# north pacific low
lat_reg_npl = lat[75:117]
lon_reg_npl = lon[106:275]

# south pacific
lat_reg_sp = lat[0:60]
lon_reg_sp = lon[131:292]

# North Atlantic West high
lat_reg_nawh = lat[116:133]
lon_reg_nawh = lon[295:360]

# North Atlantic West low
lat_reg_nawl = lat[77:117]
lon_reg_nawl = lon[274:360]

# North Atlantic East
lat_reg_nae = lat[77:133]
lon_reg_nae = lon[0:24]

# South Atlantic West
lat_reg_saw = lat[0:68]
lon_reg_saw = lon[291:360]

# South Atlantic East
lat_reg_sae = lat[0:68]
lon_reg_sae = lon[0:24]

# Masking all phase values outside the SWAR condition for each partition
############# Northern Hemisphere #############
###### Pacific Ocean #######
n_pacific_h_m = np.ma.masked_inside(n_pacific_h, -0.52, 2.09)
n_pacific_l_m = np.ma.masked_inside(n_pacific_l, -0.52, 2.09)

###### Atlantic Ocean #######
n_atlantic_w_h_m = np.ma.masked_inside(n_atlantic_w_h, -0.52, 2.09)
n_atlantic_w_l_m = np.ma.masked_inside(n_atlantic_w_l, -0.52, 2.09)
n_atlantic_e_m = np.ma.masked_inside(n_atlantic_e, -0.52, 2.09)

############# Southern Hemisphere #############
###### Indian Ocean #######
n_indian_h_m = np.ma.masked_outside(n_indian_h, -1.05, 2.62)
n_indian_l_m = np.ma.masked_outside(n_indian_l, -1.05, 2.62)
m_indian_h_m = np.ma.masked_outside(m_indian_h, -1.05, 2.62)
m_indian_l_m = np.ma.masked_outside(m_indian_l, -1.05, 2.62)
s_indian_m = np.ma.masked_outside(s_indian, -1.05, 2.62)

###### Pacific Ocean #######
s_pacific_m = np.ma.masked_outside(s_pacific, -1.05, 2.62)

###### Atlantic Ocean #######
s_atlantic_w_m = np.ma.masked_outside(s_atlantic_w, -1.05, 2.62)
s_atlantic_e_m = np.ma.masked_outside(s_atlantic_e, -1.05, 2.62)

# Initialize variables for plotting
projection = ccrs.PlateCarree(central_longitude=180.0)
xticks = [0, 30, 60, 90, 120, 150, -180, -150, -120, -90, -60, -30, -0]
yticks = [-60, -45, -30, -15, 0, 15, 30, 45, 60]
resolution = "50m"

# Create figure
fig, axes = plt.subplots(1, 1, figsize=(24, 20), subplot_kw={"projection": projection})

# Initialize variables for discrete phase map
# Specify colormap
cmap = cmo.phase

# extract all colors from the phase map
cmaplist = [cmap(i) for i in range(cmap.N)]

# create the new map
cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
    "Custom cmap", cmaplist, cmap.N
)

# define the bins and normalize
bounds = np.linspace(-np.pi, np.pi, 15)
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

######## Suplot 1 ########
cart.set_subplots(
    axes, projection, resolution, lon_min=-180, lon_max=179, lat_min=-66, lat_max=66
)

# North indian ocean high
cs1 = axes.pcolormesh(
    lon_reg_nih,
    lat_reg_nih,
    n_indian_h_m,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# North indian ocean low
cs2 = axes.pcolormesh(
    lon_reg_nil,
    lat_reg_nil,
    n_indian_l_m,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# Mid indian ocean high
cs3 = axes.pcolormesh(
    lon_reg_mih,
    lat_reg_mih,
    m_indian_h_m,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# Mid indian Ocean low
cs4 = axes.pcolormesh(
    lon_reg_mil,
    lat_reg_mil,
    m_indian_l_m,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# South Indian Ocean
cs5 = axes.pcolormesh(
    lon_reg_si,
    lat_reg_si,
    s_indian_m,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# North pacific high
cs6 = axes.pcolormesh(
    lon_reg_nph,
    lat_reg_nph,
    n_pacific_h_m,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# North pacific low
cs7 = axes.pcolormesh(
    lon_reg_npl,
    lat_reg_npl,
    n_pacific_l_m,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# South pacific
cs8 = axes.pcolormesh(
    lon_reg_sp,
    lat_reg_sp,
    s_pacific_m,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# North Atlantic West high
cs9 = axes.pcolormesh(
    lon_reg_nawh,
    lat_reg_nawh,
    n_atlantic_w_h_m,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# North Atlantic West low
cs10 = axes.pcolormesh(
    lon_reg_nawl,
    lat_reg_nawl,
    n_atlantic_w_l_m,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# North Atlantic East
cs11 = axes.pcolormesh(
    lon_reg_nae,
    lat_reg_nae,
    n_atlantic_e_m,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# South Atlantic West
cs12 = axes.pcolormesh(
    lon_reg_saw,
    lat_reg_saw,
    s_atlantic_w_m,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# South Atlantic East
cs13 = axes.pcolormesh(
    lon_reg_sae,
    lat_reg_sae,
    s_atlantic_e_m,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)

# Create Discrete Phase map
ax2 = fig.add_axes([0.93, 0.35, 0.02, 0.3])
cb = matplotlib.colorbar.ColorbarBase(
    ax2,
    cmap=cmap,
    norm=norm,
    spacing="proportional",
    ticks=bounds,
    boundaries=bounds,
    extend="both",
    extendfrac="auto",
    orientation="vertical",
)
cb.ax.set_ylabel("$months$", fontsize=16)
cb.ax.set_yticklabels(
    [
        "Jun",
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
)
cb.ax.tick_params(labelsize=15)

# Plot regional climatology grid boxes
axes.add_patch(
    mpatches.Rectangle(
        xy=[-165, 20],
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
axes.add_patch(
    mpatches.Rectangle(
        xy=[-104, 12],
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
axes.add_patch(
    mpatches.Rectangle(
        xy=[-58, 17],
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
axes.add_patch(
    mpatches.Rectangle(
        xy=[-10, -4],
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
axes.add_patch(
    mpatches.Rectangle(
        xy=[6, -8],
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
axes.add_patch(
    mpatches.Rectangle(
        xy=[37, -27],
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
axes.add_patch(
    mpatches.Rectangle(
        xy=[-123, -20],
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
axes.add_patch(
    mpatches.Rectangle(
        xy=[67, -1],
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
axes.add_patch(
    mpatches.Rectangle(
        xy=[155, -26],
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

# Plot labels for grid boxes
cart.subplot_label(
    axes,
    xdist_label=0.515,
    ydist_label=0.71,
    subplot_label="S6A",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="tab:blue",
)
cart.subplot_label(
    axes,
    xdist_label=0.69,
    ydist_label=0.57,
    subplot_label="S6B",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="tab:blue",
)
cart.subplot_label(
    axes,
    xdist_label=0.83,
    ydist_label=0.70,
    subplot_label="S6C",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="tab:blue",
)
cart.subplot_label(
    axes,
    xdist_label=0.96,
    ydist_label=0.44,
    subplot_label="S7A",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="tab:blue",
)
cart.subplot_label(
    axes,
    xdist_label=0.05,
    ydist_label=0.510,
    subplot_label="S7B",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="tab:blue",
)
cart.subplot_label(
    axes,
    xdist_label=0.12,
    ydist_label=0.25,
    subplot_label="S7C",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="tab:blue",
)
cart.subplot_label(
    axes,
    xdist_label=0.65,
    ydist_label=0.30,
    subplot_label="S7D",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="tab:blue",
)
cart.subplot_label(
    axes,
    xdist_label=0.165,
    ydist_label=0.47,
    subplot_label="S7E",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="tab:blue",
)
cart.subplot_label(
    axes,
    xdist_label=0.45,
    ydist_label=0.25,
    subplot_label="S7F",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="tab:blue",
)

cart.set_grid_ticks(
    axes,
    projection=ccrs.PlateCarree(),
    xticks=xticks,
    yticks=yticks,
    xlabels=True,
    ylabels=True,
    grid=True,
    fontsize=22,
    color="black",
)

# save figure
fig.savefig(fname="../figs/figS05", bbox_inches="tight", dpi=300)
