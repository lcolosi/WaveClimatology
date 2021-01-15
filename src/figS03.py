"""
Figure S3 Caption: 
Maps of the annual cycle phase of CCMP2 wind speed highlighting SWARs using 
three different criteria: (A) least restrictive, (B) moderately restrictive, and
(C) most restrictive. White pixels correspond to points that are not categorized
as anomalous phase or not statistically significant.
"""

# Path to access python functions
import sys

sys.path.append("../tools/")

# Path to access intermediate data
data_path = "../data/lsf_parameters/"

# libraries
import numpy as np
from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
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

# Mask all phase values outside each of the SWAR condition for each partition
######### Least restrictive case #########
###### Pacific Ocean #######
n_pacific_h_m_sf = np.ma.masked_inside(n_pacific_h, -0.52, 2.09)
n_pacific_l_m_sf = np.ma.masked_inside(n_pacific_l, -0.52, 2.09)

###### Atlantic Ocean #######
n_atlantic_w_h_m_sf = np.ma.masked_inside(n_atlantic_w_h, -0.52, 2.09)
n_atlantic_w_l_m_sf = np.ma.masked_inside(n_atlantic_w_l, -0.52, 2.09)
n_atlantic_e_m_sf = np.ma.masked_inside(n_atlantic_e, -0.52, 2.09)

###### Indian Ocean #######
n_indian_h_m_sf = np.ma.masked_outside(n_indian_h, -1.05, 2.62)
n_indian_l_m_sf = np.ma.masked_outside(n_indian_l, -1.05, 2.62)
m_indian_h_m_sf = np.ma.masked_outside(m_indian_h, -1.05, 2.62)
m_indian_l_m_sf = np.ma.masked_outside(m_indian_l, -1.05, 2.62)
s_indian_m_sf = np.ma.masked_outside(s_indian, -1.05, 2.62)

###### Pacific Ocean #######
s_pacific_m_sf = np.ma.masked_outside(s_pacific, -1.05, 2.62)

###### Atlantic Ocean #######
s_atlantic_w_m_sf = np.ma.masked_outside(s_atlantic_w, -1.05, 2.62)
s_atlantic_e_m_sf = np.ma.masked_outside(s_atlantic_e, -1.05, 2.62)

# Moderately restrictive case
###### Pacific Ocean #######
n_pacific_h_m_msf = np.ma.masked_inside(n_pacific_h, -1.05, 2.62)
n_pacific_l_m_msf = np.ma.masked_inside(n_pacific_l, -1.05, 2.62)

###### Atlantic Ocean #######
n_atlantic_w_h_m_msf = np.ma.masked_inside(n_atlantic_w_h, -1.05, 2.62)
n_atlantic_w_l_m_msf = np.ma.masked_inside(n_atlantic_w_l, -1.05, 2.62)
n_atlantic_e_m_msf = np.ma.masked_inside(n_atlantic_e, -1.05, 2.62)

###### Indian Ocean #######
n_indian_h_m_msf = np.ma.masked_outside(n_indian_h, -0.52, 2.09)
n_indian_l_m_msf = np.ma.masked_outside(n_indian_l, -0.52, 2.09)
m_indian_h_m_msf = np.ma.masked_outside(m_indian_h, -0.52, 2.09)
m_indian_l_m_msf = np.ma.masked_outside(m_indian_l, -0.52, 2.09)
s_indian_m_msf = np.ma.masked_outside(s_indian, -0.52, 2.09)

###### Pacific Ocean #######
s_pacific_m_msf = np.ma.masked_outside(s_pacific, -0.52, 2.09)

###### Atlantic Ocean #######
s_atlantic_w_m_msf = np.ma.masked_outside(s_atlantic_w, -0.52, 2.09)
s_atlantic_e_m_msf = np.ma.masked_outside(s_atlantic_e, -0.52, 2.09)

######### Most restrictive case #########
###### Pacific Ocean #######
n_pacific_h_m_sum = np.ma.masked_greater_equal(n_pacific_h, -np.pi / 2)
n_pacific_l_m_sum = np.ma.masked_greater_equal(n_pacific_l, -np.pi / 2)

###### Atlantic Ocean #######
n_atlantic_w_h_m_sum = np.ma.masked_greater_equal(n_atlantic_w_h, -np.pi / 2)
n_atlantic_w_l_m_sum = np.ma.masked_greater_equal(n_atlantic_w_l, -np.pi / 2)
n_atlantic_e_m_sum = np.ma.masked_greater_equal(n_atlantic_e, -np.pi / 2)

###### Indian Ocean #######
n_indian_h_m_sum = np.ma.masked_outside(n_indian_h, 0, (np.pi / 2))
n_indian_l_m_sum = np.ma.masked_outside(n_indian_l, 0, (np.pi / 2))
m_indian_h_m_sum = np.ma.masked_outside(m_indian_h, 0, (np.pi / 2))
m_indian_l_m_sum = np.ma.masked_outside(m_indian_l, 0, (np.pi / 2))
s_indian_m_sum = np.ma.masked_outside(s_indian, 0, (np.pi / 2))

###### Pacific Ocean #######
s_pacific_m_sum = np.ma.masked_outside(s_pacific, 0, (np.pi / 2))

###### Atlantic Ocean #######
s_atlantic_w_m_sum = np.ma.masked_outside(s_atlantic_w, 0, (np.pi / 2))
s_atlantic_e_m_sum = np.ma.masked_outside(s_atlantic_e, 0, (np.pi / 2))

# Initialize variables for plotting
projection = ccrs.PlateCarree(central_longitude=180.0)
xticks = [0, 30, 60, 90, 120, 150, -360, -180, -150, -120, -90, -60, -30, -0]
yticks = [-60, -45, -30, -15, 0, 15, 30, 45, 60]
resolution = "110m"

# Create figure and axes
fig, axes = plt.subplots(3, 1, figsize=(24, 20), subplot_kw={"projection": projection})
ax1, ax2, ax3 = axes.flatten()

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

############## Least Restrictive  #################
cart.set_subplots(
    ax1, projection, resolution, lon_min=-180, lon_max=179, lat_min=-66, lat_max=66
)

# North indian ocean high
cs1 = ax1.pcolor(
    lon_reg_nih,
    lat_reg_nih,
    n_indian_h_m_sf,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# North indian ocean low
cs2 = ax1.pcolor(
    lon_reg_nil,
    lat_reg_nil,
    n_indian_l_m_sf,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# Mid indian ocean high
cs3 = ax1.pcolor(
    lon_reg_mih,
    lat_reg_mih,
    m_indian_h_m_sf,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# Mid indian Ocean low
cs4 = ax1.pcolor(
    lon_reg_mil,
    lat_reg_mil,
    m_indian_l_m_sf,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# South Indian Ocean
cs5 = ax1.pcolor(
    lon_reg_si,
    lat_reg_si,
    s_indian_m_sf,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# north pacific high
cs6 = ax1.pcolor(
    lon_reg_nph,
    lat_reg_nph,
    n_pacific_h_m_sf,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# north pacific low
cs7 = ax1.pcolor(
    lon_reg_npl,
    lat_reg_npl,
    n_pacific_l_m_sf,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# south pacific
cs8 = ax1.pcolor(
    lon_reg_sp,
    lat_reg_sp,
    s_pacific_m_sf,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# North Atlantic West high
cs9 = ax1.pcolor(
    lon_reg_nawh,
    lat_reg_nawh,
    n_atlantic_w_h_m_sf,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# North Atlantic West low
cs10 = ax1.pcolor(
    lon_reg_nawl,
    lat_reg_nawl,
    n_atlantic_w_l_m_sf,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# North Atlantic East
cs11 = ax1.pcolor(
    lon_reg_nae,
    lat_reg_nae,
    n_atlantic_e_m_sf,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# South Atlantic West
cs12 = ax1.pcolor(
    lon_reg_saw,
    lat_reg_saw,
    s_atlantic_w_m_sf,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# South Atlantic East
cs13 = ax1.pcolor(
    lon_reg_sae,
    lat_reg_sae,
    s_atlantic_e_m_sf,
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
    grid=True,
    fontsize=20,
    color="black",
)
cart.subplot_label(
    ax1,
    xdist_label=0.1,
    ydist_label=0.88,
    subplot_label="A",
    form="box",
    fs_shade=35,
    fs_main=25,
    color=None,
)

############## Moderately Restrictive #################
cart.set_subplots(
    ax2, projection, resolution, lon_min=-180, lon_max=179, lat_min=-66, lat_max=66
)

# North indian ocean high
cs1 = ax2.pcolor(
    lon_reg_nih,
    lat_reg_nih,
    n_indian_h_m_msf,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# North indian ocean low
cs2 = ax2.pcolor(
    lon_reg_nil,
    lat_reg_nil,
    n_indian_l_m_msf,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# Mid indian ocean high
cs3 = ax2.pcolor(
    lon_reg_mih,
    lat_reg_mih,
    m_indian_h_m_msf,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# Mid indian Ocean low
cs4 = ax2.pcolor(
    lon_reg_mil,
    lat_reg_mil,
    m_indian_l_m_msf,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# South Indian Ocean
cs5 = ax2.pcolor(
    lon_reg_si,
    lat_reg_si,
    s_indian_m_msf,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# north pacific high
cs6 = ax2.pcolor(
    lon_reg_nph,
    lat_reg_nph,
    n_pacific_h_m_msf,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# north pacific low
cs7 = ax2.pcolor(
    lon_reg_npl,
    lat_reg_npl,
    n_pacific_l_m_msf,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# south pacific
cs8 = ax2.pcolor(
    lon_reg_sp,
    lat_reg_sp,
    s_pacific_m_msf,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# North Atlantic West high
cs9 = ax2.pcolor(
    lon_reg_nawh,
    lat_reg_nawh,
    n_atlantic_w_h_m_msf,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# North Atlantic West low
cs10 = ax2.pcolor(
    lon_reg_nawl,
    lat_reg_nawl,
    n_atlantic_w_l_m_msf,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# North Atlantic East
cs11 = ax2.pcolor(
    lon_reg_nae,
    lat_reg_nae,
    n_atlantic_e_m_msf,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# South Atlantic West
cs12 = ax2.pcolor(
    lon_reg_saw,
    lat_reg_saw,
    s_atlantic_w_m_msf,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# South Atlantic East
cs13 = ax2.pcolor(
    lon_reg_sae,
    lat_reg_sae,
    s_atlantic_e_m_msf,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)

cart.set_grid_ticks(
    ax2,
    projection=ccrs.PlateCarree(),
    xticks=xticks,
    yticks=yticks,
    xlabels=False,
    ylabels=True,
    grid=True,
    fontsize=20,
    color="black",
)
cart.subplot_label(
    ax2,
    xdist_label=0.1,
    ydist_label=0.88,
    subplot_label="B",
    form="box",
    fs_shade=35,
    fs_main=25,
    color=None,
)

############## Most Restrictive #################
cart.set_subplots(
    ax3, projection, resolution, lon_min=-180, lon_max=179, lat_min=-66, lat_max=66
)

# North indian ocean high
cs1 = ax3.pcolor(
    lon_reg_nih,
    lat_reg_nih,
    n_indian_h_m_sum,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# North indian ocean low
cs2 = ax3.pcolor(
    lon_reg_nil,
    lat_reg_nil,
    n_indian_l_m_sum,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# Mid indian ocean high
cs3 = ax3.pcolor(
    lon_reg_mih,
    lat_reg_mih,
    m_indian_h_m_sum,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# Mid indian Ocean low
cs4 = ax3.pcolor(
    lon_reg_mil,
    lat_reg_mil,
    m_indian_l_m_sum,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# South Indian Ocean
cs5 = ax3.pcolor(
    lon_reg_si,
    lat_reg_si,
    s_indian_m_sum,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# north pacific high
cs6 = ax3.pcolor(
    lon_reg_nph,
    lat_reg_nph,
    n_pacific_h_m_sum,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# north pacific low
cs7 = ax3.pcolor(
    lon_reg_npl,
    lat_reg_npl,
    n_pacific_l_m_sum,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# south pacific
cs8 = ax3.pcolor(
    lon_reg_sp,
    lat_reg_sp,
    s_pacific_m_sum,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# North Atlantic West high
cs9 = ax3.pcolor(
    lon_reg_nawh,
    lat_reg_nawh,
    n_atlantic_w_h_m_sum,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# North Atlantic West low
cs10 = ax3.pcolor(
    lon_reg_nawl,
    lat_reg_nawl,
    n_atlantic_w_l_m_sum,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# North Atlantic East
cs11 = ax3.pcolor(
    lon_reg_nae,
    lat_reg_nae,
    n_atlantic_e_m_sum,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# South Atlantic West
cs12 = ax3.pcolor(
    lon_reg_saw,
    lat_reg_saw,
    s_atlantic_w_m_sum,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
# South Atlantic East
cs13 = ax3.pcolor(
    lon_reg_sae,
    lat_reg_sae,
    s_atlantic_e_m_sum,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)

cart.set_grid_ticks(
    ax3,
    projection=ccrs.PlateCarree(),
    xticks=xticks,
    yticks=yticks,
    xlabels=True,
    ylabels=True,
    grid=True,
    fontsize=20,
    color="black",
)
cart.subplot_label(
    ax3,
    xdist_label=0.1,
    ydist_label=0.88,
    subplot_label="C",
    form="box",
    fs_shade=35,
    fs_main=25,
    color=None,
)

# Create Discrete Phase map
ax4 = fig.add_axes([0.80, 0.35, 0.02, 0.3])
cb = matplotlib.colorbar.ColorbarBase(
    ax4,
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

# adjust spacing for the entire figure
fig.subplots_adjust(wspace=0.35, hspace=0.05)

# save figure:
plt.savefig(fname="../figs/figS03", bbox_inches="tight", dpi=300)
