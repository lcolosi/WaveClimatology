"""
Figure 5 Caption: 
Difference between the annual cycle phases, $\phi_\mathrm{wsp}$ and 
$\phi_\mathrm{swh}$. White pixels correspond to points where phase differences 
are not statistically significant, or where the annual cycle amplitude is small
enough that phase is not well defined.
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
filename_ifremer = data_path + "IFREMER_swh_lsf_parameters.nc"
filename_ccmp2 = data_path + "CCMP2_wsp_lsf_parameters.nc"

######### Ifremer #########
nc_i = Dataset(filename_ifremer, "r")
lon = nc_i.variables["lon"][:]
lat = nc_i.variables["lat"][:]
a_amp_i = nc_i.variables["a_amp"][:]
a_phase_i = nc_i.variables["a_phase"][:]
a_amp_unc_i = nc_i.variables["a_amp_unc"][:]

######### CCMP2 #########
nc_c = Dataset(filename_ccmp2, "r")
a_amp_c = nc_c.variables["a_amp"][:]
a_phase_c = nc_c.variables["a_phase"][:]
a_amp_unc_c = nc_c.variables["a_amp_unc"][:]

# Set noise to signal ratio criteria
ns = 5 / 10

# Mask not statistical significance ifremer and ccmp2 data
######### Ifremer #########
# Compute the relative uncertainty
a_ratio_i = a_amp_unc_i / a_amp_i
# Mask not statistically significant grid points
a_mask_i = np.ma.getmask(np.ma.masked_greater_equal(a_ratio_i, (ns)))
# Apply statistical significance masks
a_phase_i_m = np.ma.masked_where(a_mask_i, a_phase_i)

######### CCMP2 #########
# Compute the relative uncertainty
a_ratio_c = a_amp_unc_c / a_amp_c
# Mask not statistically significant grid points
a_mask_c = np.ma.getmask(np.ma.masked_greater_equal(a_ratio_c, (ns)))
# Apply statistical significance masks
a_phase_c_m = np.ma.masked_where(a_mask_c, a_phase_c)

# Compute the difference between wsp and swh annual cycle phase
phase_diff = a_phase_c_m - a_phase_i_m

# For the phase to range from +- 6 months, subtract (add) 2pi and multiply by -1 for all phase values greater (less) than pi (-pi).
phase_diff[np.where(phase_diff > np.pi)] = -1 * (
    phase_diff[np.where(phase_diff > np.pi)] - 2 * np.pi
)
phase_diff[np.where(phase_diff < -np.pi)] = -1 * (
    phase_diff[np.where(phase_diff < -np.pi)] + 2 * np.pi
)

# Initialize variables for plotting
projection = ccrs.PlateCarree(central_longitude=180.0)
levels = np.linspace(-np.pi, np.pi, 50)
xticks = [0, 45, 90, 135, -360, -180, -135, -90, -45]
yticks = [-60, -40, -20, 0, 20, 40, 60]
resolution = "110m"

# Create figure
fig, axes = plt.subplots(1, 1, figsize=(24, 20), subplot_kw={"projection": projection})

# Initialize variables for discrete phase map
# Specify colormap
cmap = cmo.diff

# extract all colors from the phase map
cmaplist = [cmap(i) for i in range(cmap.N)]

# create the new map
cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
    "Custom cmap", cmaplist, cmap.N
)

# define the bins and normalize
bounds = np.linspace(-np.pi, np.pi, 13)
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

######## Suplot 1 ########
cart.set_subplots(
    axes,
    projection,
    resolution,
    lon_min=-180,
    lon_max=179,
    lat_min=-66,
    lat_max=66,
)
cs1 = axes.pcolor(
    lon,
    lat,
    phase_diff,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmap,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
cart.set_grid_ticks(
    axes,
    projection=ccrs.PlateCarree(),
    xticks=xticks,
    yticks=yticks,
    xlabels=True,
    ylabels=True,
    grid=False,
    fontsize=22,
    color="black",
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
    extend="neither",
    extendfrac="auto",
    orientation="vertical",
)
cb.ax.set_ylabel("$months$", fontsize=22)
cb.ax.set_yticklabels(np.arange(-6, 7, 1).tolist())
cb.ax.tick_params(labelsize=22)

# save figure
fig.savefig(fname="../figs/fig05", bbox_inches="tight", dpi=300)
