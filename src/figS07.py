"""
Figure S7 Caption: 
(left column) Southern Hemisphere highlighted SWARs in Figure S05 with (right 
column) IFREMER SWH (solid blue) and CCMP2 WSP (solid red) climatologies 
extracted from the outlined 4$^{\circ}$ by 4$^{\circ}$ boxes within SWARs. 
Shading, dotted lines, and solid black are as in Figure S6. SWARs include 
Central South Atlantic Ocean (A and G), Central West coast of Africa (B and H),
Southern Mozambique Channel (C and I), mid-latitude Southern Pacific (D and J),
North Indian Ocean (E and K), and Eastern Australian coast. (F and L).
"""

# Path to access python functions
import sys

sys.path.append("../tools/")

# Path to access intermediate data
data_path_i = "../data/ifremer_swh/"
data_path_c = "../data/ccmp2_wsp/"
data_path_decor = "../data/decor_scales/"
data_path_lsf = "../data/lsf_parameters/"

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

# call IFREMER SWH and CCMP2 WSP processed data:
swh, time_s, lat_s, lon_s = import_data("IFREMER_swh", data_path_i)
wsp, time_w, lat_w, lon_w = import_data("CCMP2_wsp", data_path_c)

# Call decorrelation time scales
###### SWH ######
nc_swh = Dataset(data_path_decor + "IFREMER_swh_decor_time_scale.nc", "r")
decor_swh = nc_swh.variables["decor_scale"][:]
time_decor_swh = num2date(nc_swh.variables["time"][:], nc_swh.variables["time"].units)
###### WSP ######
nc_wsp = Dataset(data_path_decor + "CCMP2_wsp_decor_time_scale.nc", "r")
decor_wsp = nc_wsp.variables["decor_scale"][:]
time_decor_wsp = num2date(nc_wsp.variables["time"][:], nc_wsp.variables["time"].units)

# Call ccmp2 lsf parameter data:
filename_ccmp2 = data_path_lsf + "CCMP2_wsp_lsf_parameters.nc"

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

# Calculate monthly climatologies for SWH and WSP data and decorrelation scales:
###### SWH ######
swh_clima_dict = clima_mean(date_time=np.ma.array(time_s), data=swh)
swh_mean_i = np.ma.array(swh_clima_dict["mean"])
swh_var_i = np.ma.array(swh_clima_dict["var"])
swh_n_i = np.ma.array(swh_clima_dict["N"])
decor_clima_dict_s = clima_mean(date_time=time_decor_swh, data=decor_swh)
decor_mean_swh_i = np.ma.array(decor_clima_dict_s["mean"])

###### WSP ######
wsp_clima_dict = clima_mean(date_time=np.ma.array(time_w), data=wsp)
wsp_mean_c = np.ma.array(wsp_clima_dict["mean"])
wsp_var_c = np.ma.array(wsp_clima_dict["var"])
wsp_n_c = np.ma.array(wsp_clima_dict["N"])
decor_clima_dict_w = clima_mean(date_time=time_decor_wsp, data=decor_wsp)
decor_mean_wsp_c = np.ma.array(decor_clima_dict_w["mean"])

## Compute basin-scale annual cycle fit ##
# Set partition indices for the world oceans excluding the equatorial regions and marginal seas:
ind_nph_lat, ind_nph_lon = [116, 133], [106, 266]
ind_npl_lat, ind_npl_lon = [75, 117], [106, 275]
ind_nawh_lat, ind_nawh_lon = [116, 133], [295, 361]
ind_nawl_lat, ind_nawl_lon = [77, 117], [274, 361]
ind_nae_lat, ind_nae_lon = [77, 133], [0, 24]
ind_sp_lat, ind_sp_lon = [0, 60], [131, 292]
ind_saw_lat, ind_saw_lon = [0, 68], [291, 360]
ind_sae_lat, ind_sae_lon = [0, 68], [0, 24]
ind_nioh_lat, ind_nioh_lon = [89, 95], [60, 105]
ind_niol_lat, ind_niol_lon = [77, 90], [51, 105]
ind_mioh_lat, ind_mioh_lon = [69, 78], [40, 98]
ind_miol_lat, ind_miol_lon = [59, 70], [59, 70]
ind_sio_lat, ind_sio_lon = [0, 60], [23, 132]

# Partition mean, variance, number of observations, and decorrelation scales
####### Mean #######
swh_nph = swh_mean_i[
    :, ind_nph_lat[0] : ind_nph_lat[1], ind_nph_lon[0] : ind_nph_lon[1]
]
swh_npl = swh_mean_i[
    :, ind_npl_lat[0] : ind_npl_lat[1], ind_npl_lon[0] : ind_nph_lon[1]
]
swh_nawh = swh_mean_i[
    :, ind_nawh_lat[0] : ind_nawh_lat[1], ind_nawh_lon[0] : ind_nawh_lon[1]
]
swh_nawl = swh_mean_i[
    :, ind_nawl_lat[0] : ind_nawl_lat[1], ind_nawl_lon[0] : ind_nawl_lon[1]
]
swh_nae = swh_mean_i[
    :, ind_nae_lat[0] : ind_nae_lat[1], ind_nae_lon[0] : ind_nae_lon[1]
]
swh_sp = swh_mean_i[:, ind_sp_lat[0] : ind_sp_lat[1], ind_sp_lon[0] : ind_sp_lon[1]]
swh_saw = swh_mean_i[
    :, ind_saw_lat[0] : ind_saw_lat[1], ind_saw_lon[0] : ind_saw_lon[1]
]
swh_sae = swh_mean_i[
    :, ind_sae_lat[0] : ind_sae_lat[1], ind_sae_lon[0] : ind_sae_lon[1]
]
swh_nioh = swh_mean_i[
    :, ind_nioh_lat[0] : ind_nioh_lat[1], ind_nioh_lon[0] : ind_nioh_lon[1]
]
swh_niol = swh_mean_i[
    :, ind_niol_lat[0] : ind_niol_lat[1], ind_niol_lon[0] : ind_niol_lon[1]
]
swh_mioh = swh_mean_i[
    :, ind_mioh_lat[0] : ind_mioh_lat[1], ind_mioh_lon[0] : ind_mioh_lon[1]
]
swh_miol = swh_mean_i[
    :, ind_miol_lat[0] : ind_miol_lat[1], ind_miol_lon[0] : ind_miol_lon[1]
]
swh_sio = swh_mean_i[
    :, ind_sio_lat[0] : ind_sio_lat[1], ind_sio_lon[0] : ind_sio_lon[1]
]
swh_reg = [
    swh_nph,
    swh_npl,
    swh_nawh,
    swh_nawl,
    swh_nae,
    swh_sp,
    swh_saw,
    swh_sae,
    swh_nioh,
    swh_niol,
    swh_mioh,
    swh_miol,
    swh_sio,
]

####### Variance #######
swh_sig_nph = swh_var_i[
    :, ind_nph_lat[0] : ind_nph_lat[1], ind_nph_lon[0] : ind_nph_lon[1]
]
swh_sig_npl = swh_var_i[
    :, ind_npl_lat[0] : ind_npl_lat[1], ind_npl_lon[0] : ind_nph_lon[1]
]
swh_sig_nawh = swh_var_i[
    :, ind_nawh_lat[0] : ind_nawh_lat[1], ind_nawh_lon[0] : ind_nawh_lon[1]
]
swh_sig_nawl = swh_var_i[
    :, ind_nawl_lat[0] : ind_nawl_lat[1], ind_nawl_lon[0] : ind_nawl_lon[1]
]
swh_sig_nae = swh_var_i[
    :, ind_nae_lat[0] : ind_nae_lat[1], ind_nae_lon[0] : ind_nae_lon[1]
]
swh_sig_sp = swh_var_i[:, ind_sp_lat[0] : ind_sp_lat[1], ind_sp_lon[0] : ind_sp_lon[1]]
swh_sig_saw = swh_var_i[
    :, ind_saw_lat[0] : ind_saw_lat[1], ind_saw_lon[0] : ind_saw_lon[1]
]
swh_sig_sae = swh_var_i[
    :, ind_sae_lat[0] : ind_sae_lat[1], ind_sae_lon[0] : ind_sae_lon[1]
]
swh_sig_nioh = swh_var_i[
    :, ind_nioh_lat[0] : ind_nioh_lat[1], ind_nioh_lon[0] : ind_nioh_lon[1]
]
swh_sig_niol = swh_var_i[
    :, ind_niol_lat[0] : ind_niol_lat[1], ind_niol_lon[0] : ind_niol_lon[1]
]
swh_sig_mioh = swh_var_i[
    :, ind_mioh_lat[0] : ind_mioh_lat[1], ind_mioh_lon[0] : ind_mioh_lon[1]
]
swh_sig_miol = swh_var_i[
    :, ind_miol_lat[0] : ind_miol_lat[1], ind_miol_lon[0] : ind_miol_lon[1]
]
swh_sig_sio = swh_var_i[
    :, ind_sio_lat[0] : ind_sio_lat[1], ind_sio_lon[0] : ind_sio_lon[1]
]
swh_reg_sig = [
    swh_sig_nph,
    swh_sig_npl,
    swh_sig_nawh,
    swh_sig_nawl,
    swh_sig_nae,
    swh_sig_sp,
    swh_sig_saw,
    swh_sig_sae,
    swh_sig_nioh,
    swh_sig_niol,
    swh_sig_mioh,
    swh_sig_miol,
    swh_sig_sio,
]

####### Number of data points #######
swh_n_nph = swh_n_i[:, ind_nph_lat[0] : ind_nph_lat[1], ind_nph_lon[0] : ind_nph_lon[1]]
swh_n_npl = swh_n_i[:, ind_npl_lat[0] : ind_npl_lat[1], ind_npl_lon[0] : ind_nph_lon[1]]
swh_n_nawh = swh_n_i[
    :, ind_nawh_lat[0] : ind_nawh_lat[1], ind_nawh_lon[0] : ind_nawh_lon[1]
]
swh_n_nawl = swh_n_i[
    :, ind_nawl_lat[0] : ind_nawl_lat[1], ind_nawl_lon[0] : ind_nawl_lon[1]
]
swh_n_nae = swh_n_i[:, ind_nae_lat[0] : ind_nae_lat[1], ind_nae_lon[0] : ind_nae_lon[1]]
swh_n_sp = swh_n_i[:, ind_sp_lat[0] : ind_sp_lat[1], ind_sp_lon[0] : ind_sp_lon[1]]
swh_n_saw = swh_n_i[:, ind_saw_lat[0] : ind_saw_lat[1], ind_saw_lon[0] : ind_saw_lon[1]]
swh_n_sae = swh_n_i[:, ind_sae_lat[0] : ind_sae_lat[1], ind_sae_lon[0] : ind_sae_lon[1]]
swh_n_nioh = swh_n_i[
    :, ind_nioh_lat[0] : ind_nioh_lat[1], ind_nioh_lon[0] : ind_nioh_lon[1]
]
swh_n_niol = swh_n_i[
    :, ind_niol_lat[0] : ind_niol_lat[1], ind_niol_lon[0] : ind_niol_lon[1]
]
swh_n_mioh = swh_n_i[
    :, ind_mioh_lat[0] : ind_mioh_lat[1], ind_mioh_lon[0] : ind_mioh_lon[1]
]
swh_n_miol = swh_n_i[
    :, ind_miol_lat[0] : ind_miol_lat[1], ind_miol_lon[0] : ind_miol_lon[1]
]
swh_n_sio = swh_n_i[:, ind_sio_lat[0] : ind_sio_lat[1], ind_sio_lon[0] : ind_sio_lon[1]]
swh_reg_n = [
    swh_nph,
    swh_n_npl,
    swh_n_nawh,
    swh_n_nawl,
    swh_n_nae,
    swh_n_sp,
    swh_n_saw,
    swh_n_sae,
    swh_n_nioh,
    swh_n_niol,
    swh_n_mioh,
    swh_n_miol,
    swh_n_sio,
]

####### Decorrelation time scale #######
decor_mean_swh_nph = decor_mean_swh_i[
    :, ind_nph_lat[0] : ind_nph_lat[1], ind_nph_lon[0] : ind_nph_lon[1]
]
decor_mean_swh_npl = decor_mean_swh_i[
    :, ind_npl_lat[0] : ind_npl_lat[1], ind_npl_lon[0] : ind_nph_lon[1]
]
decor_mean_swh_nawh = decor_mean_swh_i[
    :, ind_nawh_lat[0] : ind_nawh_lat[1], ind_nawh_lon[0] : ind_nawh_lon[1]
]
decor_mean_swh_nawl = decor_mean_swh_i[
    :, ind_nawl_lat[0] : ind_nawl_lat[1], ind_nawl_lon[0] : ind_nawl_lon[1]
]
decor_mean_swh_nae = decor_mean_swh_i[
    :, ind_nae_lat[0] : ind_nae_lat[1], ind_nae_lon[0] : ind_nae_lon[1]
]
decor_mean_swh_sp = decor_mean_swh_i[
    :, ind_sp_lat[0] : ind_sp_lat[1], ind_sp_lon[0] : ind_sp_lon[1]
]
decor_mean_swh_saw = decor_mean_swh_i[
    :, ind_saw_lat[0] : ind_saw_lat[1], ind_saw_lon[0] : ind_saw_lon[1]
]
decor_mean_swh_sae = decor_mean_swh_i[
    :, ind_sae_lat[0] : ind_sae_lat[1], ind_sae_lon[0] : ind_sae_lon[1]
]
decor_mean_swh_nioh = decor_mean_swh_i[
    :, ind_nioh_lat[0] : ind_nioh_lat[1], ind_nioh_lon[0] : ind_nioh_lon[1]
]
decor_mean_swh_niol = decor_mean_swh_i[
    :, ind_niol_lat[0] : ind_niol_lat[1], ind_niol_lon[0] : ind_niol_lon[1]
]
decor_mean_swh_mioh = decor_mean_swh_i[
    :, ind_mioh_lat[0] : ind_mioh_lat[1], ind_mioh_lon[0] : ind_mioh_lon[1]
]
decor_mean_swh_miol = decor_mean_swh_i[
    :, ind_miol_lat[0] : ind_miol_lat[1], ind_miol_lon[0] : ind_miol_lon[1]
]
decor_mean_swh_sio = decor_mean_swh_i[
    :, ind_sio_lat[0] : ind_sio_lat[1], ind_sio_lon[0] : ind_sio_lon[1]
]
swh_reg_dcor = [
    decor_mean_swh_nph,
    decor_mean_swh_npl,
    decor_mean_swh_nawh,
    decor_mean_swh_nawl,
    decor_mean_swh_nae,
    decor_mean_swh_sp,
    decor_mean_swh_saw,
    decor_mean_swh_sae,
    decor_mean_swh_nioh,
    decor_mean_swh_niol,
    decor_mean_swh_mioh,
    decor_mean_swh_miol,
    decor_mean_swh_sio,
]

# For each partition:
# Compute Mean
swh_reg_m = [np.ma.mean(np.ma.mean(ireg, axis=1), axis=1) for ireg in swh_reg]
swh_reg_m_n = [
    np.ma.mean(np.ma.array([swh_reg_m[0], swh_reg_m[1]]), axis=0),
    np.ma.mean(np.ma.array([swh_reg_m[2], swh_reg_m[3], swh_reg_m[4]]), axis=0),
    swh_reg_m[5],
    np.ma.mean(np.ma.array([swh_reg_m[6], swh_reg_m[7]]), axis=0),
    np.ma.mean(
        np.ma.array(
            [swh_reg_m[8], swh_reg_m[9], swh_reg_m[10], swh_reg_m[11], swh_reg_m[12]]
        ),
        axis=0,
    ),
]

# Compute average Variance
swh_reg_var_m = [np.ma.mean(np.ma.mean(ireg, axis=1), axis=1) for ireg in swh_reg_sig]
swh_reg_var_m_n = [
    np.ma.mean(np.ma.array([swh_reg_var_m[0], swh_reg_var_m[1]]), axis=0),
    np.ma.mean(
        np.ma.array([swh_reg_var_m[2], swh_reg_var_m[3], swh_reg_var_m[4]]), axis=0
    ),
    swh_reg_var_m[5],
    np.ma.mean(np.ma.array([swh_reg_var_m[6], swh_reg_var_m[7]]), axis=0),
    np.ma.mean(
        np.ma.array(
            [
                swh_reg_var_m[8],
                swh_reg_var_m[9],
                swh_reg_var_m[10],
                swh_reg_var_m[11],
                swh_reg_var_m[12],
            ]
        ),
        axis=0,
    ),
]

# compute the average of the decorrelation time scale and number of data points for each region:
swh_reg_dcor_m = [np.ma.mean(np.ma.mean(ireg, axis=1), axis=1) for ireg in swh_reg_dcor]
swh_reg_dcor_m_n = [
    np.ma.mean(np.ma.array([swh_reg_dcor_m[0], swh_reg_dcor_m[1]]), axis=0),
    np.ma.mean(
        np.ma.array([swh_reg_dcor_m[2], swh_reg_dcor_m[3], swh_reg_dcor_m[4]]), axis=0
    ),
    swh_reg_dcor_m[5],
    np.ma.mean(np.ma.array([swh_reg_dcor_m[6], swh_reg_dcor_m[7]]), axis=0),
    np.ma.mean(
        np.ma.array(
            [
                swh_reg_dcor_m[8],
                swh_reg_dcor_m[9],
                swh_reg_dcor_m[10],
                swh_reg_dcor_m[11],
                swh_reg_dcor_m[12],
            ]
        ),
        axis=0,
    ),
]

# compute the average of number of data points for each region:
swh_reg_n_m = [np.ma.mean(np.ma.mean(ireg, axis=1), axis=1) for ireg in swh_reg_n]
swh_reg_n_m_n = [
    np.ma.mean(np.ma.array([swh_reg_n_m[0], swh_reg_n_m[1]]), axis=0),
    np.ma.mean(np.ma.array([swh_reg_n_m[2], swh_reg_n_m[3], swh_reg_n_m[4]]), axis=0),
    swh_reg_n_m[5],
    np.ma.mean(np.ma.array([swh_reg_n_m[6], swh_reg_n_m[7]]), axis=0),
    np.ma.mean(
        np.ma.array(
            [
                swh_reg_n_m[8],
                swh_reg_n_m[9],
                swh_reg_n_m[10],
                swh_reg_n_m[11],
                swh_reg_n_m[12],
            ]
        ),
        axis=0,
    ),
]

# Compute degrees of freedom or N_effective for stdm calculation:
n_eff_reg = []
for ireg in range(5):
    n_eff_reg.append(swh_reg_n_m_n[ireg] / swh_reg_dcor_m_n[0])

# compute the standard error of the mean
swh_reg_stdm = []
for ireg in range(5):
    swh_reg_stdm.append(np.sqrt(swh_reg_var_m_n[ireg]) / np.sqrt(n_eff_reg[ireg]))

# Perform a weighted least squares annual cycle fit to the mean regional data
hfit_basin_swh, x_data_basin, x_data_sigma_basin = [], [], []
for ireg in range(5):

    # Perform weigthed least squares fit
    hfit, x_data, x_data_sigma = weighted_least_square_fit(
        data=np.ma.copy(swh_reg_m_n[ireg]),
        sigma=np.ma.copy(swh_reg_stdm[ireg]),
        trend="sinusoidal",
        parameters=3,
        period=12,
    )

    # save data
    hfit_basin_swh.append(hfit)
    x_data_basin.append(x_data)
    x_data_sigma_basin.append(x_data_sigma)

# Compute regional climatologies:
####### Central South Atlantic #######
# intialize variable
ngrid = 4
lon_grid = 350
lat_grid = 62
location = ["SH", "west"]

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
    lon_s,
    lat_s,
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
    lon_w,
    lat_w,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_wsp_c,
    location,
    lsf="weighted",
    parameters=5,
)

# Perform weigthed least squares fit for basin-scale climatology projected onto regional climatology
hfit_reg_sa, x_data_reg, x_data_sigma_reg = weighted_least_square_fit(
    swh_mean_sa,
    swh_stdm_sa,
    trend="sinusoidal",
    parameters=2,
    period=12,
    phase=np.arctan2(x_data_basin[3][1], x_data_basin[3][2]),
)

####### Central African Coast #######
# intialize variable
ngrid = 4
lon_grid = 6
lat_grid = 58
location = ["SH", "east"]

# compute SWH regional climatology
(
    swh_mean_ac,
    swh_stdm_ac,
    swh_hfit,
    swh_x_data,
    swh_residual,
    grid_cor_ac,
) = regional_clima(
    swh_mean_i,
    swh_var_i,
    swh_n_i,
    lon_s,
    lat_s,
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
    wsp_mean_ac,
    wsp_stdm_ac,
    wsp_hfit,
    wsp_x_data,
    wsp_residual,
    grid_cor_ac,
) = regional_clima(
    wsp_mean_c,
    wsp_var_c,
    wsp_n_c,
    lon_w,
    lat_w,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_wsp_c,
    location,
    lsf="weighted",
    parameters=5,
)

# Perform weigthed least squares fit for basin-scale climatology projected onto regional climatology
hfit_reg_ac, x_data_reg, x_data_sigma_reg = weighted_least_square_fit(
    swh_mean_ac,
    swh_stdm_ac,
    trend="sinusoidal",
    parameters=2,
    period=12,
    phase=np.arctan2(x_data_basin[3][1], x_data_basin[3][2]),
)

####### Western Madagascar Coast #######
# intialize variable
ngrid = 4
lon_grid = 37
lat_grid = 39
location = ["SH", "east"]

# compute SWH regional climatology
(
    swh_mean_mc,
    swh_stdm_mc,
    swh_hfit,
    swh_x_data,
    swh_residual,
    grid_cor_mc,
) = regional_clima(
    swh_mean_i,
    swh_var_i,
    swh_n_i,
    lon_s,
    lat_s,
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
    wsp_mean_mc,
    wsp_stdm_mc,
    wsp_hfit,
    wsp_x_data,
    wsp_residual,
    grid_cor_mc,
) = regional_clima(
    wsp_mean_c,
    wsp_var_c,
    wsp_n_c,
    lon_w,
    lat_w,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_wsp_c,
    location,
    lsf="weighted",
    parameters=5,
)

# Perform weigthed least squares fit for basin-scale climatology projected onto regional climatology
hfit_reg_mc, x_data_reg, x_data_sigma_reg = weighted_least_square_fit(
    swh_mean_mc,
    swh_stdm_mc,
    trend="sinusoidal",
    parameters=2,
    period=12,
    phase=np.arctan2(x_data_basin[4][1], x_data_basin[4][2]),
)

####### Central South Pacific #######
# intialize variable
ngrid = 4
lon_grid = 237
lat_grid = 46
location = ["SH", "west"]

# compute SWH regional climatology
(
    swh_mean_sp,
    swh_stdm_sp,
    swh_hfit,
    swh_x_data,
    swh_residual,
    grid_cor_sp,
) = regional_clima(
    swh_mean_i,
    swh_var_i,
    swh_n_i,
    lon_s,
    lat_s,
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
    wsp_mean_sp,
    wsp_stdm_sp,
    wsp_hfit,
    wsp_x_data,
    wsp_residual,
    grid_cor_sp,
) = regional_clima(
    wsp_mean_c,
    wsp_var_c,
    wsp_n_c,
    lon_w,
    lat_w,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_wsp_c,
    location,
    lsf="weighted",
    parameters=5,
)

# Perform weigthed least squares fit for basin-scale climatology projected onto regional climatology
hfit_reg_sp, x_data_reg, x_data_sigma_reg = weighted_least_square_fit(
    swh_mean_sp,
    swh_stdm_sp,
    trend="sinusoidal",
    parameters=2,
    period=12,
    phase=np.arctan2(x_data_basin[2][1], x_data_basin[2][2]),
)

####### North Indian Ocean #######
# intialize variable
ngrid = 4
lon_grid = 67
lat_grid = 65
location = ["SH", "east"]

# compute SWH regional climatology
(
    swh_mean_io,
    swh_stdm_io,
    swh_hfit,
    swh_x_data,
    swh_residual,
    grid_cor_io,
) = regional_clima(
    swh_mean_i,
    swh_var_i,
    swh_n_i,
    lon_s,
    lat_s,
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
    wsp_mean_io,
    wsp_stdm_io,
    wsp_hfit,
    wsp_x_data,
    wsp_residual,
    grid_cor_io,
) = regional_clima(
    wsp_mean_c,
    wsp_var_c,
    wsp_n_c,
    lon_w,
    lat_w,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_wsp_c,
    location,
    lsf="weighted",
    parameters=5,
)

# Perform weigthed least squares fit for basin-scale climatology projected onto regional climatology
hfit_reg_io, x_data_reg, x_data_sigma_reg = weighted_least_square_fit(
    swh_mean_io,
    swh_stdm_io,
    trend="sinusoidal",
    parameters=2,
    period=12,
    phase=np.arctan2(x_data_basin[4][1], x_data_basin[4][2]),
)

####### Eastern Australia Coast #######
# intialize variable
ngrid = 4
lon_grid = 155
lat_grid = 40
location = ["SH", "east"]

# compute SWH regional climatology
(
    swh_mean_ea,
    swh_stdm_ea,
    swh_hfit,
    swh_x_data,
    swh_residual,
    grid_cor_ea,
) = regional_clima(
    swh_mean_i,
    swh_var_i,
    swh_n_i,
    lon_s,
    lat_s,
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
    wsp_mean_ea,
    wsp_stdm_ea,
    wsp_hfit,
    wsp_x_data,
    wsp_residual,
    grid_cor_ea,
) = regional_clima(
    wsp_mean_c,
    wsp_var_c,
    wsp_n_c,
    lon_w,
    lat_w,
    lon_grid,
    lat_grid,
    ngrid,
    decor_mean_wsp_c,
    location,
    lsf="weighted",
    parameters=5,
)

# Perform weigthed least squares fit for basin-scale climatology projected onto regional climatology
hfit_reg_ea, x_data_reg, x_data_sigma_reg = weighted_least_square_fit(
    swh_mean_ea,
    swh_stdm_ea,
    trend="sinusoidal",
    parameters=2,
    period=12,
    phase=np.arctan2(x_data_basin[2][1], x_data_basin[2][2]),
)

# Set variables for plotting
projection = ccrs.PlateCarree(central_longitude=180.0)
time = np.arange(1, 13, 1)
time_ticks = [
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
ylim = [[-1, 4], [3, 10]]
resolution = "10m"

# initialize subplot axes:
fig = plt.figure(figsize=(22, 28))

############## Subplot 1  #################
# Central South Atlantic
ax1 = fig.add_subplot(621, projection=ccrs.PlateCarree(central_longitude=180.0))
cart.set_subplots(
    ax1, projection, resolution, lon_min=154, lon_max=179, lat_min=-12, lat_max=0
)
cs1 = ax1.pcolormesh(
    lon_reg_saw,
    lat_reg_saw,
    s_atlantic_w_m,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
cart.grid_lines_rc(
    ax1,
    xticks=[-26, -22, -18, -14, -10, -6, -2],
    yticks=[-12, -8, -4, 0],
    fontsize=15,
    linewidth=1,
    color="gray",
    alpha=0.7,
    linestyle="--",
    grid=True,
)
cax1 = plt.axes([0.1, 0.06, 0.4, 0.03])
cart.set_cbar(
    cs1,
    cax1,
    fig,
    orientation="horizontal",
    extend="both",
    cbar_label="$months$",
    nbins=None,
    fontsize=18,
    cbar_ticks=[
        np.arange(-np.pi, np.pi + 0.5, (np.pi + np.pi) / 6).tolist(),
        ["Jun", "Aug", "Oct", "Dec", "Feb", "Apr", "Jun"],
    ],
    task="custom ticks",
)
cart.subplot_label(
    ax1,
    xdist_label=0.1,
    ydist_label=0.88,
    subplot_label="A",
    form="box",
    fs_shade=28,
    fs_main=20,
    color="black",
)

# hind every other ticklabel
for label in ax1.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)

# Set regional climatology grid box
ax1.add_patch(
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

######## Subplot 2 ########
# Central South Atlantic
ax2 = fig.add_subplot(622)
regional_clima_plot(
    ax2,
    swh_mean_sa,
    swh_stdm_sa,
    hfit_reg_sa,
    wsp_mean_sa,
    wsp_stdm_sa,
    None,
    None,
    None,
    None,
    None,
    time,
    time_ticks,
    xlim,
    ylim,
    subplot_label="G",
    fontsize=20,
    linewidth=1.5,
    task="C_I_res",
    grid=False,
)


############## Subplot 3  #################
# Central African Coast
ax3 = fig.add_subplot(623, projection=ccrs.PlateCarree(central_longitude=180.0))
cart.set_subplots(
    ax3, projection, resolution, lon_min=-180, lon_max=-160, lat_min=-20, lat_max=0
)
cs3 = ax3.pcolormesh(
    lon_reg_sae,
    lat_reg_sae,
    s_atlantic_e_m,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
cart.grid_lines_rc(
    ax3,
    xticks=[2, 6, 10, 14, 18, 22],
    yticks=[-20, -16, -12, -8, -4, 0],
    fontsize=15,
    linewidth=1,
    color="gray",
    alpha=0.7,
    linestyle="--",
    grid=True,
)
cart.subplot_label(
    ax3,
    xdist_label=0.1,
    ydist_label=0.88,
    subplot_label="B",
    form="box",
    fs_shade=28,
    fs_main=20,
    color="black",
)

# hind every other ticklabel
for label in ax3.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)

# Set regional climatology grid box
ax3.add_patch(
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

########### Subplot 4 ###########
# Central African Coast
ax4 = fig.add_subplot(624)
regional_clima_plot(
    ax4,
    swh_mean_ac,
    swh_stdm_ac,
    hfit_reg_ac,
    wsp_mean_ac,
    wsp_stdm_ac,
    None,
    None,
    None,
    None,
    None,
    time,
    time_ticks,
    xlim,
    ylim,
    subplot_label="H",
    fontsize=20,
    linewidth=1.5,
    task="C_I_res",
    grid=False,
)

############## Subplot 5  #################
# Western Madagascar Coast
ax5 = fig.add_subplot(625, projection=ccrs.PlateCarree(central_longitude=180.0))
cart.set_subplots(
    ax5, projection, resolution, lon_min=-150, lon_max=-130, lat_min=-34, lat_max=-18
)
cs5 = ax5.pcolormesh(
    lon_reg_si,
    lat_reg_si,
    s_indian_m,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
cart.grid_lines_rc(
    ax5,
    xticks=[25, 29, 33, 37, 41, 45, 49],
    yticks=[-35, -31, -27, -23, -19, -15],
    fontsize=15,
    linewidth=1,
    color="gray",
    alpha=0.7,
    linestyle="--",
    grid=True,
)
cart.subplot_label(
    ax5,
    xdist_label=0.1,
    ydist_label=0.88,
    subplot_label="C",
    form="box",
    fs_shade=28,
    fs_main=20,
    color="black",
)

# hind every other ticklabel
for label in ax5.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)

# Set regional climatology grid box
ax5.add_patch(
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

########### Subplot 6 ###########
# Western Madagascar Coast
ax6 = fig.add_subplot(626)
regional_clima_plot(
    ax6,
    swh_mean_mc,
    swh_stdm_mc,
    hfit_reg_mc,
    wsp_mean_mc,
    wsp_stdm_mc,
    None,
    None,
    None,
    None,
    None,
    time,
    time_ticks,
    xlim,
    ylim,
    subplot_label="I",
    fontsize=20,
    linewidth=1.5,
    task="C_I_res",
    grid=False,
)

############## Subplot 7  #################
# Central South Pacific
ax7 = fig.add_subplot(627, projection=ccrs.PlateCarree(central_longitude=180.0))
cart.set_subplots(
    ax7, projection, resolution, lon_min=55, lon_max=100, lat_min=-36, lat_max=-11
)
cs7 = ax7.pcolormesh(
    lon_reg_sp,
    lat_reg_sp,
    s_pacific_m,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
cart.grid_lines_rc(
    ax7,
    xticks=[-123, -119, -115, -111, -107, -103, -99, -95, -91, -87, -83, -79, -75, -71],
    yticks=[-36, -32, -28, -24, -20, -16, -12, -8],
    fontsize=15,
    linewidth=1,
    color="gray",
    alpha=0.7,
    linestyle="--",
    grid=True,
)
cart.subplot_label(
    ax7,
    xdist_label=0.1,
    ydist_label=0.88,
    subplot_label="D",
    form="box",
    fs_shade=28,
    fs_main=20,
    color="black",
)

# hind every other ticklabel
for label in ax7.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)

# Set regional climatology grid box
ax7.add_patch(
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

########### Subplot 8 ###########
# Central South Pacific
ax8 = fig.add_subplot(628)
regional_clima_plot(
    ax8,
    swh_mean_sp,
    swh_stdm_sp,
    hfit_reg_sp,
    wsp_mean_sp,
    wsp_stdm_sp,
    None,
    None,
    None,
    None,
    None,
    time,
    time_ticks,
    xlim,
    ylim,
    subplot_label="J",
    fontsize=20,
    linewidth=1.5,
    task="C_I_res",
    grid=False,
)

############## Subplot 9  #################
# North Indian Ocean
ax9 = fig.add_subplot(629, projection=ccrs.PlateCarree(central_longitude=180.0))
cart.set_subplots(
    ax9, projection, resolution, lon_min=-106, lon_max=-85, lat_min=7, lat_max=3
)
cs9 = ax9.pcolormesh(
    lon_reg_mih,
    lat_reg_mih,
    m_indian_h_m,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
cs9 = ax9.pcolormesh(
    lon_reg_mil,
    lat_reg_mil,
    m_indian_l_m,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
cart.grid_lines_rc(
    ax9,
    xticks=[63, 67, 71, 75, 79, 83, 87, 91, 95],
    yticks=[-5, -1, 3],
    fontsize=15,
    linewidth=1,
    color="gray",
    alpha=0.7,
    linestyle="--",
    grid=True,
)
cart.subplot_label(
    ax9,
    xdist_label=0.1,
    ydist_label=0.88,
    subplot_label="E",
    form="box",
    fs_shade=28,
    fs_main=20,
    color="black",
)

# hind every other ticklabel
for label in ax9.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)

# Set regional climatology grid box
ax9.add_patch(
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

########### Subplot 10 ###########
# North Indian Ocean
ax10 = fig.add_subplot(6, 2, 10)
regional_clima_plot(
    ax10,
    swh_mean_io,
    swh_stdm_io,
    hfit_reg_io,
    wsp_mean_io,
    wsp_stdm_io,
    None,
    None,
    None,
    None,
    None,
    time,
    time_ticks,
    xlim,
    ylim,
    subplot_label="K",
    fontsize=20,
    linewidth=1.5,
    task="C_I_res",
    grid=False,
)

############## Subplot 11  #################
# South East Pacific
ax11 = fig.add_subplot(6, 2, 11, projection=ccrs.PlateCarree(central_longitude=180.0))
cart.set_subplots(
    ax11, projection, resolution, lon_min=-25, lon_max=-0, lat_min=-34, lat_max=-15
)
cs11 = ax11.pcolormesh(
    lon_reg_sp,
    lat_reg_sp,
    s_pacific_m,
    vmin=-np.pi,
    vmax=np.pi,
    cmap=cmo.phase,
    transform=ccrs.PlateCarree(central_longitude=0.0),
)
cart.grid_lines_rc(
    ax11,
    xticks=[147, 151, 155, 159, 163, 167, 171, 175, 179],
    yticks=[-34, -30, -26, -22, -18, -14],
    fontsize=15,
    linewidth=1,
    color="gray",
    alpha=0.7,
    linestyle="--",
    grid=True,
)
cart.subplot_label(
    ax11,
    xdist_label=0.1,
    ydist_label=0.22,
    subplot_label="F",
    form="box",
    fs_shade=28,
    fs_main=18,
    color="black",
)

# hind every other ticklabel
for label in ax11.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)

# Set regional climatology grid box
ax11.add_patch(
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

# Set the aspect ratio
ax11.set_aspect(aspect="equal", adjustable="box")

########### Subplot 12 ###########
# North indian ocean
ax12 = fig.add_subplot(6, 2, 12)
regional_clima_plot(
    ax12,
    swh_mean_ea,
    swh_stdm_ea,
    hfit_reg_ea,
    wsp_mean_ea,
    wsp_stdm_ea,
    None,
    None,
    None,
    None,
    None,
    time,
    time_ticks,
    xlim,
    ylim,
    subplot_label="L",
    fontsize=20,
    linewidth=1.5,
    task="C_I_res",
    grid=False,
)


# adjust spacing for the entire figure (not the subplot)
fig.subplots_adjust(wspace=0.3, hspace=0.15)

# save figure
fig.savefig(fname="../figs/figS07", bbox_inches="tight", dpi=300)
