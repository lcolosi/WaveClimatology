"""
Figure S4 Caption: 
Basin-scale climatologies with Ifremer SWH (Solid blue) and CCMP2 WSP(Solid red).
Blue shading represents the standard error of the mean. Dotted blue is the annual
cycle weighted least-squares fitted to monthly climatology for mean SWH of the 
hemisphere ocean basin. Basins include (A) North Pacific, (B) North Atlantic, (C)
South Pacific, (D) South Atlantic, and (E) Indian Ocean where marginal seas and 
the equatorial regions across the Pacific and Atlantic oceans not considered.
"""

# Path to access python functions
import sys

sys.path.append("../tools/")

# Path to access intermediate data
data_path_i = "../data/ifremer_swh/"
data_path_c = "../data/ccmp2_wsp/"
data_path_decor = "../data/decor_scales/"

# libraries
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date
import cmocean.cm as cmo

# my functions
from data_processing import import_data
from averaging_stats import clima_mean
from lsf import weighted_least_square_fit, LSF_parameters
from regional_clima_figs import regional_clima, regional_clima_plot

# call IFREMER SWH and CCMP2 WSP processed data:
swh, time_s, lat_s, lon_s = import_data("IFREMER_swh", data_path_i)
wsp, time_w, lat_w, lon_w = import_data("CCMP2_wsp", data_path_c)

# Call decorrelation time scales
nc_swh = Dataset(data_path_decor + "IFREMER_swh_decor_time_scale.nc", "r")
decor_swh = nc_swh.variables["decor_scale"][:]
time_decor_swh = num2date(nc_swh.variables["time"][:], nc_swh.variables["time"].units)
nc_wsp = Dataset(data_path_decor + "CCMP2_wsp_decor_time_scale.nc", "r")
decor_wsp = nc_wsp.variables["decor_scale"][:]
time_decor_wsp = num2date(nc_wsp.variables["time"][:], nc_wsp.variables["time"].units)

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

########### CCMP2 WSP ###########
# Calculate monthly climatologies for SWH and WSP data and decorrelation scales:
wsp_clima_dict = clima_mean(date_time=np.ma.array(time_w), data=wsp)
wsp_mean_c = np.ma.array(wsp_clima_dict["mean"])
wsp_var_c = np.ma.array(wsp_clima_dict["var"])
wsp_n_c = np.ma.array(wsp_clima_dict["N"])
decor_clima_dict_w = clima_mean(date_time=time_decor_wsp, data=decor_wsp)
decor_mean_wsp_c = np.ma.array(decor_clima_dict_w["mean"])

# Partition mean, variance, number of observations, and decorrelation scales
####### Mean #######
wsp_nph = wsp_mean_c[
    :, ind_nph_lat[0] : ind_nph_lat[1], ind_nph_lon[0] : ind_nph_lon[1]
]
wsp_npl = wsp_mean_c[
    :, ind_npl_lat[0] : ind_npl_lat[1], ind_npl_lon[0] : ind_nph_lon[1]
]
wsp_nawh = wsp_mean_c[
    :, ind_nawh_lat[0] : ind_nawh_lat[1], ind_nawh_lon[0] : ind_nawh_lon[1]
]
wsp_nawl = wsp_mean_c[
    :, ind_nawl_lat[0] : ind_nawl_lat[1], ind_nawl_lon[0] : ind_nawl_lon[1]
]
wsp_nae = wsp_mean_c[
    :, ind_nae_lat[0] : ind_nae_lat[1], ind_nae_lon[0] : ind_nae_lon[1]
]
wsp_sp = wsp_mean_c[:, ind_sp_lat[0] : ind_sp_lat[1], ind_sp_lon[0] : ind_sp_lon[1]]
wsp_saw = wsp_mean_c[
    :, ind_saw_lat[0] : ind_saw_lat[1], ind_saw_lon[0] : ind_saw_lon[1]
]
wsp_sae = wsp_mean_c[
    :, ind_sae_lat[0] : ind_sae_lat[1], ind_sae_lon[0] : ind_sae_lon[1]
]
wsp_nioh = wsp_mean_c[
    :, ind_nioh_lat[0] : ind_nioh_lat[1], ind_nioh_lon[0] : ind_nioh_lon[1]
]
wsp_niol = wsp_mean_c[
    :, ind_niol_lat[0] : ind_niol_lat[1], ind_niol_lon[0] : ind_niol_lon[1]
]
wsp_mioh = wsp_mean_c[
    :, ind_mioh_lat[0] : ind_mioh_lat[1], ind_mioh_lon[0] : ind_mioh_lon[1]
]
wsp_miol = wsp_mean_c[
    :, ind_miol_lat[0] : ind_miol_lat[1], ind_miol_lon[0] : ind_miol_lon[1]
]
wsp_sio = wsp_mean_c[
    :, ind_sio_lat[0] : ind_sio_lat[1], ind_sio_lon[0] : ind_sio_lon[1]
]
wsp_reg = [
    wsp_nph,
    wsp_npl,
    wsp_nawh,
    wsp_nawl,
    wsp_nae,
    wsp_sp,
    wsp_saw,
    wsp_sae,
    wsp_nioh,
    wsp_niol,
    wsp_mioh,
    wsp_miol,
    wsp_sio,
]

####### Variance #######
wsp_sig_nph = wsp_var_c[
    :, ind_nph_lat[0] : ind_nph_lat[1], ind_nph_lon[0] : ind_nph_lon[1]
]
wsp_sig_npl = wsp_var_c[
    :, ind_npl_lat[0] : ind_npl_lat[1], ind_npl_lon[0] : ind_nph_lon[1]
]
wsp_sig_nawh = wsp_var_c[
    :, ind_nawh_lat[0] : ind_nawh_lat[1], ind_nawh_lon[0] : ind_nawh_lon[1]
]
wsp_sig_nawl = wsp_var_c[
    :, ind_nawl_lat[0] : ind_nawl_lat[1], ind_nawl_lon[0] : ind_nawl_lon[1]
]
wsp_sig_nae = wsp_var_c[
    :, ind_nae_lat[0] : ind_nae_lat[1], ind_nae_lon[0] : ind_nae_lon[1]
]
wsp_sig_sp = wsp_var_c[:, ind_sp_lat[0] : ind_sp_lat[1], ind_sp_lon[0] : ind_sp_lon[1]]
wsp_sig_saw = wsp_var_c[
    :, ind_saw_lat[0] : ind_saw_lat[1], ind_saw_lon[0] : ind_saw_lon[1]
]
wsp_sig_sae = wsp_var_c[
    :, ind_sae_lat[0] : ind_sae_lat[1], ind_sae_lon[0] : ind_sae_lon[1]
]
wsp_sig_nioh = wsp_var_c[
    :, ind_nioh_lat[0] : ind_nioh_lat[1], ind_nioh_lon[0] : ind_nioh_lon[1]
]
wsp_sig_niol = wsp_var_c[
    :, ind_niol_lat[0] : ind_niol_lat[1], ind_niol_lon[0] : ind_niol_lon[1]
]
wsp_sig_mioh = wsp_var_c[
    :, ind_mioh_lat[0] : ind_mioh_lat[1], ind_mioh_lon[0] : ind_mioh_lon[1]
]
wsp_sig_miol = wsp_var_c[
    :, ind_miol_lat[0] : ind_miol_lat[1], ind_miol_lon[0] : ind_miol_lon[1]
]
wsp_sig_sio = wsp_var_c[
    :, ind_sio_lat[0] : ind_sio_lat[1], ind_sio_lon[0] : ind_sio_lon[1]
]
wsp_reg_sig = [
    wsp_sig_nph,
    wsp_sig_npl,
    wsp_sig_nawh,
    wsp_sig_nawl,
    wsp_sig_nae,
    wsp_sig_sp,
    wsp_sig_saw,
    wsp_sig_sae,
    wsp_sig_nioh,
    wsp_sig_niol,
    wsp_sig_mioh,
    wsp_sig_miol,
    wsp_sig_sio,
]

####### Number of data points #######
wsp_n_nph = wsp_n_c[:, ind_nph_lat[0] : ind_nph_lat[1], ind_nph_lon[0] : ind_nph_lon[1]]
wsp_n_npl = wsp_n_c[:, ind_npl_lat[0] : ind_npl_lat[1], ind_npl_lon[0] : ind_nph_lon[1]]
wsp_n_nawh = wsp_n_c[
    :, ind_nawh_lat[0] : ind_nawh_lat[1], ind_nawh_lon[0] : ind_nawh_lon[1]
]
wsp_n_nawl = wsp_n_c[
    :, ind_nawl_lat[0] : ind_nawl_lat[1], ind_nawl_lon[0] : ind_nawl_lon[1]
]
wsp_n_nae = wsp_n_c[:, ind_nae_lat[0] : ind_nae_lat[1], ind_nae_lon[0] : ind_nae_lon[1]]
wsp_n_sp = wsp_n_c[:, ind_sp_lat[0] : ind_sp_lat[1], ind_sp_lon[0] : ind_sp_lon[1]]
wsp_n_saw = wsp_n_c[:, ind_saw_lat[0] : ind_saw_lat[1], ind_saw_lon[0] : ind_saw_lon[1]]
wsp_n_sae = wsp_n_c[:, ind_sae_lat[0] : ind_sae_lat[1], ind_sae_lon[0] : ind_sae_lon[1]]
wsp_n_nioh = wsp_n_c[
    :, ind_nioh_lat[0] : ind_nioh_lat[1], ind_nioh_lon[0] : ind_nioh_lon[1]
]
wsp_n_niol = wsp_n_c[
    :, ind_niol_lat[0] : ind_niol_lat[1], ind_niol_lon[0] : ind_niol_lon[1]
]
wsp_n_mioh = wsp_n_c[
    :, ind_mioh_lat[0] : ind_mioh_lat[1], ind_mioh_lon[0] : ind_mioh_lon[1]
]
wsp_n_miol = wsp_n_c[
    :, ind_miol_lat[0] : ind_miol_lat[1], ind_miol_lon[0] : ind_miol_lon[1]
]
wsp_n_sio = wsp_n_c[:, ind_sio_lat[0] : ind_sio_lat[1], ind_sio_lon[0] : ind_sio_lon[1]]
wsp_reg_n = [
    wsp_nph,
    wsp_n_npl,
    wsp_n_nawh,
    wsp_n_nawl,
    wsp_n_nae,
    wsp_n_sp,
    wsp_n_saw,
    wsp_n_sae,
    wsp_n_nioh,
    wsp_n_niol,
    wsp_n_mioh,
    wsp_n_miol,
    wsp_n_sio,
]

####### Decorrelation time scale #######
decor_mean_wsp_nph = decor_mean_wsp_c[
    :, ind_nph_lat[0] : ind_nph_lat[1], ind_nph_lon[0] : ind_nph_lon[1]
]
decor_mean_wsp_npl = decor_mean_wsp_c[
    :, ind_npl_lat[0] : ind_npl_lat[1], ind_npl_lon[0] : ind_nph_lon[1]
]
decor_mean_wsp_nawh = decor_mean_wsp_c[
    :, ind_nawh_lat[0] : ind_nawh_lat[1], ind_nawh_lon[0] : ind_nawh_lon[1]
]
decor_mean_wsp_nawl = decor_mean_wsp_c[
    :, ind_nawl_lat[0] : ind_nawl_lat[1], ind_nawl_lon[0] : ind_nawl_lon[1]
]
decor_mean_wsp_nae = decor_mean_wsp_c[
    :, ind_nae_lat[0] : ind_nae_lat[1], ind_nae_lon[0] : ind_nae_lon[1]
]
decor_mean_wsp_sp = decor_mean_wsp_c[
    :, ind_sp_lat[0] : ind_sp_lat[1], ind_sp_lon[0] : ind_sp_lon[1]
]
decor_mean_wsp_saw = decor_mean_wsp_c[
    :, ind_saw_lat[0] : ind_saw_lat[1], ind_saw_lon[0] : ind_saw_lon[1]
]
decor_mean_wsp_sae = decor_mean_wsp_c[
    :, ind_sae_lat[0] : ind_sae_lat[1], ind_sae_lon[0] : ind_sae_lon[1]
]
decor_mean_wsp_nioh = decor_mean_wsp_c[
    :, ind_nioh_lat[0] : ind_nioh_lat[1], ind_nioh_lon[0] : ind_nioh_lon[1]
]
decor_mean_wsp_niol = decor_mean_wsp_c[
    :, ind_niol_lat[0] : ind_niol_lat[1], ind_niol_lon[0] : ind_niol_lon[1]
]
decor_mean_wsp_mioh = decor_mean_wsp_c[
    :, ind_mioh_lat[0] : ind_mioh_lat[1], ind_mioh_lon[0] : ind_mioh_lon[1]
]
decor_mean_wsp_miol = decor_mean_wsp_c[
    :, ind_miol_lat[0] : ind_miol_lat[1], ind_miol_lon[0] : ind_miol_lon[1]
]
decor_mean_wsp_sio = decor_mean_wsp_c[
    :, ind_sio_lat[0] : ind_sio_lat[1], ind_sio_lon[0] : ind_sio_lon[1]
]
wsp_reg_dcor = [
    decor_mean_wsp_nph,
    decor_mean_wsp_npl,
    decor_mean_wsp_nawh,
    decor_mean_wsp_nawl,
    decor_mean_wsp_nae,
    decor_mean_wsp_sp,
    decor_mean_wsp_saw,
    decor_mean_wsp_sae,
    decor_mean_wsp_nioh,
    decor_mean_wsp_niol,
    decor_mean_wsp_mioh,
    decor_mean_wsp_miol,
    decor_mean_wsp_sio,
]

# For each partition:
# Compute Mean
wsp_reg_m = [np.ma.mean(np.ma.mean(ireg, axis=1), axis=1) for ireg in wsp_reg]
wsp_reg_m_n = [
    np.ma.mean(np.ma.array([wsp_reg_m[0], wsp_reg_m[1]]), axis=0),
    np.ma.mean(np.ma.array([wsp_reg_m[2], wsp_reg_m[3], wsp_reg_m[4]]), axis=0),
    wsp_reg_m[5],
    np.ma.mean(np.ma.array([wsp_reg_m[6], wsp_reg_m[7]]), axis=0),
    np.ma.mean(
        np.ma.array(
            [wsp_reg_m[8], wsp_reg_m[9], wsp_reg_m[10], wsp_reg_m[11], wsp_reg_m[12]]
        ),
        axis=0,
    ),
]

# Compute average Variance
wsp_reg_var_m = [np.ma.mean(np.ma.mean(ireg, axis=1), axis=1) for ireg in wsp_reg_sig]
wsp_reg_var_m_n = [
    np.ma.mean(np.ma.array([wsp_reg_var_m[0], wsp_reg_var_m[1]]), axis=0),
    np.ma.mean(
        np.ma.array([wsp_reg_var_m[2], wsp_reg_var_m[3], wsp_reg_var_m[4]]), axis=0
    ),
    wsp_reg_var_m[5],
    np.ma.mean(np.ma.array([wsp_reg_var_m[6], wsp_reg_var_m[7]]), axis=0),
    np.ma.mean(
        np.ma.array(
            [
                wsp_reg_var_m[8],
                wsp_reg_var_m[9],
                wsp_reg_var_m[10],
                wsp_reg_var_m[11],
                wsp_reg_var_m[12],
            ]
        ),
        axis=0,
    ),
]

# compute the average of the decorrelation time scale and number of data points for each region:
wsp_reg_dcor_m = [np.ma.mean(np.ma.mean(ireg, axis=1), axis=1) for ireg in wsp_reg_dcor]
wsp_reg_dcor_m_n = [
    np.ma.mean(np.ma.array([wsp_reg_dcor_m[0], wsp_reg_dcor_m[1]]), axis=0),
    np.ma.mean(
        np.ma.array([wsp_reg_dcor_m[2], wsp_reg_dcor_m[3], wsp_reg_dcor_m[4]]), axis=0
    ),
    wsp_reg_dcor_m[5],
    np.ma.mean(np.ma.array([wsp_reg_dcor_m[6], wsp_reg_dcor_m[7]]), axis=0),
    np.ma.mean(
        np.ma.array(
            [
                wsp_reg_dcor_m[8],
                wsp_reg_dcor_m[9],
                wsp_reg_dcor_m[10],
                wsp_reg_dcor_m[11],
                wsp_reg_dcor_m[12],
            ]
        ),
        axis=0,
    ),
]

# compute the average of number of data points for each region:
wsp_reg_n_m = [np.ma.mean(np.ma.mean(ireg, axis=1), axis=1) for ireg in wsp_reg_n]
wsp_reg_n_m_n = [
    np.ma.mean(np.ma.array([wsp_reg_n_m[0], wsp_reg_n_m[1]]), axis=0),
    np.ma.mean(np.ma.array([wsp_reg_n_m[2], wsp_reg_n_m[3], wsp_reg_n_m[4]]), axis=0),
    wsp_reg_n_m[5],
    np.ma.mean(np.ma.array([wsp_reg_n_m[6], wsp_reg_n_m[7]]), axis=0),
    np.ma.mean(
        np.ma.array(
            [
                wsp_reg_n_m[8],
                wsp_reg_n_m[9],
                wsp_reg_n_m[10],
                wsp_reg_n_m[11],
                wsp_reg_n_m[12],
            ]
        ),
        axis=0,
    ),
]

# Compute degrees of freedom or N_effective for stdm calculation:
n_eff_reg = []
for ireg in range(5):
    n_eff_reg.append(wsp_reg_n_m_n[ireg] / wsp_reg_dcor_m_n[0])

# compute the standard error of the mean
wsp_reg_stdm = []
for ireg in range(5):
    wsp_reg_stdm.append(np.sqrt(wsp_reg_var_m_n[ireg]) / np.sqrt(n_eff_reg[ireg]))

# Perform a weighted least squares annual cycle fit to the mean regional data
hfit_reg_wsp = []
for ireg in range(5):

    # Perform weigthed least squares fit
    hfit, x_data, x_data_sigma = weighted_least_square_fit(
        data=np.ma.copy(wsp_reg_m_n[ireg]),
        sigma=np.ma.copy(wsp_reg_stdm[ireg]),
        trend="sinusoidal",
        parameters=3,
        period=12,
    )

    # save data
    hfit_reg_wsp.append(hfit)

# Change order of climatologies variables and fits for SH partitions:
####### Climatology #######
wsp_reg_m_n[2] = np.reshape(
    np.ma.array([wsp_reg_m_n[2][6:13], wsp_reg_m_n[2][0:6]]), (1, 12)
)[0]
wsp_reg_m_n[3] = np.reshape(
    np.ma.array([wsp_reg_m_n[3][6:13], wsp_reg_m_n[3][0:6]]), (1, 12)
)[0]
wsp_reg_m_n[4] = np.reshape(
    np.ma.array([wsp_reg_m_n[4][6:13], wsp_reg_m_n[4][0:6]]), (1, 12)
)[0]
####### stdm #######
wsp_reg_stdm[2] = np.reshape(
    np.ma.array([wsp_reg_stdm[2][6:13], wsp_reg_stdm[2][0:6]]), (1, 12)
)[0]
wsp_reg_stdm[3] = np.reshape(
    np.ma.array([wsp_reg_stdm[3][6:13], wsp_reg_stdm[3][0:6]]), (1, 12)
)[0]
wsp_reg_stdm[4] = np.reshape(
    np.ma.array([wsp_reg_stdm[4][6:13], wsp_reg_stdm[4][0:6]]), (1, 12)
)[0]
####### fits #######
hfit_reg_wsp[2] = np.reshape(
    np.ma.array([hfit_reg_wsp[2][6:13], hfit_reg_wsp[2][0:6]]), (1, 12)
)[0]
hfit_reg_wsp[3] = np.reshape(
    np.ma.array([hfit_reg_wsp[3][6:13], hfit_reg_wsp[3][0:6]]), (1, 12)
)[0]
hfit_reg_wsp[4] = np.reshape(
    np.ma.array([hfit_reg_wsp[4][6:13], hfit_reg_wsp[4][0:6]]), (1, 12)
)[0]

########### IFREMER SWH ###########
# Calculate monthly climatologies for SWH and WSP data and decorrelation scales:
swh_clima_dict = clima_mean(date_time=np.ma.array(time_s), data=swh)
swh_mean_i = np.ma.array(swh_clima_dict["mean"])
swh_var_i = np.ma.array(swh_clima_dict["var"])
swh_n_i = np.ma.array(swh_clima_dict["N"])
decor_clima_dict_s = clima_mean(date_time=time_decor_swh, data=decor_swh)
decor_mean_swh_i = np.ma.array(decor_clima_dict_s["mean"])

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
hfit_reg_swh = []
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
    hfit_reg_swh.append(hfit)

# Change order of climatologies variables and fits for SH partitions:
####### Climatology #######
swh_reg_m_n[2] = np.reshape(
    np.ma.array([swh_reg_m_n[2][6:13], swh_reg_m_n[2][0:6]]), (1, 12)
)[0]
swh_reg_m_n[3] = np.reshape(
    np.ma.array([swh_reg_m_n[3][6:13], swh_reg_m_n[3][0:6]]), (1, 12)
)[0]
swh_reg_m_n[4] = np.reshape(
    np.ma.array([swh_reg_m_n[4][6:13], swh_reg_m_n[4][0:6]]), (1, 12)
)[0]
####### stdm #######
swh_reg_stdm[2] = np.reshape(
    np.ma.array([swh_reg_stdm[2][6:13], swh_reg_stdm[2][0:6]]), (1, 12)
)[0]
swh_reg_stdm[3] = np.reshape(
    np.ma.array([swh_reg_stdm[3][6:13], swh_reg_stdm[3][0:6]]), (1, 12)
)[0]
swh_reg_stdm[4] = np.reshape(
    np.ma.array([swh_reg_stdm[4][6:13], swh_reg_stdm[4][0:6]]), (1, 12)
)[0]
####### fits #######
hfit_reg_swh[2] = np.reshape(
    np.ma.array([hfit_reg_swh[2][6:13], hfit_reg_swh[2][0:6]]), (1, 12)
)[0]
hfit_reg_swh[3] = np.reshape(
    np.ma.array([hfit_reg_swh[3][6:13], hfit_reg_swh[3][0:6]]), (1, 12)
)[0]
hfit_reg_swh[4] = np.reshape(
    np.ma.array([hfit_reg_swh[4][6:13], hfit_reg_swh[4][0:6]]), (1, 12)
)[0]

# Set variables for plotting
time = np.arange(1, 13, 1)
time_ticks_nh = [
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
    "Jul",
    "Aug",
    "Sep",
    "Oct",
    "Nov",
    "Dec",
]
time_ticks_sh = [
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
ylim = [[1, 4], [4, 10]]

# Create figure and axes
fig, axes = plt.subplots(3, 2, figsize=(16, 20))
ax1, ax2, ax3, ax4, ax5, ax6 = axes.flatten()

######## Subplot 1 ########
# North Pacific
regional_clima_plot(
    ax1,
    swh_reg_m_n[0],
    swh_reg_stdm[0],
    hfit_reg_swh[0],
    wsp_reg_m_n[0],
    wsp_reg_stdm[0],
    hfit_reg_wsp[0],
    None,
    None,
    None,
    None,
    time,
    time_ticks_nh,
    xlim,
    ylim,
    subplot_label="A",
    fontsize=20,
    linewidth=1.5,
    task="IC",
    grid=False,
)

######## Subplot 2 ########
# North Atlantic
regional_clima_plot(
    ax2,
    swh_reg_m_n[1],
    swh_reg_stdm[1],
    hfit_reg_swh[1],
    wsp_reg_m_n[1],
    wsp_reg_stdm[1],
    hfit_reg_wsp[1],
    None,
    None,
    None,
    None,
    time,
    time_ticks_nh,
    xlim,
    ylim,
    subplot_label="B",
    fontsize=20,
    linewidth=1.5,
    task="IC",
    grid=False,
)

########### Subplot 3 ###########
# South pacific
regional_clima_plot(
    ax3,
    swh_reg_m_n[2],
    swh_reg_stdm[2],
    hfit_reg_swh[2],
    wsp_reg_m_n[2],
    wsp_reg_stdm[2],
    hfit_reg_wsp[2],
    None,
    None,
    None,
    None,
    time,
    time_ticks_sh,
    xlim,
    ylim,
    subplot_label="C",
    fontsize=20,
    linewidth=1.5,
    task="IC",
    grid=False,
)

########### Subplot 4 ###########
# South Atlantic
regional_clima_plot(
    ax4,
    swh_reg_m_n[3],
    swh_reg_stdm[3],
    hfit_reg_swh[3],
    wsp_reg_m_n[3],
    wsp_reg_stdm[3],
    hfit_reg_wsp[3],
    None,
    None,
    None,
    None,
    time,
    time_ticks_sh,
    xlim,
    ylim,
    subplot_label="D",
    fontsize=20,
    linewidth=1.5,
    task="IC",
    grid=False,
)

########### Subplot 5 ###########
# Indian ocean
regional_clima_plot(
    ax5,
    swh_reg_m_n[4],
    swh_reg_stdm[4],
    hfit_reg_swh[4],
    wsp_reg_m_n[4],
    wsp_reg_stdm[4],
    hfit_reg_wsp[4],
    None,
    None,
    None,
    None,
    time,
    time_ticks_sh,
    xlim,
    ylim,
    subplot_label="E",
    fontsize=20,
    linewidth=1.5,
    task="IC",
    grid=False,
)

########### Subplot 6 ###########
# suppress output for axis 6
ax6.axis("off")

# adjust layout of subplots:
fig.tight_layout()

# adjust spacing for the entire figure
fig.subplots_adjust(wspace=0.35, hspace=0.2)

# save figure:
fig.savefig(fname="../figs/figS04", bbox_inches="tight", dpi=300)
