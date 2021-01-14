# Regional Climatology Functions
## Luke Colosi | lcolosi@ucsd.edu | December 22nd, 2020

########### Regional Climatology Function ###########
def regional_clima(
    data_mean,
    data_var,
    data_n,
    lon,
    lat,
    lon_grid,
    lat_grid,
    ngrid,
    dcor,
    loc,
    lsf,
    parameters,
):

    """
    regional_clima(data_mean, data_var, data_n, lon, lat, lon_grid, lat_grid, ngrid, dcor, loc, lsf, parameters)

        Function to compute the regional climatology from monthly climatological data with specified grid and location
        in the ocean.

        Parameters
        ----------
        data_mean : Monthly Climatology data with the following dimensions in 2D masked geospatial arrays:
                    (ntime, nlat, nlon) = (12, 133, 360)
        date_var : Monthly Climatology data's variance with the following dimensions in 2D masked
                   geospatial arrays: (ntime, nlat, nlon) = (12, 133, 360)
        data_n : Number of observations used in computing Monthly Climatology mean and variance with the
                 following dimensions: (ntime, nlat, nlon) = (12, 133, 360). Used to compute the standard
                 error of the mean for error bars.
        lon : Longitude vector
            ex: lon = np.arange(0, 360, 1)
        lat : Latitude vector
            ex: lon = np.arange(-66, 66, 1)
        lon_grid : Initial longitude grid point to compute the regional climatology
            ex: lon_grid = 230
        lat_grid : Initial latitude grid point to compute the regional climatology
            ex: lon_grid = 230
        ngrid : Specifies the size of the n by n grid box that the regional climatology will computed in.
            ex: n_grid = 2
        dcor : Decorrelation scale (climatological resolution) used for computing degrees of freedom or dof
               (used to compute standard error of the mean) where dof is defined as:
               dof = n_eff = nobs/dcor where nobs = number of observations in time series. dcor must be an array.
        loc : Specifies if the regional climatology is in the northern or southern hemisphere or if it is east or
              west of the prime meridian.
            ex: loc = [loc_lat, loc_lon] loc[0] = 'NH' or loc[0] = 'SH' and loc[1] = 'west' or loc[1] = 'east'
        lsf : Specifies whether the least squares fit model is weighted or unweighted. Options
              include: lsf = 'weighted' or 'unweighted'
        parameters : Specifies the amount of paramaters for the model. Look at weighted_least_squares_fit.py
                     documentation for details.

        Returns
        -------
        data_reg_mean : Regional climatology mean
            ex: data_reg_mean.shape = (1,12)
        data_reg_stdm : Regional climatology standard error of the mean
            ex: data_reg_stdm.shape = (1,12)
        hfit : Least square fit model.
        x_data : Least squares fit model coefficients.
        residual : Difference between the model and data.
        grid_coordinates : Indices from longitude and latitude which the climatology is computed for.

        Libraries necessary to run function
        -----------------------------------
        import numpy as np
        from unweighted_least_square_fit import least_square_fit
        from weighted_least_square_fit import weighted_least_square_fit
    """

    # Set path to my functions:
    import sys

    sys.path.append("../tools/")

    # import libraries
    import numpy as np

    # import my functions
    from lsf import least_square_fit, weighted_least_square_fit

    # case 1: west of prime meridian
    if loc[1] == "west":
        loc_lon = 360
    # case 2: east of prime meridian
    elif loc[1] == "east":
        loc_lon = 0

    # latitude and longitude grid points that will be averaged over:
    lat_grid_i = lat[lat_grid]
    lat_grid_f = lat[lat_grid + ngrid - 1]
    lon_grid_i = lon[lon_grid] - loc_lon
    lon_grid_f = lon[lon_grid + ngrid - 1] - loc_lon
    grid_cor = [lat_grid_i, lat_grid_f, lon_grid_i, lon_grid_f]

    # call mean, variance, number of observations, and decorrelation scale from grid box indices:
    data_grid_mean = data_mean[
        :, lat_grid : (lat_grid + ngrid), lon_grid : (lon_grid + ngrid)
    ]
    data_grid_var = data_var[
        :, lat_grid : (lat_grid + ngrid), lon_grid : (lon_grid + ngrid)
    ]
    data_grid_n = data_n[
        :, lat_grid : (lat_grid + ngrid), lon_grid : (lon_grid + ngrid)
    ]
    data_grid_dcor = dcor[
        :, lat_grid : (lat_grid + ngrid), lon_grid : (lon_grid + ngrid)
    ]

    # Compute the mean and average variance, decorrelation scale number of observations for the region
    data_reg_mean = np.ma.mean(
        np.ma.mean(data_grid_mean, axis=1, dtype=np.float64), axis=1
    )
    data_reg_var = np.ma.mean(
        np.ma.mean(data_grid_var, axis=1, dtype=np.float64), axis=1
    )
    dcor_reg = np.ma.mean(np.ma.mean(data_grid_dcor, axis=1, dtype=np.float64), axis=1)
    data_reg_n_mean = np.ma.mean(
        np.ma.mean(data_grid_n, axis=1, dtype=np.float64), axis=1
    )

    # Compute N_eff:
    n_eff = data_grid_n / data_grid_dcor

    # Compute the average number degrees of freedom (n_eff):
    n_eff_mean = np.ma.mean(np.ma.mean(n_eff, axis=1, dtype=np.float64), axis=1)

    # compute the standard error of the mean
    data_reg_stdm = np.sqrt(data_reg_var) / np.sqrt(n_eff_mean)

    # compute the least square fit:
    if lsf == "weighted":
        hfit, x_data, x_data_sigma = weighted_least_square_fit(
            data=np.ma.copy(data_reg_mean),
            sigma=np.ma.copy(data_reg_stdm),
            trend="sinusoidal",
            parameters=parameters,
            period=12,
        )
    elif lsf == "unweighted":
        hfit, x_data = least_square_fit(
            data=np.ma.copy(data_reg_mean),
            trend="sinusoidal",
            parameters=parameters,
            period=12,
        )

    # compute the residue between the model and regional climatology:
    residual = hfit - data_reg_mean

    # For SH regional climatologies, shift the time series such that austral summer months are center in the figure
    if loc[0] == "SH":

        # Shift
        data_reg_mean = np.reshape(
            np.ma.array([data_reg_mean[6:13], data_reg_mean[0:6]]), (1, 12)
        )[0]
        data_reg_stdm = np.reshape(
            np.ma.array([data_reg_stdm[6:13], data_reg_stdm[0:6]]), (1, 12)
        )[0]
        hfit = np.reshape(np.ma.array([hfit[6:13], hfit[0:6]]), (1, 12))[0]
        residual = np.reshape(np.ma.array([residual[6:13], residual[0:6]]), (1, 12))[0]

    return data_reg_mean, data_reg_stdm, hfit, x_data, residual, grid_cor


########### Regional Climatology Plotting Function ###########
def regional_clima_plot(
    ax,
    swh_mean,
    swh_stdm,
    swh_hfit,
    wsp_mean,
    wsp_stdm,
    wsp_hfit,
    swh_model_mean,
    swh_model_stdm,
    wsp_model_mean,
    wsp_model_stdm,
    time,
    time_ticks,
    xlim,
    ylim,
    subplot_label,
    fontsize,
    linewidth,
    task,
    grid=True,
):

    """
    region_clima_plot(ax, swh_mean, swh_stdm, swh_hfit, wsp_mean, wsp_stdm, wsp_hfit, swh_model_mean, swh_model_stdm, wsp_model_mean, wsp_model_stdm, time, time_ticks, xlim, ylim, subplot_label, fontsize, linewidth, task, grid=True)

        Function for plotting regional climatologies on a subplot.

        Parameters
        ----------
        ax : Geospatial axes for the subplot
            ex: fig, axes = plt.subplots(3, 2, figsize=(16,12),subplot_kw={'projection': projection})
                ax1, ax2, ax3, ax4, ax5, ax6 = axes.flatten()
        swh_mean : Ifremer SWH regional climatology mean
        swh_stdm : Ifremer SWH regional climatology standard error of the mean
        swh_hfit : Ifremer SWH regional climatology lsf model
        wsp_mean : CCMP v2 regional climatology mean
        wsp_stdm : CCMP v2 regional climatology standard error of the mean
        wsp_hfit : CCMP v2 regional climatology lsf model
        swh_model_mean : WW3 swh regional climatology mean
        swh_model_stdm : WW3 swh regional climatology standard error of the mean
        wsp_model_mean : WW3 wsp regional climatology mean
        wsp_model_stdm : WW3 wsp regional climatology standard error of the mean
        time : time vector for plotting regional climatologies
            ex: time = np.arange(0,13)
        time_ticks : tick labels for subplot
            ex: time_ticks = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        xlim : x_axis limits in list form
            ex: xlim = [0,13]
        ylim : y axis lims for both axes in the form of a list: ylim = [[swh_limits], [wsp_limits]]
            ex: ylim = [[1,4], [6, 10]]
        subplot_label : subplot label for paper
            ex: subplot_label = 'A'
        fontsize : Specifies the font size of the x and y tickmarks as well as the x and y axes labels.
            ex: fontsize = 15
        linewidth : Specifies the line width for the line plot.
            ex: linewidth = 2
        task : Specifies whether the plot will be displaying WW3, Ifremer, and CCMP2 comparision (ww3), Ifremer and CCMP2 climatology only with models (IC), residual of SWH with CCMP2 climatology (residual), or Ifremer and CCMP2 comparison with swh residual (C_I_res). Options include: task = 'ww3' or task = 'IC', task = 'residual', task = 'C_I_res'.
        grid : Specifies if the subplot has a grid or not. By default, grid == True meaning a grid is plotted.
               grid == False removes the grid on the figure.

        Returns
        -------
        Plots regional climatology subplot from specified region

        Libraries necessary to run function
        -----------------------------------
        import numpy as np
        import matplotlib.pyplot as plt
        import cartopy_fig_module as cart

    """

    # Path to access python functions
    import sys

    sys.path.append("../tools/")

    # import libraries:
    import numpy as np
    import matplotlib.pyplot as plt

    # Import my function
    import cartopy_figs as cart

    ########## SWH ##########
    # Set axis:
    ax_twin = ax.twinx()

    # set axis color:
    color = "tab:blue"
    ax_twin.tick_params(axis="y", labelcolor=color)

    # WW3 comparison
    if task == "ww3":

        # Plot WW3
        ax_twin.plot(
            time, swh_model_mean, ".--b", color="tab:blue", linewidth=linewidth
        )
        ax_twin.plot(time, swh_model_mean, ".b", color="tab:blue", markersize=2)
        ax_twin.fill_between(
            time,
            swh_model_mean - swh_model_stdm,
            swh_model_mean + swh_model_stdm,
            color="tab:blue",
            alpha=0.3,
        )

        # Plot Ifremer SWH
        ax_twin.plot(time, swh_mean, ".-b", linewidth=linewidth)
        ax_twin.plot(time, swh_mean, ".b", markersize=2)
        ax_twin.fill_between(
            time, swh_mean - swh_stdm, swh_mean + swh_stdm, color="blue", alpha=0.3
        )

    # Ifremer-CCMP2 comparison
    elif task == "IC":

        # Plot Ifremer SWH and lsf model
        ax_twin.plot(time, swh_mean, ".-b", linewidth=linewidth)
        ax_twin.plot(time, swh_mean, ".b", markersize=2)
        ax_twin.plot(time, swh_hfit, ".--b", linewidth=linewidth)
        ax_twin.plot(time, swh_hfit, ".b", markersize=2)
        ax_twin.fill_between(
            time, swh_mean - swh_stdm, swh_mean + swh_stdm, color=color, alpha=0.3
        )

    # Ifremer-CCMP2 residual comparison
    elif task == "residual":

        # Compute residual
        res_swh = swh_hfit - swh_mean

        # Plot Ifremer SWH residual
        ax_twin.plot(time, res_swh, ".-b", linewidth=linewidth)
        ax_twin.plot(time, res_swh, ".b", markersize=2)
        ax_twin.fill_between(time, 0, res_swh, facecolor="blue", alpha=0.2)

    # Ifremer-CCMP2 comparison and residual
    elif task == "C_I_res":

        # Plot Ifremer SWH and lsf model
        ax_twin.plot(time, swh_mean, ".-b", linewidth=linewidth)
        ax_twin.plot(time, swh_mean, ".b", markersize=2)
        ax_twin.plot(time, swh_hfit, ".--b", linewidth=linewidth)
        ax_twin.plot(time, swh_hfit, ".b", markersize=2)
        ax_twin.fill_between(
            time, swh_mean - swh_stdm, swh_mean + swh_stdm, color=color, alpha=0.3
        )

        # Compute residual
        res_swh = swh_mean - swh_hfit

        # Plot Ifremer SWH residual
        ax_twin.plot(time, res_swh, ".-k", linewidth=linewidth)
        ax_twin.plot(time, res_swh, ".k", markersize=2)
        ax_twin.fill_between(time, 0, res_swh, facecolor="black", alpha=0.2)

    # set SWH y-axis attributes:
    if task == "residual":
        ax_twin.set_ylabel("$Residual\;SWH\;m$", color=color, fontsize=fontsize)
        ax_twin.set_ylim(ylim[0])
    elif task == "C_I_res":
        ax_twin.set_ylabel("$m$", color=color, fontsize=fontsize)
        ax_twin.set_yticks(np.arange(ylim[0][0], ylim[0][1] + 0.5, 0.5))
        for label in ax_twin.yaxis.get_ticklabels()[1::2]:
            label.set_visible(False)
    else:
        ax_twin.set_ylabel("$SWH\;(m)$", color=color, fontsize=fontsize)
        ax_twin.set_ylim(ylim[0])

    ax_twin.tick_params(axis="y", labelsize=fontsize)

    # Set grid lines
    if grid == True:
        ax_twin.grid(color=color, linestyle="-.", linewidth=1, alpha=0.35)

    ########## WSP ##########
    # set color:
    color = "tab:red"
    ax.tick_params(axis="y", labelcolor=color)

    # WW3 comparison
    if task == "ww3":

        # plot WW3 WSP
        ax.plot(time, wsp_model_mean, ".--", color="tab:red", linewidth=linewidth)
        ax.plot(time, wsp_model_mean, ".", color="tab:red", markersize=2)
        ax.fill_between(
            time,
            wsp_model_mean - wsp_model_stdm,
            wsp_model_mean + wsp_model_stdm,
            color="tab:red",
            alpha=0.3,
        )

        # plot CCMP2 WSP
        ax.plot(time, wsp_mean, ".-r", linewidth=linewidth)
        ax.plot(time, wsp_mean, ".r", markersize=2)
        ax.fill_between(
            time, wsp_mean - wsp_stdm, wsp_mean + wsp_stdm, color="red", alpha=0.3
        )

    # Ifremer-CCMP2 comparison and residual
    elif task == "C_I_res":

        # plot CCMP2 WSP
        ax.plot(time, wsp_mean, ".-r", linewidth=linewidth)
        ax.plot(time, wsp_mean, ".r", markersize=2)
        ax.fill_between(
            time, wsp_mean - wsp_stdm, wsp_mean + wsp_stdm, color=color, alpha=0.3
        )

    # Ifremer-CCMP2 comparison
    else:

        # plot CCMP2 WSP
        ax.plot(time, wsp_mean, ".-r", linewidth=linewidth)
        ax.plot(time, wsp_mean, ".r", markersize=2)
        ax.plot(time, wsp_hfit, ".--r", linewidth=linewidth)
        ax.plot(time, wsp_hfit, ".r", markersize=2)
        ax.fill_between(
            time, wsp_mean - wsp_stdm, wsp_mean + wsp_stdm, color=color, alpha=0.3
        )

    # set WSP y-axis attributes:
    if task == "C_I_res":
        ax.set_ylabel("$m\,s^{-1}$", color=color, fontsize=fontsize)
        ax.set_yticks(np.arange(ylim[1][0], ylim[1][1] + 0.5, 0.5))
        for label in ax.yaxis.get_ticklabels()[1::2]:
            label.set_visible(False)
    else:
        ax.set_ylabel("$WSP\;(m\,s^{-1})$", color=color, fontsize=fontsize)
        ax.set_ylim(ylim[1])

    ax.tick_params(axis="y", labelsize=fontsize)

    # set grid lines
    if grid == True:
        ax.grid(color=color, linestyle="-.", linewidth=1, alpha=0.35)

    # set attributes of rest of subplot:
    ax.set_xlim([0, 13])
    ax.set_xticklabels(time_ticks, fontsize=fontsize)
    ax.tick_params(axis="x", labelsize=fontsize)

    # hind every other ticklabel
    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(np.arange(start + 1, end + 1, 1))
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)

    # set labels on figure
    cart.subplot_label(
        ax,
        xdist_label=0.08,
        ydist_label=0.92,
        subplot_label=subplot_label,
        form="box",
        fs_shade=28,
        fs_main=20,
        color="black",
    )
