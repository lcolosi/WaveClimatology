# Cartopy functions
## Luke Colosi | lcolosi@ucsd.edu | November 2nd, 2020

########### Initializing Subplot Function ###########
def set_subplots(ax, projection, resolution, lon_min, lon_max, lat_min, lat_max):

    """
    set_subplots(ax, projection)

        Function for placing x and y axes labels for longitude and latitude respectively

        Parameters
        ----------
        ax : geospatial axes for the subplot (cartopy object)
            ex: fig, axes = plt.subplots(3, 2, figsize=(16,12), subplot_kw={'projection': projection})
                ax1, ax2, ax3, ax4, ax5, ax6 = axes.flatten()
                ax = ax1 or ax2 or ...
        lon_min : minimum extent for longitude on the scale from -180 to 179
            ex: lon_min = -180
        lon_max : maximum extent for longitude on the scale from -180 to 179
            ex: lon_max = 179
        lat_min : minimum extent for latitude on the scale from -90 to 89
            ex: lat_min = -66
        lat_max : maximum extent for latitude on the scale from -90 to 89
            ex: lat_max = 66

        Returns
        -------
        No objects returned. A geospatial map with desired longitude and latitude extent with coastlines and land.

        Libraries necessary to run function
        -----------------------------------
        import cartopy.feature as cfeature

    """
    # import libraries
    import cartopy.feature as cfeature

    # Set extents of map
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], projection)

    # Plot coastlines and land
    ax.coastlines(resolution=resolution)
    ax.add_feature(
        cfeature.NaturalEarthFeature("physical", "land", resolution, facecolor="Gray")
    )

    return


########### Setting Axes Tickmarks Function ###########
def set_grid_ticks(
    ax, projection, xticks, yticks, xlabels, ylabels, grid, fontsize, color
):

    """
    set_grid_ticks(ax, projection, xticks, yticks, xlabels, ylabels, grid, fontsize, color)

        Function for plotting geospatial data with tick marks pointing out of the figure and tick labels at degrees.

        Paramters
        ---------
        ax: Geospatial axes for the subplot (cartopy object)
            ex: fig, axes = plt.subplots(3, 2, figsize=(16,12), subplot_kw={'projection': projection})
                ax1, ax2, ax3, ax4, ax5, ax6 = axes.flatten()
                ax = ax1 or ax2 or ...
        projection : Projection map (Rectangular only). For the gridlines, avoid adding arguments to the project.
            ex : projection = ccrs.PlateCarree()
        xticks : List of longitudinal tick marks
            ex: xticks = [0, 60, 120, 180, -120, -60]
        yticks : List of latitudinal tick marks
            ex: yticks = [-60, -30, 0, 30, 60]
        xlabels : Specify if you want x axis labels left axis. True means longitude labels
                are present.
            ex: xlabels = True
        ylabels : Specify if you want y axis labels on the bottom axis. True means latitude labels
                are present.
            ex: ylabels = True
        grid : Specify if you want grid lines. True means grid is present.
        fontsize : Specifies the fontsize of the x and y tickmarks.
        color : Specifies color of grid lines and tickmarks.

        Returns
        -------
        No objects returned. Properly formatted geospatial tick marks.

        Libraries necessary to run function
        -----------------------------------
        import cartopy.crs as ccrs
        from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

    """
    import cartopy.crs as ccrs
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

    # Set x and y tick marks
    ax.set_xticks(xticks, crs=ccrs.PlateCarree())
    ax.set_yticks(yticks, crs=ccrs.PlateCarree())

    # Set format of x and y tick labels
    lon_formatter = LongitudeFormatter()
    lat_formatter = LatitudeFormatter()

    # Apply formatting
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

    # Set size of tick marks
    ax.xaxis.set_tick_params(labelsize=fontsize, colors=color)
    ax.yaxis.set_tick_params(labelsize=fontsize, colors=color)

    # Hind tickmarks
    if not xlabels:
        for label in ax.xaxis.get_ticklabels():
            label.set_visible(False)
    if not ylabels:
        for label in ax.yaxis.get_ticklabels():
            label.set_visible(False)

    # Set grid lines
    if grid == True:
        ax.grid(linewidth=2, color=color, alpha=0.3, linestyle="--")

    return


########### Grid lines for Regional Climatology Function ###########
def grid_lines_rc(
    ax, xticks, yticks, fontsize, linewidth, color, alpha, linestyle, grid
):

    """
    grid_labels_lines(ax, xticks, yticks, fontsize, linewidth, color, alpha, linestyle)

    Function for placing x- and y-axis tick marks and grid lines for regional climatologies.

        Parameters
        ----------
        ax : geospatial axes for the subplot (cartopy object)
            ex: fig, axes = plt.subplots(3, 2, figsize=(16,12), subplot_kw={'projection': projection})
                ax1, ax2, ax3, ax4, ax5, ax6 = axes.flatten()
        xticks : List of longitudinal tick marks
            ex: xticks = [0, 60, 120, 180, -120, -60, -0]
        yticks : List of latitudinal tick marks
            ex: yticks = [-60, -30, 0, 30, 60]
        fontsize : Specifies the font size of the tickmarks on the x and y axes
            ex: fontsize = 20
        linewidth : Specifies linewidth for grid lines.
        color : Specifies color for grid lines.
        alpha : Specifies the degree of transparency of grid lines.
        linestyle : Specifies the line type of grid lines.
        grid : Specify if you want grid lines. True means grid is present.

        Returns
        -------
        Plots with gridline and and tick marks on the left and bottom axes.

        Libraries necessary to run function
        -----------------------------------
        import cartopy.crs as ccrs
        import cartopy.mpl.ticker as cticker

    """

    # import libraries:
    import cartopy.crs as ccrs
    import cartopy.mpl.ticker as cticker

    # set x and y tick labels
    ax.set_xticks(xticks, crs=ccrs.PlateCarree())
    ax.set_xticklabels(xticks, fontsize=fontsize)
    ax.set_yticks(yticks, crs=ccrs.PlateCarree())
    ax.set_yticklabels(yticks, fontsize=fontsize)

    # Place degrees symbol on x and y tick marks
    lon_formatter = cticker.LongitudeFormatter()
    lat_formatter = cticker.LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

    # Set grid lines
    if grid == True:
        ax.grid(linewidth=linewidth, color=color, alpha=alpha, linestyle=linestyle)


########### Colorbar Function ###########
def set_cbar(
    cs, cax, fig, orientation, extend, cbar_label, nbins, fontsize, cbar_ticks, task
):

    """
    set_cbar(cs, cax, fig, orientation, extend, cbar_label, nbins, fontsize, cbar_ticks, task)

        Function for placing a color bar on a plot. Two types of colorbar labels:
            1. Default color bar tick marks.
            2. Customized color bar tick marks.
        Many other colorbar keyword arguments can be found at:
        https://matplotlib.org/3.2.1/tutorials/colors/colorbar_only.html

        Parameters
        ----------
        cs : Map of data on subplot axis using cartopy projection.
            ex: cs = ax.pcolor(lon, lat, swh_phase, vmin=-np.pi, vmax=np.pi, cmap=cmo.phase, transform=projection)
        cax : Color bar axis with positioning vector of the colorbar with the folowing
              parameters: cax = plt.axes([left, bottom, width, height]).
            ex: cax = plt.axes([.47, .17, 0.01, 0.16])
        fig : Figure object which the colorbar will attached to.
            ex: fig, axes = plt.subplots(3, 2, figsize=(16,12),
                                         subplot_kw={'projection': projection})
        orientation : Specifies if the color bar is vertical or horizontal. Options for
                      keyword argument includes: orientation = 'horizontal' or
                      orientation = 'vertical'.
        extend : Specifies whether the colorbar will have extension towards high or low
                 values. Options include: extend = 'neither', 'both', 'min', or 'max'.
        cbar_label : Color bar label.
            ex: cbar_label = '$m$'
        fontsize : Fontsize of color bar label and tickmarks.
            ex: fontsize = 20
        nbins : Number of tick marks on colorbar axis
            ex: nbins = 5
        cbar_ticks : A list of tick marks that will be placed on colorbar (note that the
                     number of tick mark labels must be equal to the number of bins on color
                     bar)
            ex: cbar_ticks = [np.arange(-np.pi, np.pi+0.5, (np.pi + np.pi)/6).tolist(),[Jun,
                             Aug, October, Dec, Feb, Apr, June]]
        task : Specifies whether the colorbar will need to be modified with custom tick
               marks. Options include: task = 'custom ticks' or task = 'regular'.

        Returns
        -------
        Plots with colorbars in desired location and orientation.

        Libraries necessary to run function
        -----------------------------------
        from matplotlib import ticker

    """

    # import libraries
    from matplotlib import ticker

    # create colorbar for plot
    if task == "regular":
        cbar = fig.colorbar(cs, cax=cax, orientation=orientation, extend=extend)

        # set number of tick marks:
        tick_locator = ticker.MaxNLocator(nbins=nbins)
        cbar.locator = tick_locator
        cbar.update_ticks()
        cbar.set_ticks(cbar.locator)

    elif task == "custom ticks":
        cbar = fig.colorbar(
            cs, cax=cax, orientation=orientation, ticks=cbar_ticks[0], extend=extend
        )

    # case 1: vertical colorbar
    if orientation == "vertical":
        cbar.ax.set_ylabel("%s" % cbar_label, fontsize=fontsize)
        if task == "custom ticks":
            cbar.ax.set_yticklabels(cbar_ticks[1])

    # case 2: horizontal colorbar
    elif orientation == "horizontal":
        cbar.ax.set_xlabel("%s" % cbar_label, fontsize=fontsize)
        if task == "custom ticks":
            cbar.ax.set_xticklabels(cbar_ticks[1])

    # set the fontsize of colorbar tickmarks
    cbar.ax.tick_params(labelsize=fontsize)

    return


########### Latitude and Logitude Axis Labels Function ###########
def set_axes_label(ax, xdist_lat, ydist_lat, xdist_lon, ydist_lon, fontsize):

    """
    set_axes_label(ax, xdist_lat, ydist_lat, xdist_lon, ydist_lon, fontsize)

        Function for placing x and y axis labels for longitude and latitude respectively

        Parameters
        ----------
        ax : Geospatial axes for the subplot (cartopy object)
            ex: fig, axes = plt.subplots(3, 2, figsize=(16,12),
                                         subplot_kw={'projection': projection})
                ax1, ax2, ax3, ax4, ax5, ax6 = axes.flatten()
                ax = ax1 or ax2 or ...
        xdist_lat : Horizontal distance for latitude label referenced from right side of
                    figure
             ex: xdist_lat = -0.1
        ydist_lat : Vertical distance for latitude label referenced from bottom of figure
            ex: yticks = 0.50
        xdist_lon : Horizontal distance for longitude label referenced from right side of
                    figure
            ex: xdist_lon = 0.5
        ydist_lon : Vertical distance for longitude label referenced from bottom of figure
            ex: ydist_lon = -0.25
        fontsize : Fontsize of label


        Returns
        -------
        A geospatial map with axis labels on the left and bottom

        Libraries necessary to run function
        -----------------------------------
        import matplotlib.pyplot as plt
    """

    # y axis label
    ax.text(
        xdist_lat,
        ydist_lat,
        "Latitude",
        va="bottom",
        ha="center",
        rotation="vertical",
        rotation_mode="anchor",
        transform=ax.transAxes,
        fontsize=fontsize,
    )

    # x axis label
    ax.text(
        xdist_lon,
        ydist_lon,
        "Longitude",
        va="bottom",
        ha="center",
        rotation="horizontal",
        rotation_mode="anchor",
        transform=ax.transAxes,
        fontsize=fontsize,
    )

    return


########### Subplot label Function ###########
def subplot_label(
    ax, xdist_label, ydist_label, subplot_label, form, fs_shade, fs_main, color
):

    """
    subplot_label(ax, xdist_label, ydist_label, subplot_label, form, fs_shade, fs_main, color)

        Function for placing subplot labels for figures that will be used in research papers.
        Two types of labeling:
            1. Shading behind letter or number.
            2. Box behind letter or number.
        Features of the labeling:
            1. Black label with gray shading.
            2. Square box with 0.8 transparency and 1 linewidth.

        Parameters
        ----------
        ax : Geospatial axes for the subplot (cartopy object)
             ex: fig, axes = plt.subplots(3, 2, figsize=(16,12), subplot_kw={'projection': projection})
                 ax1, ax2, ax3, ax4, ax5, ax6 = axes.flatten()
                 ax = ax1 or ax2 or ...
        xdist_label : Horizontal distance for subplot label referenced from right side of
                      figure.
            ex: xdist_label = 0.2
        ydist_label : Vertical distance for subplot label referenced from bottom of figure
             ex: ydist_label = 0.8
        subplot_label : String of words for label
            ex: subplot_label = 'A'
        form : Specifies the format of the subplot label. Options for keyword argument:
               form = 'box' or 'shading'.
        fs_shade : Fontsize of shading label
            ex: fs_shade = 28
        fs_main : Fontsize of main label
            ex: fs_main = 18
        color : Specifies color of shading.

        Returns
        -------
        A geospatial map with a subplot label in specified location.

        Libraries necessary to run function
        -----------------------------------
        import matplotlib.pyplot as plt
    """

    # Shading label
    if form == "shading":
        ax.text(
            xdist_label,
            ydist_label,
            "%s" % subplot_label,
            va="center",
            ha="center",
            transform=ax.transAxes,
            fontsize=fs_shade,
            fontweight="bold",
            color=color,
            alpha=0.5,
        )
        ax.text(
            xdist_label,
            ydist_label,
            "%s" % subplot_label,
            va="center",
            ha="center",
            transform=ax.transAxes,
            fontsize=fs_main,
            fontweight="bold",
        )

    # Boxed label
    elif form == "box":
        ax.text(
            xdist_label,
            ydist_label,
            "%s" % subplot_label,
            va="center",
            ha="center",
            transform=ax.transAxes,
            fontsize=fs_main,
            bbox=dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=1, alpha=0.8),
            fontweight="bold",
        )
    return
