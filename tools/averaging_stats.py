# Monthly, Climatological and Running mean functions and Statistical Moment Function
## Luke Colosi | lcolosi@ucsd.edu | December 22nd, 2020

########### Monthly Averaging Function ###########
def monthly_average(date_time, data):

    """
    monthly_average(date_time, data)

        Function to compute monthly means.

        Parameters
        ----------
        data : Numpy 1D or 2D masked array of data.
        date_time : Numpy 2D array of datetime array.

        Returns
        -------
        monthly_data : A dictionary contain the folloing key variables (all are lists type):
            a) 'time' : time step of the month averaged
            b) 'data' : month data 3D array
            c) 'mean' : data monthly mean
            d) 'median' : data monthly median
            e) 'std' : data monthly standard deviation
            f) 'var' : data monthly variance
            g) 'N' : amount of observations averaged over (used for calculating the standard
                     deviation of the mean)

        Libraries necessary to run function
        -----------------------------------
        Numpy : import numpy as np
        datetime : import datetime

        Function is originally created by Bia Villas Bôas with small changes by Luke Colosi.
    """

    # import library:
    import datetime
    import numpy as np

    # make an assertion that date_time variable must be an array:
    assert isinstance(date_time, np.ndarray), "date_time should be a numpy array"

    # set year and month 2D time arrays
    years = np.array([y.year for y in date_time])
    months = np.array([m.month for m in date_time])

    # initialize dictionary
    monthly_data = {}
    monthly_data = {
        "time": [],
        "data": [],
        "mean": [],
        "median": [],
        "std": [],
        "var": [],
        "N": [],
    }

    # initialize a year loop
    for year in np.unique(years):
        for month in np.unique(months):

            # initialize time indices
            ind_year = years == year
            ind_month = months == month
            ind = ind_year * ind_month

            # apply indice to swh and time arrays and calculate mean, std, median, and nobs
            # using a try block
            try:
                # Call data and time
                # Case 1: 1d array
                if data.ndim == 1:
                    tmp = data[ind]
                # Case 2: 3d array
                elif data.ndim == 3:
                    tmp = data[ind, :, :]
                time = date_time[ind]

                # Set data and time variable
                delta_t = datetime.timedelta(
                    seconds=np.mean([(t - time[0]).total_seconds() for t in time])
                )
                monthly_data["time"].append(time[0] + delta_t)
                monthly_data["data"].append(tmp)

                # compute mean, std, median, and nobs
                mean = np.ma.mean(tmp, axis=0)
                std = np.ma.std(tmp, ddof=1, axis=0)
                var = np.ma.var(tmp, ddof=1, axis=0)
                median = np.ma.median(tmp, axis=0)

                # Mask all mean values that are compute with one data point (std is undefined)
                mean_n = np.ma.masked_where(np.ma.getmask(std), mean)

                # save data in dictionary
                monthly_data["mean"].append(mean_n)
                monthly_data["median"].append(median)
                monthly_data["std"].append(std)
                monthly_data["var"].append(var)
                monthly_data["N"].append(tmp.count(axis=0))

            except:
                pass

    return monthly_data


########### Climatological Averaging Function ###########
def clima_mean(date_time, data):

    """
    clima_mean(date_time, data)

        Function to computes a monthly climatology.

        Parameters
        ----------
        data : Numpy 1D or 2D masked array of data.
        date_time : Numpy 2D array of datetime array.

        Returns
        -------
        monthly_data : A dictionary contain the folloing key variables (all are lists type):
            a) 'month' : time step of the month averaged
            b) 'data' : month data 3D array
            c) 'mean' : data monthly mean
            d) 'median' : data monthly median
            e) 'std' : data monthly standard deviation
            f) 'var' : data monthly variance
            g) 'N' : amount of observations averaged over (used for calculating the standard
                     error of the mean)

        Libraries necessary to run function
        -----------------------------------
        Numpy: import numpy as np
    """

    # import library
    import numpy as np

    # make an assertion that date_time variable must be an array:
    assert isinstance(date_time, np.ndarray), "date_time should be a numpy array"

    # set month 2D time arrays
    months = np.array([m.month for m in date_time])

    # initialize dictionary
    monthly_data = {}
    monthly_data = {
        "month": [],
        "data": [],
        "mean": [],
        "median": [],
        "std": [],
        "var": [],
        "N": [],
    }

    # initialize a monthloop
    for m in range(1, 13):

        # initialize time indices
        ind = months == m

        # Call data and time
        # Case 1: 1d array
        if data.ndim == 1:
            tmp = np.ma.array(data[ind])
        # Case 2: 3d array
        elif data.ndim == 3:
            tmp = np.ma.array(data[ind, :, :])
        time = date_time[ind]

        # Set data and time variable
        monthly_data["data"].append(tmp)
        monthly_data["month"].append(m)

        # compute the mean, median, std, and the number of observations for each month:
        monthly_data["mean"].append(np.ma.mean(tmp, axis=0))
        monthly_data["median"].append(np.ma.median(tmp, axis=0))
        monthly_data["std"].append(np.ma.std(tmp, ddof=1, axis=0))
        monthly_data["var"].append(np.ma.var(tmp, ddof=1, axis=0))
        monthly_data["N"].append(tmp.count(axis=0))

    return monthly_data


########### Statistical Moments Function ###########
def stat_moments_temporal(data, date_time, time_int, sample_size, skew_kurt=False):

    """
    stat_moments_monthly(date, data_time, time_int, sample_size, skew_kurt=False)

        Function to compute the first four statistical moments for each month or each season of the year. Statistical moments include mean, variance, skewness, and kurtosis.

        Parameters
        ----------
        data : Numpy 3d or 1d masked array with time in the first index of the array.
        date_time : Numpy 1d array for the time series of data with datetime valued elements.
        time_int : Specifies whether statistical moments are computed seasonally or monthly. Options include: time_int = 'monthly' or 'seasonally'.
        sample_size : Specifies whether the data set represents the entire population or a sample of the population. Options include: sample_size = 'population' or 'sample'.
        skew_kurt : Specifies if skewness or kurtosis are computed. Options inlcude skew_kurt = True or False. Default is skew_kurt=False.

        Returns
        -------
        stats : A dictionary containg:
                    1. Data partitioned monthly or seasonally
                    2. Time (numerical value corresponding to the month or season that the statistical moment is computed)
                    3. Number of observations
                    2. Mean
                    3. Variance (Documentation: https://en.wikipedia.org/wiki/Variance)
                    4. Skewness (Documentation: https://en.wikipedia.org/wiki/Skewness)
                    5. Excess Kurtosis (Documentation: https://en.wikipedia.org/wiki/Kurtosis)

        Libraries necessary to run function
        ---------------------------
        Numpy: import numpy as np
    """

    # import library
    import numpy as np

    # set month 2D time arrays
    months = np.array([m.month for m in date_time])

    # initialize the dictionary
    stats = {}
    stats = {
        "data": [],
        "time": [],
        "N": [],
        "mean": [],
        "var": [],
        "skew": [],
        "kurt": [],
    }

    ####### Monthly Statistics #######
    if time_int == "monthly":

        # initialize a month loop
        for m in range(1, 13):

            # initialize time indices
            ind = months == m

            # Call data
            # Case 1: 1d array
            if data.ndim == 1:
                tmp = data[ind]
            # Case 2: 3d array
            elif data.ndim == 3:
                tmp = data[ind, :, :]

            # save data, time, and nobs variables
            stats["data"].append(tmp)
            stats["time"].append(m)
            stats["N"].append(tmp.count(axis=0))

            # set number of observations:
            n = tmp.count(axis=0)

            ####### Sample of Population #######
            if sample_size == "sample":

                # compute sample mean and variance
                monthly_mean = np.ma.sum(tmp, axis=0) / n
                monthly_var = (1 / (n - 1)) * np.ma.sum(
                    (tmp - monthly_mean) ** 2, axis=0
                )

                # compute third, and fourth moments:
                if skew_kurt == True:
                    m2 = (1 / n) * np.ma.sum((tmp - monthly_mean) ** 2, axis=0)
                    m3 = (1 / n) * np.ma.sum((tmp - monthly_mean) ** 3, axis=0)
                    m4 = (1 / n) * np.ma.sum((tmp - monthly_mean) ** 4, axis=0)

                    # compute sample skewness, and kurtosis:
                    monthly_skew = (m3) / ((monthly_var) ** (3 / 2))
                    monthly_kurt = ((m4) / ((m2) ** 2)) - 3

            ####### Full Population #######
            elif sample_size == "population":

                # compute first moment:
                monthly_mean = np.ma.sum(tmp, axis=0) / n

                # Compute second moment:
                # compute second power of data
                tmp_square = tmp ** 2

                # sum powers of data and divide by nobs:
                monthly_var_mean = np.ma.sum(tmp_square, axis=0) / n

                # compute population variance:
                monthly_var = monthly_var_mean - monthly_mean ** 2

                # compute third and fourth moments:
                if skew_kurt == True:
                    # compute the third and fourth power of data
                    tmp_cube = tmp ** 3
                    tmp_quad = tmp ** 4

                    # sum powers of data and divide by nobs:
                    monthly_skew_mean = np.ma.sum(tmp_cube, axis=0) / n
                    monthly_kurt_mean = np.ma.sum(tmp_quad, axis=0) / n

                    # compute population skewness and kurtosis:
                    monthly_std = (monthly_var) ** (1 / 2)
                    monthly_skew = (
                        monthly_skew_mean
                        - 3 * monthly_mean * monthly_var
                        - monthly_mean ** 3
                    ) / (monthly_std) ** 3
                    monthly_kurt = (
                        (
                            monthly_kurt_mean
                            - 4 * monthly_mean * monthly_skew_mean
                            + 6 * monthly_var_mean * monthly_mean ** 2
                            - 3 * monthly_mean ** 4
                        )
                        / (monthly_std ** 4)
                    ) - 3

            # Save mean and variance:
            stats["mean"].append(monthly_mean)
            stats["var"].append(monthly_var)

            # Save skewness, and kurtosis:
            if skew_kurt == True:
                stats["skew"].append(monthly_skew)
                stats["kurt"].append(monthly_kurt)

    ####### Seasonal Statistics #######
    elif time_int == "seasonally":

        # initialize variables
        seasons = [[12, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]]

        # initialize a seasonal loop
        for s in range(0, 4):

            # call season:
            season = seasons[s]

            # initialize time indices
            bool_1 = months == season[0]
            bool_2 = months == season[1]
            bool_3 = months == season[2]
            ind_1 = bool_1 | bool_2
            ind = ind_1 | bool_3

            # Call data
            # Case 1: 1d array
            if data.ndim == 1:
                tmp = data[ind]
            # Case 2: 3d array
            elif data.ndim == 3:
                tmp = data[ind, :, :]

            # save data, time, and nobs variables
            stats["data"].append(tmp)
            stats["time"].append(s)
            stats["N"].append(tmp.count(axis=0))

            # set number of observations:
            n = tmp.count(axis=0)

            ####### Sample of Population #######
            if sample_size == "sample":

                # compute sample mean and variance
                season_mean = np.ma.sum(tmp, axis=0) / n
                season_var = (1 / (n - 1)) * np.ma.sum((tmp - season_mean) ** 2, axis=0)

                # compute second, third, and fourth moments:
                if skew_kurt == True:
                    m2 = (1 / n) * np.ma.sum((tmp - season_mean) ** 2, axis=0)
                    m3 = (1 / n) * np.ma.sum((tmp - season_mean) ** 3, axis=0)
                    m4 = (1 / n) * np.ma.sum((tmp - season_mean) ** 4, axis=0)

                    # compute sample skewness, and kurtosis:
                    season_skew = (m3) / ((season_var) ** (3 / 2))
                    season_kurt = ((m4) / ((m2) ** 2)) - 3

            ####### Full Population #######
            elif sample_size == "population":

                # compute first moment:
                season_mean = np.ma.sum(tmp, axis=0) / n

                # Compute second moment:
                # compute second power of data
                tmp_square = tmp ** 2

                # sum powers of data and divide by nobs:
                season_var_mean = np.ma.sum(tmp_square, axis=0) / n

                # compute population variance:
                season_var = season_var_mean - season_mean ** 2

                # compute third and fourth moments:
                if skew_kurt == True:
                    # compute the third and fourth power of data
                    tmp_cube = tmp ** 3
                    tmp_quad = tmp ** 4

                    # sum powers of data and divide by nobs:
                    season_skew_mean = np.ma.sum(tmp_cube, axis=0) / n
                    season_kurt_mean = np.ma.sum(tmp_quad, axis=0) / n

                    # compute variance, skewness, and kurtosis:
                    season_std = (season_var) ** (1 / 2)
                    season_skew = (
                        season_skew_mean
                        - 3 * season_mean * season_var
                        - season_mean ** 3
                    ) / (season_std) ** 3
                    season_kurt = (
                        (
                            season_kurt_mean
                            - 4 * season_mean * season_skew_mean
                            + 6 * season_var_mean * season_mean ** 2
                            - 3 * season_mean ** 4
                        )
                        / (season_std ** 4)
                    ) - 3

            # Save mean and variance:
            stats["mean"].append(season_mean)
            stats["var"].append(season_var)

            # Save skewness, and kurtosis:
            if skew_kurt == True:
                stats["skew"].append(season_skew)
                stats["kurt"].append(season_kurt)

    return stats


########### Rnning Mean (Box Car Filter) Function ###########
def running_mean(data, k_dim, task):

    """
    running_mean(data, k_dim, task)

        Function for computing the running mean of 1 or 2 dimensional arrays.

        Parameters
        ----------
        data : Numpy masked array of any 2d or 1d dimensional size. If the data set is 1d, then the data must be a
               column vector.
        k_dim : Dimensions of the kernal matrix in list format. If the data set is 1d, the kernal dimension must be
                k_dim = [ndim,1] such that the kernal matrix will slide down the coloumn vector.
                    ex: k_dim = [4, 4]
        task : Specifies the purpose of the box car filter. Options for this input include task = 'running_mean'
               (used for smoothing out data sets) and task = 'deresolve' (used for decreasing the resolution of a
               data set or image)

        Returns
        -------
        data_rm : 1D or 2D numpy array with full convolution sampled using running mean (same output as output using
                  valid flag with np.convolve or scipy.signal.convolve2d) or deresolved schemes.
        w_conv : Full linear convolution (same output as output using full flag with np.convolve or
                 scipy.signal.convolve2d)

        Libraries necessary to run function
        -----------------------------------
        Scipy: from scipy import signal
        Numpy: import numpy as np
    """

    # import libraries
    from scipy import signal
    import numpy as np

    # Check if data is a masked array
    assert type(data) == np.ma.core.MaskedArray, "Data is not a masked array"

    # create kernal matrix:
    w = np.ones((k_dim[0], k_dim[1]))

    # Normalize kernal matrix:
    w_norm = w / np.sum(w)

    ######### 1d convolution #########
    if k_dim[1] == 1:

        # flatten kernal matrix and data:
        w_norm = w_norm.flatten()
        data = data.flatten()

        # convolve data and kernal matrix:
        w_conv = np.ma.convolve(data, w_norm, mode="full", propagate_mask=False)

        # case 1: deresolve
        if task == "deresolve":
            # extract lower resolution averaged elements in w_conv matrix:
            data_rm = w_conv[(k_dim[0] - 1) :: k_dim[0]]

        # case 2: running mean (same output as valid)
        elif task == "running_mean":
            # extract running mean elements in w_conv matrix:
            data_rm = w_conv[(k_dim[0] - 1) : len(w_conv) - k_dim[0] + 1]

    ######### 2d convolution #########
    elif k_dim[1] > 1:

        # convolve data and kernal matrix:
        w_conv = signal.convolve2d(data, w_norm)
        w_row, w_coln = w_conv.shape

        # case 1: deresolve
        if task == "deresolve":
            # extract lower resolution averaged elements in w_conv matrix:
            data_rm = w_conv[(k_dim[0] - 1) :: k_dim[0], (k_dim[1] - 1) :: k_dim[1]]

        # case 2: running mean
        elif task == "running_mean":

            # extract running mean elements in w_conv matrix:
            data_rm = w_conv[
                (k_dim[0] - 1) : w_row - k_dim[0] + 1,
                (k_dim[1] - 1) : w_coln - k_dim[1] + 1,
            ]

    return data_rm, w_conv
