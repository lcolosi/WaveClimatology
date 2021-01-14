# Decorrelation Time Scale Functions
## Luke Colosi | lcolosi@ucsd.edu | December 22nd, 2020

########### Autocovariance Function ###########
def autocov(data, window, lag, task, bias):

    """
    autocov(data, window, lag, task, bias)

        Function for computing the autocovariance and autocorrelation function for positive and negative lag.

        Parameters
        ----------
        data : Time series of data that you wish to compute the autocovariance function for. This data must
               be preprocessed which means:
                   1) Detrend to remove any signal that are undesirable in autocovariance function
                   2) Missing data gaps are filled with masked values in order to have a continuous time series
                   3) Flagged data should be replace with mask values.
               Data must be a masked array!
        window : Time interval which the fixed time series will be set to for computing autocovariance. Options
                 include the full time series or varying lengths from the begining of the time series to a
                 desired point (i.e. day, week, month, etc.). This depends on the sampling frequency on your data.
                     ex: window = len(data) (full time series) or window = 7 (week when data collected daily)
        lag : The desired amount of lags where the correlation is coputed for. The specified amount of lags is
              highly dependent on the window length of the time series. You want to set the amount of lags to
              a value where the correlation coefficent is for the proper amount of iterations along to fixed
              time series.
                     ex: lag_dt = len(data) (compute correlation coefficient at lag decreasing by one measurement)
                                             at a time)
        task : Specifies whether you want to compute the autocovariance function for the entire time series or
               just windows of the time series. Options include;
                   task = 'full' or 'window'
        bias : Specifies where the covariance is unbiased (normalized by 1/n-m) or biased (normalized by 1/n)
               where n = total number of data points in time series and m is the lag time from 0 lag. Furthermore,
               specifies whether the correlation coefficent is biased or unbaised using the same normalizations in
               numerator (unbiased (normalized by 1/n-m) or biased (normalized by 1/n)) and the normalization 1/n for
               both cases in the demominator.
                   ex : bias = 'biased' or 'unbiased'

        Returns
        -------
        coef_pos : positive lag autocorrelation function
        coef_neg : negative lag autocorrelation function
        cov_pos : positive lag autocovariance function
        cov_neg : negative lag autocovariance function
        upper_95_CI: Upper 95% confidence interval value
        lower_95_CI: Lower 95% confidence interval value

        Libraries necessary to run function
        -----------------------------------
        import numpy as np

    """

    # import libraries:
    import numpy as np

    # Choose interval length n which the correlation coefficient will be computed
    n = window

    # set fixed data segments
    if task == "window":
        fix = data[:n]

    # set autocorrelation functions
    coef = np.empty((lag,))
    cov = np.empty((lag,))
    upper_95_CI = np.empty((lag,))
    lower_95_CI = np.empty((lag,))

    # loop through each lag time to compute the correlation coefficient
    for i in range(lag):

        # set running and fixed window for full time series case
        if task == "full":
            # set running data segment fowards
            running = data[i : i + n]

            # set running data segment backwards
            fix = data[: n - i]

        # Set running window for window case:
        elif task == "window":
            # set running data segment fowards
            ind_run_i = i * n
            ind_run_f = ind_run_i + n
            running = data[ind_run_i:ind_run_f]

        # Remove mean from each segment before computing covariance and correlation:
        fix = fix - np.ma.mean(data)
        running = running - np.ma.mean(data)

        # Compute correlation coefficient terms at lag time i
        inner_product = np.ma.dot(fix.T, running)

        # compute correlation coefficient and covariance
        # case 1: unbiased
        if bias == "unbiased":
            coef[i] = ((1 / len(running)) * inner_product) / (
                (1 / n) * ((data - np.ma.mean(data)) ** 2).sum()
            )
            cov[i] = (1 / len(running)) * inner_product
        # case 2: biased
        elif bias == "biased":
            coef[i] = inner_product / ((data - np.ma.mean(data)) ** 2).sum()
            cov[i] = (1 / n) * inner_product

        # compute number of data point used in the correlation coefficient calculation at each lag time:
        npoints_nonmasked = running + fix
        n_eff = np.count_nonzero(npoints_nonmasked)

        # Compute upper and lower 95% confidence interval values using: delta_r = erf^{-1}(p = 0.95) * sqrt{2/N}
        upper_95_CI[i] = 1.96 / np.sqrt(n_eff)
        lower_95_CI[i] = -1.96 / np.sqrt(n_eff)

    # Create positive and negative lag autocorrelation and autocovariance function
    coef_pos = coef
    coef_neg = coef[: -len(coef) : -1]
    cov_pos = cov
    cov_neg = cov[: -len(cov) : -1]

    return coef_pos, coef_neg, cov_pos, cov_neg, upper_95_CI, lower_95_CI


########### Decorrelation Time Scale Function ###########
def decor_scale(data, window, lag, method, dt, units):

    """
    decor_scale(data, window, lag, method, dt, units)

        Function for computing decorrelation time scales for a time series.

        Parameters
        ----------
        data : 1D time series of data assuming time gaps are filled with masked values. Signal that would cause any
               unwanted correlation are removed (demeaning the data is preformed in function). Data must be
               a masked array.
        window : Time interval which the fixed time series will be set to compute
                 autocovariance and autocorrelation. Options include the full time series or varying
                 lengths from the begining of the time series to a desired point (i.e. day, week, month,
                 etc.). This depends on the sampling frequency on your data.
            ex: window = len(data) (full time series) or window = 7 (week when data
                         collected daily)
        lag : The desired amount of lags where the correlation is computed. The specified
              amount of lags is highly dependent on the window length of the time series.
            ex: lag_dt = len(data) (compute correlation coefficient at lag decreasing
                by one measurement at a time)
        method : Specifies the method in which the decorrelation time scale is computed.
                 Options include: method = 'zero_crossing', '95_sig_level', 'half_autocor_unbias',
                 'integral_bias_coef', 'integral_unbiased_coef'. Integral time scale estimates follow SIO 221C dof
                 lecture note. 'zero_crossing', '95_sig_level', and 'half_autocor_unbias' methods output a integer.
        dt : time interval between observations.
                ex: dt = 1 #units: day
        units : Specifies the units which the original data has and desired final units of the
                decorrelation time scale. Options include:
                    units = [data_units, decor_units]
                where data_units are the time steps between observations.
                    ex: data_units = 'sec' and decor_units = 'day'

        Returns
        -------
        ds : decorrelation time scale of the data in units specified above.
        ds_N : decorrelation time scale as a function of N
        coef_pos: Positive lag autocorrelation function
        coef_neg: Negative lag autocorrelation function
        upper_95_CI: Upper 95% confidence interval value
        lower_95_CI: Lower 95% confidence interval value

        Libraries necessary to run function
        -----------------------------------
        Numpy : import numpy as np
        autocovariance function: from autocovariance_temporal import autocov

    """

    # Set path to autocovariance function:
    import sys

    sys.path.append("/zdata/home/lcolosi/python_functions/")

    # Import libraries:
    import numpy as np
    from autocovariance_temporal import autocov

    # demean time series:
    data_dm = data - np.ma.mean(data)

    # Compute autocovariance:

    if method == "integral_unbiased_coef" or method == "half_autocor_unbias":
        coef_pos, coef_neg, cov_pos, cov_neg, upper_95_CI, lower_95_CI = autocov(
            data_dm, window, lag, task="full", bias="unbiased"
        )
    else:
        coef_pos, coef_neg, cov_pos, cov_neg, upper_95_CI, lower_95_CI = autocov(
            data_dm, window, lag, task="full", bias="biased"
        )
    # case 1: zero crossing
    if method == "zero_crossing":

        # set variables:
        ds_N = None

        # compute the indices of all zero crossing and pick indices corresponding to the
        # point before crossing the x-axis
        icross = np.diff(np.sign(coef_pos)) == -2

        # compute decorrelation scale
        ds = np.where(icross == True)[0][0] * dt

    # case 2: Statistically significant correlation coefficient threshold method (at 95% significance level)
    if method == "95_sig_level":

        # set variables:
        ds_N = None

        # compute the index of the first crossing of the 95% confidence interval
        icross = np.diff(np.sign(coef_pos - upper_95_CI)) == -2

        # compute decorrelation scale
        ds = np.where(icross == True)[0][0] * dt

    # case 3: 2 times the time it takes lagged covariance p(n) to reach p(0)/2
    if method == "half_autocor_unbias":

        # set variables:
        ds_N = None

        # compute the index of the first crossing of the covariance function (p(n)) when p(n) = p(0)/2
        icross = np.diff(np.sign(coef_pos - (coef_pos[0] / 2))) == -2

        # compute decorrelation scale as 2 times the time is takes p(n) = p(0)/2
        ds = 2 * np.where(icross == True)[0][0] * dt

    # case 4: Interal time scale method using biased autocorrelation (Maximum Area underneath the autocorrelation curve)
    if method == "integral_bias_coef":

        # Initialize variable:
        ds_N = np.zeros((len(coef_pos),))

        # compute decorrelation scale for a range of N values:
        for i_N in range(1, len(coef_pos) + 1):

            # Combine positive and negative lag to one array
            if i_N == 1:
                coef = coef_pos[:i_N]
            else:
                coef = np.hstack((coef_neg[-(i_N - 1) :], coef_pos[:i_N]))

            # compute decorrelation scale by suming i_N points
            ds_i = np.nansum(coef) * dt

            # Save ith decorrelation scale:
            ds_N[(i_N - 1)] = ds_i

        # set the decorrelation time scale to the maximum of ds_N:
        ds = np.nanmax(ds_N)

    # case 5: Interal time scale method using biased or normalized autocovariance (Maximum Area underneath the autocorrelation curve)
    if method == "integral_unbiased_coef":

        # Initialize variable:
        ds_N = np.zeros((len(cov_pos),))

        # compute decorrelation scale for a range of N values:
        for i_N in range(1, len(coef_pos) + 1):

            # Combine positive and negative lag to one array and compute decorrelation time scale
            if i_N == 1:
                coef = coef_pos[:i_N]
                ds_i = 1
            else:
                coef = np.hstack((coef_neg[-(i_N - 1) :], coef_pos[:i_N]))
                l = len(coef_pos[: (i_N - 1)])
                n = np.arange(-l, l + 1, 1)
                ds_i = np.nansum((1 - (abs(n) / l)) * coef)

            # Save ith decorrelation scale:
            ds_N[(i_N - 1)] = ds_i

        # set the decorrelation time scale to the maximum of ds_N:
        ds = np.nanmax(ds_N)

        # Convert units:
        if units[0] != units[1]:

            # secs to days
            if units[0] == "secs" and units[1] == "days":
                ds = ds * (1 / 60 * 60 * 24)

            # mins to days
            elif units[0] == "mins" and units[1] == "days":
                ds = ds * (1 / 60 * 24)

            # days to hours
            elif units[0] == "days" and units[1] == "hours":
                ds = ds * (24)

    return ds, ds_N, coef_pos, coef_neg, upper_95_CI, lower_95_CI
