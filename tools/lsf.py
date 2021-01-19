# Least Squares-Fit functions
## Luke Colosi | lcolosi@ucsd.edu | December 22nd, 2020

########### Detrend Function ###########
def detrend(data):

    """
    detrend(data)

        Function for removing a linear trend from a 1 dimensional data series, but retains the mean.

        Parameters
        ----------
        data : signal you wish to detrend

        Returns
        -------
        data_detrend : detrended signal

        Libraries necessary to run function
        -----------------------------------
        import numpy as np
        from unweighted_least_square_fit import least_square_fit

    """

    # import libraries
    import numpy as np
    from lsf import least_square_fit

    # Check if data is a 1 dimensional array
    assert data.ndim == 1, "Data is not a one dimensional array"

    # fit a linear trend at grid point:
    data_trend, x_trend = least_square_fit(
        data, trend="linear", parameters=2, period=12
    )

    # initialize time vector and linear trend:
    time = np.arange(1, len(data) + 1, 1)
    linear_trend = x_trend[0] + x_trend[1] * time

    # remove linear trend:
    data_detrend = data - linear_trend + np.mean(linear_trend)

    return data_detrend


########### Unweighted least squares fit Function ###########
def least_square_fit(data, trend, parameters, period):

    """
    least_square_fit(data, trend, parameters, period)

        Function to compute the unweighted least square fit of temporal data.

        Parameters
        ----------
        data : numpy 1D array of temporal data
        trend : Type of function least square is fitting. Options for this variables include:
               a) trend = 'linear'
               b) trend = 'sinusoidal'
        parameters : number of parameters in model
               a) if trend = 'linear', then options for parameters includes:
                    i) parameters = 1 : fits mean
                    ii) parameters = 2 : fits mean and linear trend
               b) if trend = 'sinusiodal', then options for paramters includes:
                    i) parameters = 3 : fits mean, and annual cycle (period = 1 year, 12 months, or 365.25 days)
                    ii) parameters = 4 : fits mean, linear trend, and annual cycle
                    iii) parameters = 5 : fits mean, annual cycle, and semi-annual cycle (period= 1/2 year,6 months,182.625 days)
                    iv) paramters = 6: fits mean, linear trend, annual cycle, and semi-annual cycle
        period : period of the sinusoidal signal

        Returns
        -------
        hfit: Unweighted least-squares fit model
        x_data: Coefficient parameters of the least-squares fit

        Libraries necessary to run function
        -----------------------------------
        Numpy : import numpy as np
        from scipy.linalg import inv
    """

    # import libraries:
    import numpy as np
    from scipy.linalg import inv

    # Check if data is a masked array
    assert type(data) == np.ma.core.MaskedArray, "Data is not a masked array"

    # initialize time vector:
    time = np.arange(1, len(data) + 1, 1)

    # consider for LSF only unmasked data:
    # case 1: Mask values exist
    if np.size(data.mask) > 1:

        # set mask
        ind = data.mask

        # Remove masked data points from time and data:
        time_n = time[~ind]
        data_n = data[~ind]

    # No masked values exist
    elif np.size(data.mask) == 1:

        # set new variables
        time_n = np.ma.copy(time)
        data_n = np.ma.copy(data)

    # set number of parameters and length of data:
    N, M = parameters, len(data_n)

    # create model matrix:
    A = np.zeros((M, N))

    # set conditional statements for each case that the lsf will be applied:
    if trend == "linear":

        # set the two cases for the linear trend parameters:
        if parameters == 1:
            # place parameters in model matrix:
            A[:, 0] = 1.0
            # Compute solution to the matrix equation Ax = b:
            x_data = np.dot(
                np.dot(inv(np.dot(A.conj().T, A)), A.conj().T), data_n.conj().T
            )
            # Compute unweighted least square fit model:
            hfit = x_data[0]

        elif parameters == 2:
            # place parameters in model matrix:
            A[:, 0] = 1.0
            A[:, 1] = time_n
            # Compute solution to the matrix equation Ax = b:
            x_data = np.dot(
                np.dot(inv(np.dot(A.conj().T, A)), A.conj().T), data_n.conj().T
            )
            # Compute unweighted least square fit model:
            hfit = x_data[0] + x_data[1] * time

    elif trend == "sinusoidal":

        # set the two cases for the linear trend parameters:
        if parameters == 3:
            # place parameters in model matrix:
            A[:, 0] = 1.0
            A[:, 1] = np.sin(time_n * 2 * np.pi / period)
            A[:, 2] = np.cos(time_n * 2 * np.pi / period)
            # Compute solution to the matrix equation Ax = b:
            x_data = np.dot(
                np.dot(inv(np.dot(A.conj().T, A)), A.conj().T), data_n.conj().T
            )
            # Compute unweighted least square fit model:
            hfit = (
                x_data[0]
                + x_data[1] * np.sin(time * 2 * np.pi / period)
                + x_data[2] * np.cos(time * 2 * np.pi / period)
            )

        elif parameters == 4:
            # place parameters in model matrix:
            A[:, 0] = 1.0
            A[:, 1] = time_n
            A[:, 2] = np.sin(time_n * 2 * np.pi / period)
            A[:, 3] = np.cos(time_n * 2 * np.pi / period)
            # Compute solution to the matrix equation Ax = b:
            x_data = np.dot(
                np.dot(inv(np.dot(A.conj().T, A)), A.conj().T), data_n.conj().T
            )
            # Compute unweighted least square fit model:
            hfit = (
                x_data[0]
                + x_data[1] * time
                + x_data[2] * np.sin(time * 2 * np.pi / period)
                + x_data[3] * np.cos(time * 2 * np.pi / period)
            )

        elif parameters == 5:
            # place parameters in model matrix:
            A[:, 0] = 1.0
            A[:, 1] = np.sin(time_n * 2 * np.pi / period)
            A[:, 2] = np.cos(time_n * 2 * np.pi / period)
            A[:, 3] = np.sin(time_n * 2 * np.pi / (period / 2))
            A[:, 4] = np.cos(time_n * 2 * np.pi / (period / 2))
            # Compute solution to the matrix equation Ax = b:
            x_data = np.dot(
                np.dot(inv(np.dot(A.conj().T, A)), A.conj().T), data_n.conj().T
            )
            # Compute unweighted least square fit model:
            hfit = (
                x_data[0]
                + x_data[1] * np.sin(time * 2 * np.pi / period)
                + x_data[2] * np.cos(time * 2 * np.pi / period)
                + x_data[3] * np.sin(time * 2 * np.pi / (period / 2))
                + x_data[4] * np.cos(time * 2 * np.pi / (period / 2))
            )

        elif parameters == 6:
            # place parameters in model matrix:
            A[:, 0] = 1.0
            A[:, 1] = time_n
            A[:, 2] = np.sin(time_n * 2 * np.pi / period)
            A[:, 3] = np.cos(time_n * 2 * np.pi / period)
            A[:, 4] = np.sin(time_n * 2 * np.pi / (period / 2))
            A[:, 5] = np.cos(time_n * 2 * np.pi / (period / 2))
            # Compute solution to the matrix equation Ax = b:
            x_data = np.dot(
                np.dot(inv(np.dot(A.conj().T, A)), A.conj().T), data_n.conj().T
            )
            # Compute unweighted least square fit model:
            hfit = (
                x_data[0]
                + x_data[1] * time
                + x_data[2] * np.sin(time * 2 * np.pi / period)
                + x_data[3] * np.cos(time * 2 * np.pi / period)
                + x_data[4] * np.sin(time * 2 * np.pi / (period / 2))
                + x_data[5] * np.cos(time * 2 * np.pi / (period / 2))
            )

    return hfit, x_data


########### Weighted least squares fit Function ###########
def weighted_least_square_fit(data, sigma, trend, parameters, period, phase=None):

    """
    weighted_least_square_fit(data, sigma, trend, parameters, period, phase)

        Function to compute the weighted least square fit of temporal data.

        Parameters
        ----------
        data : Numpy 1D array of temporal data
        sigma : Numpy 1D array of temporal Uncertainty (standard deviation or standard error of the mean). This
                program assumes that uncertainties are independent and uncorrelated.
        trend : Type of function least square is fitting. Options for this variables include:
               a) trend = 'linear'
               b) trend = 'sinusoidal'
        parameters : Number of parameters in model
               a) if trend = 'linear', then options for parameters includes:
                    i) parameters = 1 : fits mean
                    ii) parameters = 2 : fits mean and linear trend
               b) if trend = 'sinusiodal', then options for paramters includes:
                    i) parameters = 2: fits mean and amplitude of annual cycle with predetermined phase.
                    i) parameters = 3 : fits mean, and annual cycle (period = 1 year, 12 months, or 365.25 days)
                    ii) parameters = 4 : fits mean, linear trend, and single cycle of chosen frequency.
                    iii) parameters = 5 : fits mean, annual cycle, and semi-annual cycle (period= 1/2 year,6 months,182.625 days)
                    iv) paramters = 6: fits mean, linear trend, annual cycle, and semi-annual cycle
        period : Period of the sinusoidal signal
        phase : Specifies the fitted phased value used in the Sinusoidal trend with 2 parameters. Default value None.

        Returns
        -------
        hfit: Weighted least-squares fit model
        x_data: Coefficient parameters of the least-squares fit
        x_data_sigma: Uncertainties of the coefficient parameters of the least-squares fit

        Libraries necessary to run function
        -----------------------------------
        import numpy as np
        from scipy.linalg import inv

    """

    # import libraries:
    import numpy as np
    from scipy.linalg import inv

    # Check if data is a masked array
    assert type(data) == np.ma.core.MaskedArray, "Data is not a masked array"

    # initialize time vector:
    time = np.arange(1, len(data) + 1, 1)

    # consider for LSF only unmasked data:
    # case 1: Mask values exist
    if np.size(data.mask) > 1:

        # set mask
        ind = data.mask

        # Remove masked data points from time and data:
        time_n = time[~ind]
        data_n = data[~ind]
        sigma_n = sigma[~ind]

    # No masked values exist
    elif np.size(data.mask) == 1:

        # set new variables
        time_n = np.ma.copy(time)
        data_n = np.ma.copy(data)
        sigma_n = np.ma.copy(sigma)

    # set number of parameters and length of data:
    N, M = parameters, len(data_n)

    # create model matrix:
    A = np.zeros((M, N))

    # Create weighted matrix
    W = np.zeros((M, M))
    for isigma in range(M):
        W[isigma, isigma] = 1 / sigma_n[isigma]

    # Weight data column vector:
    data_w = np.dot(W, data_n)

    # set conditional statements for each case that the lsf will be applied:
    if trend == "linear":

        # set the two cases for the linear trend parameters:
        if parameters == 1:
            # place parameters in model matrix:
            A[:, 0] = 1.0
            # Weight model matrix
            A_w = np.dot(W, A)
            # Compute solution to the matrix equation W*Ax = W*b:
            x_data = np.dot(
                np.dot(inv(np.dot(A_w.conj().T, A_w)), A_w.conj().T), data_w.conj().T
            )
            # Compute weighted least square fit model:
            hfit = x_data[0]

        elif parameters == 2:
            # place parameters in model matrix:
            A[:, 0] = 1.0
            A[:, 1] = time_n
            # Weight model matrix
            A_w = np.dot(W, A)
            # Compute solution to the matrix equation W*Ax = W*b:
            x_data = np.dot(
                np.dot(inv(np.dot(A_w.conj().T, A_w)), A_w.conj().T), data_w.conj().T
            )
            # Compute weighted least square fit model:
            hfit = x_data[0] + x_data[1] * time

    elif trend == "sinusoidal":

        if parameters == 2:
            # place parameters in model matrix:
            A[:, 0] = 1.0
            A[:, 1] = np.cos((time_n * 2 * np.pi / period) - phase)
            # Weight model matrix
            A_w = np.dot(W, A)
            # Compute solution to the matrix equation W*Ax = W*b:
            x_data = np.dot(
                np.dot(inv(np.dot(A_w.conj().T, A_w)), A_w.conj().T), data_w.conj().T
            )
            # Compute weighted least square fit model:
            hfit = x_data[0] + x_data[1] * np.cos((time * 2 * np.pi / period) - phase)

        if parameters == 3:
            # place parameters in model matrix:
            A[:, 0] = 1.0
            A[:, 1] = np.sin(time_n * 2 * np.pi / period)
            A[:, 2] = np.cos(time_n * 2 * np.pi / period)
            # Weight model matrix
            A_w = np.dot(W, A)
            # Compute solution to the matrix equation W*Ax = W*b:
            x_data = np.dot(
                np.dot(inv(np.dot(A_w.conj().T, A_w)), A_w.conj().T), data_w.conj().T
            )
            # Compute weighted least square fit model:
            hfit = (
                x_data[0]
                + x_data[1] * np.sin(time * 2 * np.pi / period)
                + x_data[2] * np.cos(time * 2 * np.pi / period)
            )

        elif parameters == 4:
            # place parameters in model matrix:
            A[:, 0] = 1.0
            A[:, 1] = time_n
            A[:, 2] = np.sin(time_n * 2 * np.pi / period)
            A[:, 3] = np.cos(time_n * 2 * np.pi / period)
            # Weight model matrix
            A_w = np.dot(W, A)
            # Compute solution to the matrix equation W*Ax = W*b:
            x_data = np.dot(
                np.dot(inv(np.dot(A_w.conj().T, A_w)), A_w.conj().T), data_w.conj().T
            )
            # Compute weighted least square fit model:
            hfit = (
                x_data[0]
                + x_data[1] * time
                + x_data[2] * np.sin(time * 2 * np.pi / period)
                + x_data[3] * np.cos(time * 2 * np.pi / period)
            )

        elif parameters == 5:
            # place parameters in model matrix:
            A[:, 0] = 1.0
            A[:, 1] = np.sin(time_n * 2 * np.pi / period)
            A[:, 2] = np.cos(time_n * 2 * np.pi / period)
            A[:, 3] = np.sin(time_n * 2 * np.pi / (period / 2))
            A[:, 4] = np.cos(time_n * 2 * np.pi / (period / 2))
            # Weight model matrix
            A_w = np.dot(W, A)
            # Compute solution to the matrix equation W*Ax = W*b:
            x_data = np.dot(
                np.dot(inv(np.dot(A_w.conj().T, A_w)), A_w.conj().T), data_w.conj().T
            )
            # Compute weighted least square fit model:
            hfit = (
                x_data[0]
                + x_data[1] * np.sin(time * 2 * np.pi / period)
                + x_data[2] * np.cos(time * 2 * np.pi / period)
                + x_data[3] * np.sin(time * 2 * np.pi / (period / 2))
                + x_data[4] * np.cos(time * 2 * np.pi / (period / 2))
            )

        elif parameters == 6:
            # place parameters in model matrix:
            A[:, 0] = 1.0
            A[:, 1] = time_n
            A[:, 2] = np.sin(time_n * 2 * np.pi / period)
            A[:, 3] = np.cos(time_n * 2 * np.pi / period)
            A[:, 4] = np.sin(time_n * 2 * np.pi / (period / 2))
            A[:, 5] = np.cos(time_n * 2 * np.pi / (period / 2))
            # Weight model matrix
            A_w = np.dot(W, A)
            # Compute solution to the matrix equation W*Ax = W*b:
            x_data = np.dot(
                np.dot(inv(np.dot(A_w.conj().T, A_w)), A_w.conj().T), data_w.conj().T
            )
            # Compute weighted least square fit model:
            hfit = (
                x_data[0]
                + x_data[1] * time
                + x_data[2] * np.sin(time * 2 * np.pi / period)
                + x_data[3] * np.cos(time * 2 * np.pi / period)
                + x_data[4] * np.sin(time * 2 * np.pi / (period / 2))
                + x_data[5] * np.cos(time * 2 * np.pi / (period / 2))
            )

    # compute covariance matrix
    C = inv(np.dot(A_w.conj().T, A_w))

    # initialize x_data uncertainty vector
    x_data_sigma = np.zeros((N, 1))

    # Assign uncertainty of parameters x_data
    for isigma in range(N):
        x_data_sigma[isigma] = np.sqrt(C[isigma, isigma])

    return hfit, x_data, x_data_sigma


########### Least squares fit parameters Function ###########
def LSF_parameters(
    data, model, x_solution, trend, parameters, lsf="unweighted", sigma=False
):

    """
    LSF_parameters(data, model, x_solution, trend, parameters, lsf = 'unweighted', sigma = False)

        Function to compute the unweighted or weighted least square fit parameters

        Parameters
        ----------
        data : numpy 2D array of temporal data
            ex: print(data.shape) => (133,)
        model : Least square fit model computed using temporal data
        x_solution : x solutions from the matrix equation Ax = b
        trend : type of function least square is fitting. Options for this variables include trend == 'sinusoidal'
        parameters : number of parameters in model
               a) if trend = 'linear', then options for parameters includes:
                    i) parameters = 1 : fits mean
                    ii) parameters = 2 : fits mean and linear trend
               b) if trend = 'exponential', then the only option for paramters includes:
                    i) parameters = 2: fits linear trend and mean offset to a linearized exponential equation
               c) if trend = 'sinusiodal', then options for paramters includes:
                    i) parameters = 3 : fits mean, and annual cycle (period = 1 year, 12 months, or 365.25 days)
                    ii) parameters = 4 : fits mean, linear trend, and annual cycle
                    iii) parameters = 5 : fits mean, annual cycle, and semi-annual cycle (period= 1/2 year,6 months,182.625 days)
                    iv) paramters = 6: fits mean, linear trend, annual cycle, and semi-annual cycle
        lsf : Specifies whether the model is weight or unweigthed. By default,
              lsf = 'unweighted'. If the model is weighted, data is weighted by the
              uncertainty in order for the diagnostics to be consistent with the
              minimization procedure. Options include: lsf = 'unweighted' or 'weighted'.
        sigma : Uncertainty of temporal data. By default, sigma = False. Mask on sigma is
                assumed to agree with data.

        Returns
        -------
        residual : Residuals (data - model)
        rms : Root mean square error of least square fit (rough estimate of how well the
              model fits the data)
        amp1 : Amplitude of the Annual Cycle
        phase1 : Phase Constant of the Annual Cycle
        amp2 : Amplitude of the Semi-annual Cycle
        phase2 : Phase Constant of the Semi-annual Cycle
        fve : Fraction of Variance Explained (diagonstic the goodness of fit of the model)

        Libraries necessary to run function
        -----------------------------------
        Numpy : import numpy as np
    """

    # import library:
    import numpy as np

    # consider only unmasked data:
    data_n = data[~data.mask]
    model_n = model[~data.mask]
    if np.any(sigma):
        sigma_n = sigma[~data.mask]

    # Case 1: weighted least squares fit
    if lsf == "weighted":

        # Create weighted matrix
        W = np.ma.diag(1 / sigma_n)
        # Weight data
        data_w = np.ma.dot(W, data_n)
        # Weight model
        model_w = np.ma.dot(W, model_n)

    # Compute residuals, root mean squared, and coefficient of determination
    if lsf == "unweighted":
        # residual:
        residual = data_n - model_n
        # root mean square error:
        rms = np.sqrt(np.ma.mean((residual) ** 2))
        # coefficient of determination
        fve = 1 - (
            np.ma.sum((residual) ** 2) / np.ma.sum((data_n - np.ma.mean(data_n)) ** 2)
        )
    elif lsf == "weighted":
        # residual:
        residual = data_w - model_w
        # root mean square error:
        rms = np.sqrt(np.ma.mean((data_n - model_n) ** 2))
        # coefficient of determination
        fve = 1 - (
            np.ma.sum((residual) ** 2) / np.ma.sum((data_w - np.ma.mean(data_w)) ** 2)
        )

    # set conditional statements for each model case:
    if trend == "linear" or "exponential":
        # Annual cycle amplitude:
        amp1 = None
        # Annual cycle phase constant:
        phase1 = None
        # semi-annual cycle amplitude:
        amp2 = None
        # semi-annual cycle phase constant:
        phase2 = None

    if trend == "sinusoidal":

        # set the two cases for the linear trend parameters:
        if parameters == 3:
            # calcluate parameters from model:
            # annual cycle amplitude:
            amp1 = np.sqrt((x_solution[1] ** 2) + (x_solution[2] ** 2))
            # annual cycle phase constant:
            phase1 = np.arctan2(x_solution[1], x_solution[2])
            # semi-annual cycle amplitude:
            amp2 = None
            # semi-annual cycle phase constant:
            phase2 = None

        elif parameters == 4:
            # calcluate parameters from model:
            # annual cycle amplitude:
            amp1 = np.sqrt((x_solution[2] ** 2) + (x_solution[3] ** 2))
            # annual cycle phase constant:
            phase1 = np.arctan2(x_solution[2], x_solution[3])
            # semi-annual cycle amplitude:
            amp2 = None
            # semi-annual cycle phase constant:
            phase2 = None

        elif parameters == 5:
            # calcluate parameters from model:
            # annual cycle amplitude:
            amp1 = np.sqrt((x_solution[1] ** 2) + (x_solution[2] ** 2))
            # annual cycle phase constant:
            phase1 = np.arctan2(x_solution[1], x_solution[2])
            # semi-annual cycle amplitude:
            amp2 = np.sqrt((x_solution[3] ** 2) + (x_solution[4] ** 2))
            # semi-annual cycle phase constant:
            phase2 = np.arctan2(x_solution[3], x_solution[4])

        elif parameters == 6:
            # calcluate parameters from model:
            # Annual cycle amplitude:
            amp1 = np.sqrt((x_solution[2] ** 2) + (x_solution[3] ** 2))
            # Annual Cycle phase constant:
            phase1 = np.arctan2(x_solution[2], x_solution[3])
            # semi-annual cycle amplitude:
            amp2 = np.sqrt((x_solution[4] ** 2) + (x_solution[5] ** 2))
            # semi-annual cycle phase constant:
            phase2 = np.arctan2(x_solution[4], x_solution[5])

    return residual, rms, amp1, phase1, amp2, phase2, fve


########### Uncertainties of Phase and Amplitude Function ###########
def uncertainty_phase_amp(parameters, parameters_sigma, cycles, trend=None):

    """
    uncertainty_phase_amp(parameters, parameters_sigma)

        Function to propagate error through phase and amplitude calculations.

        Parameters
        ----------
        parameters: Parameters from Weighted least-squares fit
            ex: parameters = [a_0, a_1, a_2].T = x_data
        parameters_sigma : Uncertainties of parameters from Weighted least-squares fit
            ex: parameters = np.sqrt([C_11, C_22, C_33]) = x_data_sigma
        cycles : Number of cycles accounted for in the model. Options for this keyword argument:
            a) cycles = 'one' (Acounts for the annual cycle)
            b) cycles = 'two' (Accounts for the annual and semi-annual cycle)
        trend : Specifies whether the model accounts for a trend in the data. Options include:
            a) trend = None (Default)
            b) trend = 'linear'

        Returns
        -------
        sigma_phase_1: Uncertainty in Annual cycle Phase
        sigma_amp_1: Uncertainty in Annual cycle Amplitude
        sigma_phase_2: Uncertainty in Semi-Annual cycle Phase
        sigma_amp_2: ncertainty in Semi-Annual cycle Amplitude

        Libraries necessary to run function
        -----------------------------------
        Numpy : import numpy as np
    """

    # import libraries:
    import numpy as np

    # Using propagation of error (assuming uncertainties are uncorrelated), compute uncertainty of phase and amp:

    # case 1: Only Annual cycle
    if cycles == "one":

        # case 1: No trend
        if trend == None:
            # assign parameter variables:
            a_1 = parameters[1]
            a_2 = parameters[2]
            sig_1 = parameters_sigma[1]
            sig_2 = parameters_sigma[2]

        else:
            # assign parameter variables:
            a_1 = parameters[2]
            a_2 = parameters[3]
            sig_1 = parameters_sigma[2]
            sig_2 = parameters_sigma[3]

        # Compute uncertainties
        sigma_amp_1 = (1 / np.sqrt(a_1 ** 2 + a_2 ** 2)) * np.sqrt(
            (sig_1 * a_1) ** 2 + (sig_2 * a_2) ** 2
        )
        sigma_phase_1 = np.sqrt(
            sig_1 ** (2) * (-a_2 / (a_1 ** (2) + a_2 ** (2))) ** (2)
            + sig_2 ** (2) * (a_1 / (a_1 ** (2) + a_2 ** (2))) ** (2)
        )
        # cycle 2:
        sigma_phase_2 = None
        sigma_amp_2 = None

    # case 2: Annual and Semi-annual Cycles
    if cycles == "two":

        # case 1: No trend
        if trend == None:
            # assign parameter variables:
            a_1 = parameters[1]
            a_2 = parameters[2]
            a_3 = parameters[3]
            a_4 = parameters[4]
            sig_1 = parameters_sigma[1]
            sig_2 = parameters_sigma[2]
            sig_3 = parameters_sigma[3]
            sig_4 = parameters_sigma[4]

        else:
            # assign parameter variables:
            a_1 = parameters[2]
            a_2 = parameters[3]
            a_3 = parameters[4]
            a_4 = parameters[5]
            sig_1 = parameters_sigma[2]
            sig_2 = parameters_sigma[3]
            sig_3 = parameters_sigma[4]
            sig_4 = parameters_sigma[5]

        # Compute uncertainties
        # cycle 1:
        sigma_amp_1 = (1 / np.sqrt(a_1 ** 2 + a_2 ** 2)) * np.sqrt(
            (sig_1 * a_1) ** 2 + (sig_2 * a_2) ** 2
        )
        sigma_phase_1 = np.sqrt(
            sig_1 ** (2) * (-a_2 / (a_1 ** (2) + a_2 ** (2))) ** (2)
            + sig_2 ** (2) * (a_1 / (a_1 ** (2) + a_2 ** (2))) ** (2)
        )
        # cycle 2:
        sigma_amp_2 = (1 / np.sqrt(a_3 ** 2 + a_4 ** 2)) * np.sqrt(
            (sig_3 * a_3) ** 2 + (sig_4 * a_4) ** 2
        )
        sigma_phase_2 = np.sqrt(
            sig_3 ** (2) * (-a_4 / (a_3 ** (2) + a_4 ** (2))) ** (2)
            + sig_4 ** (2) * (a_3 / (a_3 ** (2) + a_4 ** (2))) ** (2)
        )

    return sigma_phase_1, sigma_amp_1, sigma_phase_2, sigma_amp_2
