import numpy as np
from scipy.linalg import pinv
from scipy.linalg import convolution_matrix
from scipy.interpolate import interp1d
from scipy.interpolate import PchipInterpolator


def deconv_gp(x, Y, TR, tseg, ntaps, k, dt_rs, ll=4, sigf=1, interp='pchip'):
    """
    Modified from MATLAB implementation:
    for deconvolution of impulse response functions using Gaussian Process priors.
    based on paper: 
    Chang et al., Neuroimage 44(3), 857-869
    
    this version 20 Apr 2019 - catie

    Hyperparameters of GPP Regression:
    ll = 4
    sigf = 1

    """
    # time axis
    tax_h = np.arange(-k, ntaps - k) * TR

    # voxelwise fmri data
    YY = Y[tseg, :]
    YYC = YY - np.mean(YY, axis=0)

    # shift fmri forward by "k" points:
    YS = np.vstack((np.zeros((k, YYC.shape[1])), YYC[:-k]))

    # append 0's to match convmtx dim
    YD = np.vstack((YS, np.zeros((ntaps - 1, YYC.shape[1]))))

    # global signal (& use for scaling below)
    glob = np.mean(Y, axis=1)
    y0 = glob[tseg] - np.mean(glob[tseg])
    ys = np.concatenate((np.zeros(k), y0[:-k]))
    yd = np.vstack((ys[:,np.newaxis], np.zeros((ntaps - 1, 1))))
    yd = (yd - np.mean(yd)) / np.std(yd)

    # normalize time series  to unit variance
    YD = (YD - np.mean(YD, axis=0)) / np.std(YD, axis=0)
    
    # set up input matrix
    X = convolution_matrix(x[tseg][:, 0], ntaps, mode='full')

    # build cov matrxs
    # hyperparameters (fixed here, but may also be estimated from the data)
    sigv = np.std(yd)  # will be fixed for all voxels
    nu = sigf ** 2
    lsq = ll ** 2
    sigvsq = sigv ** 2
    # build covariance matrix
    Kp = np.zeros((ntaps, ntaps))
    for i in range(ntaps):
        for j in range(ntaps):
            Kp[i, j] = nu * np.exp(-0.5 * (1 / lsq) * (i - j) ** 2)
    K = Kp + 1e-5 * np.eye(ntaps)  # add small values to diagonal
    R = np.linalg.inv(K)
    # set up boundary conditions (both impulse responses begin and ends @ 0)
    b = ntaps
    S = np.zeros((2, b))
    S[0, 0] = 1
    S[1, -1] = 1

    # deconv voxelwise
    P = 2 * ((1 / sigvsq) * np.dot(X.T, X) + R)
    q = (-2 / sigvsq) * np.dot(X.T, YD)
    A = S
    bb = np.zeros((2, YD.shape[1]))
    M = np.block([[P, A.T], [A, np.zeros((2, 2))]])
    m = np.block([[-q], [bb]])
    xxstar = np.linalg.solve(M, m)
    H = xxstar[:ntaps, :]

    # also OLS
    B = pinv(X).dot(YD)

    # resample outputs to desired resolution
    tax_in = tax_h
    tax_rs = np.arange(tax_h[0], tax_h[-1], dt_rs)
    if interp == 'linear':
        H_rs = interp1d(tax_in, H, axis=0, kind='linear', fill_value='extrapolate')(tax_rs)
        B_rs = interp1d(tax_in, B, axis=0, kind='linear', fill_value='extrapolate')(tax_rs)
    elif interp == 'pchip':
        H_rs = PchipInterpolator(tax_in, H, axis=0)(tax_rs)
        B_rs = PchipInterpolator(tax_in, B, axis=0)(tax_rs)

    # return outputs (& params)
    OUT = {'H': H, 'B': B, 'H_rs': H_rs, 'B_rs': B_rs,
           'll': ll, 'tseg': tseg, 'k': k, 'x': x, 'dt_rs': dt_rs, 'tax_rs': tax_rs}
    return OUT
