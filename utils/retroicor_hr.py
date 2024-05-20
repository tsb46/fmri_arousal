import numpy as np
import scipy.signal as signal

def retroicor_hr(slice_order, TR, cardiac_trig_times, nframes):
    """
    Python implementation of cardiac component of retroicor_rvhr function from MATLAB 
    created by Catie Chang.

    Parameters:
    - slice_order: list indicating order of slice acquisition
    - TR: repetition time in seconds
    - cardiac_trig_times: list of cardiac (R-wave peak) times in seconds relative to the start of the fMRI scan
    - nframes: number of volumes of fmri scan

    Returns:
    - PHASES: list of cardiac phases for each slice
    - REGRESSORS: the cardiac retroicor regressors for each slice
    """
    # Process inputs
    nslc = len(slice_order)

    # Get the time "slot" at which each slice was acquired
    slice_slots = np.argsort(slice_order)

    # Cardiac input
    etrig = np.array(cardiac_trig_times)


    # Find cardiac and respiratory phase vectors
    PHASES_0 = []
    for jj in range(nslc):
        slice_times = np.arange((TR/nslc)*(jj+1-0.5), TR*nframes, TR)
        phases_thisSlice = []
        for st in slice_times:
            prev_trigs = np.where(etrig < st)[0]
            t1 = etrig[prev_trigs[-1]] if prev_trigs.size > 0 else 0
            next_trigs = np.where(etrig > st)[0]
            t2 = etrig[next_trigs[0]] if next_trigs.size > 0 else nframes * TR
            phi_cardiac = (st - t1) * 2 * np.pi / (t2 - t1)
            phases_thisSlice.append(phi_cardiac)
        
        PHASES_0.append(phases_thisSlice)


    PHASES = [PHASES_0[i] for i in slice_slots]

    # Generate slice-specific retroicor regressors
    REGRESSORS_RET = np.empty((nframes, 4, nslc))
    for jj in range(nslc):
        phi_c = np.array(PHASES[jj])
        covs = [np.cos(phi_c), np.sin(phi_c), np.cos(2 * phi_c), np.sin(2 * phi_c)]
        REGRESSORS_RET[:, :, jj] = np.column_stack(covs)

    # quadratic detrending
    REGRESSORS = np.empty((nframes, 4, nslc))
    for jj in range(nslc):
        REGRESSORS[:,:,jj] = quaddetrend_cols(REGRESSORS_RET[:,:,jj])

    return PHASES, REGRESSORS


def quaddetrend_cols(Y):
    """
    Remove a quadratic trend from each column of the data matrix.

    Parameters:
    - Y: 2D numpy array where each column is a time series

    Returns:
    - Yq: 2D numpy array with the quadratic trend removed from each column
    """
    x = np.arange(1, Y.shape[0] + 1)
    Yq = np.copy(Y)
    for j in range(Y.shape[1]):
        y = Y[:, j]
        p = np.polyfit(x, y, 2)
        ytrend = np.polyval(p, x)
        y_detrended = y - ytrend
        Yq[:, j] = y_detrended
    return Yq
