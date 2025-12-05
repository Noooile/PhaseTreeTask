import numpy as np
from numba import njit, prange



@njit(parallel=True)
def mean_squared_displacement(coordinates:np.ndarray):
    """Compute the **Mean Squared Displacement** of a trajectory over all possible time interval.

    The **MSD** is time-averaged and ensemble-averaged: the computation uses the whole trajectory as well as all available atoms.

    Parameters
    ----------
    coordinates : ndarray
        3D array of the `(x, y, z)` coordinates of each atom in each frame.

    Returns
    -------
    msd : ndarray
        1D array of the time and ensemble averaged **Mean Squared Displacement** for all possible time intervals.

    Notes
    -----
    This function is compiled "just-in-time" using `numba`.
    """

    n_frames = coordinates.shape[0]
    n_atoms = coordinates.shape[1]
    
    msd = np.zeros(n_frames, dtype=float)

    # loop over all possible time lag
    # assume a constant time-step
    for m in prange(n_frames):

        # loop over all atoms (ensemble-average)
        for i in range(n_atoms):

            # loop over all equal time intervals (time-average)
            for k in range(n_frames - m):

                msd[m] += (coordinates[k+m, i, 0] - coordinates[k, i, 0])**2
                msd[m] += (coordinates[k+m, i, 1] - coordinates[k, i, 1])**2
                msd[m] += (coordinates[k+m, i, 2] - coordinates[k, i, 2])**2

        msd[m] = msd[m] / (n_atoms * (n_frames - m))

    return msd
