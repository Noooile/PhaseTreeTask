import numpy as np
from numba import njit, prange



@njit(parallel=True)
def radial_distribution_function(r_cut:float, n_bins:int, box_size:np.ndarray, coordinates:np.ndarray):
    """Compute the **Radial Distribution Function** of a trajectory averaged over time.

    The **RDF** is time-averaged and ensemble-averaged: the computation uses the whole trajectory as well as all available atoms.

    Parameters
    ----------
    r_cut : float
        Outer cutoff distance.
    n_bins : int
        Number of bins used to evenly split the distance interval.
    box_size: ndarray
        1D array containing the dimensions of the simulation box, assumed to always be cubic.
    coordinates : ndarray
        3D array of the `(x, y, z)` coordinates of each atom in each frame.

    Returns
    -------
    distribution : ndarray
        1D array of the time and ensemble averaged **RDF**.

    Notes
    -----
    This function is compiled "just-in-time" using `numba`.
    """

    n_frames = coordinates.shape[0]
    n_atoms = coordinates.shape[1]

    rsq_cut = r_cut ** 2

    lattice = np.diag(box_size)
    lattice_inv = np.diag(1/box_size)

    density = n_atoms / np.prod(box_size)
    r_bins = np.linspace(0, r_cut, n_bins+1)
    volume_bins = 4 * np.pi * r_bins[1:]**2 * r_bins[1]

    counts = np.zeros(n_bins)

    # average over all frames
    for frame in prange(n_frames):

        # average over all atoms
        for i, pos_i in enumerate(coordinates[frame]):

            # find neighbors of atom i
            # iterate for indices j > i (taking advantage of r_ji = -r_ij)
            for j, pos_j in enumerate(coordinates[frame, (i+1):]):

                if j==i:
                    continue

                pos_ij = pos_j - pos_i
                pos_ij = pos_ij - lattice @ np.round(lattice_inv @ pos_ij) # minimum image convention for the PBC

                rsq_ij = np.dot(pos_ij, pos_ij) # squared norm

                if rsq_ij < rsq_cut:

                    counts[int(np.floor(np.sqrt(rsq_ij)/r_bins[1]))] += 2 # counts twice, one for r_ij and one for r_ji

    distribution = counts / (density * volume_bins * n_frames * n_atoms)

    return distribution
