import numpy as np
import matplotlib.pyplot as plt

from fastmsd import mean_squared_displacement as fast_msd
from fastrdf import radial_distribution_function as fast_rdf


class Trajectory:
    """Molecular dynamics trajectory.

    Read and parse a MD trajectory file, generate an animated plot of the simulation and compute **Mean Squared Displacement** as well as the **Radial Distribution Function**.

    Parameters
    ----------
    path : str
        Path to the `.xyz` file containing the trajectory data.
    verbose : bool, default=False

    Notes
    -----
    .. important::
        Most instance attributes are **protected** (names with a single leading underscore).
        Those should not be modified outside of class methods.
    """


    def __init__(self, path:str, verbose=False):

        self._path = path
        self.verbose = verbose

        self._n_frames, self._n_atoms, self._box_size, self._timestamps, self._coordinates, self._elements = Trajectory.read_xyz(self._path)

        self._r_cut = self._box_size.min() * 0.5 # define outer cutoff distance as half of the smallest side of the box
        self._n_bins = 200 # number of bins to split the distance interval


    def __repr__(self):
        
        unique_elements = str.join(', ', np.unique(self._elements))
        string = 'Trajectory information:\n'
        string += f'File path: "{self._path}"\n'
        string += f'Number of atoms: {self._n_atoms}, number of frames: {self._n_frames}\n'
        string += f'Box size: {self._box_size}, element(s): {unique_elements}'

        return string
    

    @staticmethod
    def read_xyz(path:str):
        """Read the `.xyz` file at a given path.

        Parameters
        ----------
        path : str
            Path of the file.

        Returns
        -------
        n_frames : int
            Number of frames.
        n_atoms : int
            Number of atoms per frame, assumed to be constant throughout the simulation.
        box_size : ndarray
            1D array containing the dimensions of the simulation box, assumed to always be cubic.
        timestamps : ndarray
            1D array of the timestamps.
        coordinates : ndarray
            3D array of the `(x, y, z)` coordinates of each atom in each frame.
        elements : ndarray
            1D array of the type of each atom.
        """

        with open(path) as f:

            lines = f.readlines()

        n_atoms = int(lines[0])
        lskip = n_atoms + 2
        n_frames = len(lines) // lskip
        box_size = np.array(lines[1].split()[9:12], dtype=float)

        timestamps = np.empty(n_frames, dtype=float)
        coordinates = np.empty((n_frames, n_atoms, 3), dtype=float)
        elements = np.array(str.join(' ', lines[2:lskip]).split()[::4], dtype='<U2')

        for frame in range(n_frames):

            timestamps[frame] = lines[frame*lskip + 1].split()[4]
            coordinates[frame] = np.array(str.join(' ', lines[frame*lskip+2:(frame+1)*lskip]).split()).reshape(-1, 4)[:, 1:].astype(float)

        return n_frames, n_atoms, box_size, timestamps, coordinates, elements
    

    def mean_squared_displacement(self):

        msd = fast_msd(self._timestamps, self._coordinates)


    def radial_distribution_function(self):

        rdf = fast_rdf(self._r_cut, self._n_bins, self._box_size, self._coordinates)
