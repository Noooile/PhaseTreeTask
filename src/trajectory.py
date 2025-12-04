import numpy as np


class Trajectory:
    """Molecular dynamics trajectory.

    Read and parse a MD trajectory file, generate an animated plot of the simulation and compute **Mean Squared Displacement** as well as the **Radial Distribution Function**.

    Parameters
    ----------
    path : str
        Path to the `.xyz` file containing the trajectory data.
    verbose : bool, default=False

    .. important::
        Most instance attributes are **protected** (names with a single leading underscore).
        Those should not be modified outside of class methods.
    """


    def __init__(self, path:str, verbose=False):

        self._path = path
        self.verbose = verbose

        self._n_config, self._n_atoms, self._box_size, self._timestamps, self._coordinates, self._elements = Trajectory.read_xyz(self._path)


    def __repr__(self):
        
        unique_elements = str.join(', ', np.unique(self._elements))
        string = 'Trajectory information:\n'
        string += f'File path: "{self._path}"\n'
        string += f'Number of atoms: {self._n_atoms}, number of configurations: {self._n_config}\n'
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
        n_config : int
            Number of configurations (or frames).
        n_atoms : int
            Number of atoms per configurations, assumed to be constant throughout the simulation.
        box_size : ndarray
            1D array containing the dimensions of the simulation box, assumed to always be cubic.
        timestamps : ndarray
            1D array of the timestamps.
        coordinates : ndarray
            3D array of the `(x, y, z)` coordinates of each atom in each configuration.
        elements : ndarray
            1D array of the type of each atom.
        """

        with open(path) as f:

            lines = f.readlines()

        n_atoms = int(lines[0])
        lskip = n_atoms + 2
        n_config = len(lines) // lskip
        box_size = np.array(lines[1].split()[9:12], dtype=float)

        timestamps = np.empty(n_config, dtype=float)
        coordinates = np.empty((n_config, n_atoms, 3), dtype=float)
        elements = np.array(str.join(' ', lines[2:lskip]).split()[::4], dtype='<U2')

        for config in range(n_config):

            timestamps[config] = lines[config*lskip + 1].split()[4]
            coordinates[config] = np.array(str.join(' ', lines[config*lskip+2:(config+1)*lskip]).split()).reshape(-1, 4)[:, 1:].astype(float)

        return n_config, n_atoms, box_size, timestamps, coordinates, elements
