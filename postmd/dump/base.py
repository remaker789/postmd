'''
Trajectory Analysis Module
=========================
Read and analyze molecular dynamics trajectories based on MDAnalysis package.

Supported Formats
----------------
xyz : XYZ trajectory format
lammpstrj : LAMMPS trajectory dump format
'''
import numpy as np
from .io import read_xyz, read_lammpstrj
import warnings

SUPPORTED_FORMATS = ['xyz', 'lammpstrj']


class Traj:
    '''
    Handle molecular dynamics trajectories.

    Parameters
    ----------
    topology_file : str
        Topology file path
    traj_file : str
        Trajectory file path
    traj_format : {'xyz', 'lammpstrj'}
        Trajectory format
    dt : float, optional
        Timestep (ps)
    traj_dt : float, optional
        Time interval between frames (ps)
    lammps_coordinate_convention : optional
        LAMMPS coordinate convention, default='auto'
    unwrap_images : bool, optional
        Whether to unwrap periodic boundary conditions. Default is False.
    **kwargs : dict
        Additional arguments for trajectory reader
    '''
    


    def __init__(self, topology_file=None, traj_file=None, traj_format=None, dt=None, 
                 traj_dt=None, lammps_coordinate_convention="auto", unwrap_images=False, **kwargs):
        self.topology_file = topology_file
        self.traj_file = traj_file
        self.dt = dt
        self.traj_dt = traj_dt
        self.traj_format = self._check_format(traj_format)
        self.lammps_coordinate_convention = lammps_coordinate_convention
        self.unwrap_images = unwrap_images
        self.kwargs = kwargs  # 存储额外的参数
        self.traj = self.read_traj()
        
    def _check_format(self, traj_format):
        '''
        Check if trajectory format is supported.

        Parameters
        ----------
        traj_format : str
            Trajectory format to check

        Returns
        -------
        str
            Validated and normalized format string

        Raises
        ------
        ValueError
            If format is None or not supported
        '''
        if traj_format is None:
            raise ValueError(f"Trajectory format required. Supported: {SUPPORTED_FORMATS}")
            
        traj_format = traj_format.lower()
        if traj_format not in SUPPORTED_FORMATS:
            raise ValueError(f"Unsupported format: {traj_format}. Supported: {SUPPORTED_FORMATS}")
            
        return traj_format

    def read_traj(self):
        '''
        Read trajectory file.

        Returns
        -------
        MDAnalysis.Universe
        '''
        if self.topology_file is None:
            warnings.warn("No topology file provided. Some functions may be unavailable. We advise to use data file to import molecule information if possible.")
            
        if self.traj_format == "xyz":
            traj = read_xyz(self.topology_file, self.traj_file, self.dt, **self.kwargs)
        elif self.traj_format == "lammpstrj":
            traj = read_lammpstrj(self.topology_file, self.traj_file, self.dt, self.lammps_coordinate_convention, self.unwrap_images, **self.kwargs)
        
        return traj
    
    def extract_atomic_vel(self, outfile="atomic_vel.csv"):
        '''
        Extract atomic velocities to CSV.

        Parameters
        ----------
        outfile : str, optional
            Output file path, default="atomic_vel.csv"

        Notes
        -----
        CSV format: time(ps),v1x,v1y,v1z,v2x,v2y,v2z,...
        '''
        try:
            with open(outfile, "w") as f:
                f.write("# time(ps),v1x,v1y,v1z,v2x,v2y,v2z,...\n")
                
                for frame in self.traj.trajectory:
                    t = frame.frame * self.traj_dt  # time in ps
                    if self.traj_format == "lammpstrj":
                        data = frame.velocities
                    elif self.traj_format == "xyz":
                        data = frame.positions
                        
                    # Flatten data and combine with time
                    row = np.concatenate(([t], data.flatten()))
                    np.savetxt(f, row.reshape(1, -1), delimiter=",")
                    
        except (IOError, OSError) as e:
            raise IOError(f"Error writing to file {outfile}: {str(e)}")
        except Exception as e:
            raise RuntimeError(f"Error extracting atomic data: {str(e)}")

    
