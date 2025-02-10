import MDAnalysis as mda

def read_lammpstrj(topology_file,traj_file,dt,lammps_coordinate_convention="auto",unwrap_images:bool=False,**kwargs):
    u = mda.Universe(topology=topology_file, coordinates=traj_file,topology_format="DATA", format="LAMMPSDUMP",dt=dt, lammps_coordinate_convention=lammps_coordinate_convention, unwrap_images=unwrap_images, **kwargs)
    return u
    
def read_xyz(topology_file, traj_file, dt, **kwargs):
    u = mda.Universe(topology=topology_file, coordinates=traj_file,topology_format="DATA", format="XYZ",dt=dt, **kwargs)
    return u