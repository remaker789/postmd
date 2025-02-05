def extract_atomic_vel(Traj):
    import MDAnalysis as mda
    u = mda.Universe("../system.data", "../msd_wrapped_coords.lammpstrj",topology_format="DATA", format="LAMMPSDUMP",dt=0.2, lammps_coordinate_convention="auto", unwrap_images=False)
    timestep = 0.002 # in ps
    steps_traj = 100 # trajectory per how many steps

    filename = "atomic_vel.csv"
    with open(filename, "w+") as f:
        f.write("# time,v1x,v1y,v1z,v2x,v2y,v2z,...")

    with open(filename, "a") as f:
        for frame in u.trajectory[:2001]: # 前400 ps的速度
            # print(frame)
            # print(type(frame))
            t = frame.frame*timestep*steps_traj # in ps
            np.savetxt(f, np.concatenate(([t],frame.velocities.flatten())).reshape(1,-1),delimiter=",")
    
    # then read the data
    # atomic_vel_data = np.loadtxt("atomic_vel.csv",delimiter=",")