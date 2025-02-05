#!/usr/bin/env python
"""LAMMPS snippets for molecular dynamics simulations."""

from .plumed_snippets import post_us

# Pre-processing snippets
def pre_us():
    print("""
# see example in enhanced_sampling\\abc(rectangle_potential)
# umbrella sampling totally using LAMMPS
dump mydump all atom 10000 dump_us.lammpstrj

variable k equal 1.5 # spring constant, Kcal/mol/Angstrom^2
variable a loop 51
label loop
variable xdes equal ${a}-26
variable xave equal xcm(topull,x)
fix mytth topull spring tether ${k} ${xdes} NULL NULL 0
run 200000
fix myat1 all ave/time 1 1 20 v_xave v_xdes file position.${xdes}.dat
run 1000000 # data acquisition
write_data ${xdes}.data

unfix myat1
next a
jump SELF loop
""")

def pre_msd():
    print("""
## ------------------ MSD -------------------
# dump unwrapped coords and post-processing
# you can specify the group of interest, like water, to decrease the size of trajectory file
# <DUMP_FREQ> : 20 , you can choose larger DUMP_FREQ to reduce the file size
dump mydump all custom <DUMP_FREQ> dump.lammpstrj id type x y z ix iy iz #vx vy vz fx fy fz
dump_modify mydump sort id
""")
    
    
def pre_vacf():
    pass
    
def pre_fric_coeff():
    print("""
# ------------- friction coefficient -------------
# freeze CNT
fix freeze cnt setforce 0.0 0.0 0.0
fix friction all ave/time 1 1 1 f_freeze[*] file friction.txt   
""")

    
    
# density profile
def pre_1d_profile():
    print("""
## ------ 1D density profile -------
compute bin_x all chunk/atom bin/1d x 0.0 0.2 units box
fix density_bin_x all ave/chunk 10 100000 1000000 bin_x density/number density/mass temp file density_x.profile
""")

def pre_radial_profile():
    print(""" 
# ------------- radial profile ------------
region CNT_center_slab_y block -0.2 0.2 -${radius_CNT} ${radius_CNT}  INF INF units box
compute bin_slab_y_water water chunk/atom bin/1d y 0.0 0.1 region CNT_center_slab_y
fix slab_y_water water ave/chunk 10 10000 100000 bin_slab_y_water density/number vx vy vz file slab_y_water.profile

compute bin_slab_y_oxygen oxygen chunk/atom bin/1d y 0.0 0.1 region CNT_center_slab_y
fix slab_y_oxygen oxygen ave/chunk 10 10000 100000 bin_slab_y_oxygen density/number vx vy vz file slab_y_oxygen.profile
""")

def pre_viscosity():
    print("""
# ------------ stress tensor for shear viscosity -------------
compute atom_stress water  stress/atom T_solution
compute sol_stress water reduce sum c_atom_stress[*]
fix mystress  all ave/time 1 1 1 c_sol_stress[*] file stress.time     

# ------------- bulk stress tensor for shear viscosity -------------
variable pxy equal pxy
variable pxz equal pxz
variable pyz equal pyz
variable pxx equal pxx
variable pyy equal pyy
variable pzz equal pzz
fix bulk_stress all ave/time 1 1 1 v_pxy v_pxz v_pyz v_pxx v_pyy v_pzz file bulk_stress.time
""")

def pre_rdf():
    print("""
## --------------- RDF -----------------
compute atom_rdf all rdf 500 1 1 1 2 2 2 
fix 2 all ave/time 10 100000 1000000 c_atom_rdf[*] file atom.rdf mode vector

# usually, we do not need to compute RDF for all atoms
# compute all_rdf all rdf 500
# fix 1 all ave/time 10 100000 1000000 c_all_rdf[*] file all.rdf mode vector
""")

def pre_dump():
    print("""
## ----- MSD, hydrogen bond, water dipole distribution -----
## dump unwrapped coords and post-processing
dump mydump solution custom 20 dump_sol.lammpstrj id type x y z ix iy iz vx vy vz fx fy fz
dump_modify mydump sort id

# ---------- trajecotry --------------
dump myDump2 all xyz 1000 dump.xyz
dump_modify myDump2 element O H C
""")

# Post-processing snippets
def post_msd():
    print("""
## ------------------ MSD -------------------
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda

# read trajectory file
dt = <YOUR_DT> # timestep, in ps. usually 0.001 or 0.002
dump_freq = <DUMP_FREQ> # the same as pre-processing
dump_dt = dt * dump_freq
u = mda.Universe("system.data", "dump.lammpstrj",topology_format="DATA", format="LAMMPSDUMP",dt=dump_dt, lammps_coordinate_convention="auto", unwrap_images=True)

# run MSD analysis
mymsd = msd.EinsteinMSD(u, select=<SELECTION>, msd_type="z",fft=True)
step = <MSD_STEP> # the step of EinsteinMSD to process the trajectory file
mymsd.run(verbose=True,step=step)

time = np.arange(mymsd.n_frames)*dump_dt*step   # in ps
msd_data = mymsd.results.timeseries             # in Ang^2

# plot MSD curve
fig,ax = plt.subplots()
ax.plot(time,msd_data, c="black", ls="-")
ax.set_xlabel("$t$ (ps)")
ax.set_ylabel("MSD ($\mathrm{\AA^2}$)")
""")

# Command handlers
def handle_pre(option):
    """Handle pre-processing snippets"""
    handlers = {
        "msd": pre_msd,
        "fric_coeff": pre_fric_coeff,
        "1d_profile": pre_1d_profile,
        "radial_profile": pre_radial_profile,
        "viscosity": pre_viscosity,
        "rdf": pre_rdf,
        "dump": pre_dump,
        "us": pre_us
    }
    
    if option.lower() == "all":
        for func in handlers.values():
            func()
            print("\n" + "="*50 + "\n")
    else:
        handler = handlers.get(option.lower())
        if handler:
            handler()
        else:
            print(f"Unknown option: {option}")
            print(f"Available options: {', '.join(handlers.keys())}")

def handle_post(option):
    """Handle post-processing snippets"""
    handlers = {
        "msd": post_msd,
        "us": post_us
    }
    
    if option.lower() == "all":
        for func in handlers.values():
            func()
            print("\n" + "="*50 + "\n")
    else:
        handler = handlers.get(option.lower())
        if handler:
            handler()
        else:
            print(f"Unknown option: {option}")
            print(f"Available options: {', '.join(handlers.keys())}")