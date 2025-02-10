#!/usr/bin/env python
"""PLUMED snippets for molecular dynamics simulations."""

def pre_us():
    print("""
# see example in enhanced_sampling\\abc(rectangle_potential)
# ------------------ us.plumed.template ------------------
a: FIXEDATOM AT=0,0,0
# DISTANCE: the coords of second atom minus the ones of first atom
d1: DISTANCE ATOMS=a,1 COMPONENTS
# make a new variable d1x, equal to d1.x but *periodic*
d1x: COMBINE ARG=d1.x PERIODIC=-2.5,2.5

# umbrella sampling needs to apply hamonic restriction in multiple windows(bins)
# we set {LT_AT} as -25, -24, -23, ...,0,..., 23, 24, 25 Angstrom, respectively.
myus: RESTRAINT ARG=d1x KAPPA=600.0 AT={LT_AT}

# output collective variable and bias to mycolvar.{LT_AT} file
PRINT ARG=d1x,myus.bias FILE=mycolvar.{LT_AT} STRIDE=20

# ------------------ python script --------------------
# Then, use following python scripts to generate multiple plumed scripts for umbrella sampling
import os
os.chdir(os.path.dirname(__file__))
import numpy as np

# substitute {LT_AT} in us.plumed.template
for i in np.arange(-2.5,2.6,0.1):
    with open(f"us.plumed.{i:.1f}", "w+", encoding='utf-8') as f:
        with open("us.plumed.template", "r", encoding='utf-8') as content:
            for line in content:
                f.write(line.replace(r"{LT_AT}", f"{i:.1f}"))
                
                
# ------------------ lammps input --------------------
# a loop of CVs
variable a loop 51
label loop

variable a1 equal (${a}-26)/10.0
variable a2 format a1 %.1f

# data collection, usually you need to exclude initially unstable data!
fix myus all plumed plumedfile us.plumed.${a2} outfile plumed.${a2}.log
dump mydump all xtc 100 dump.${a2}.xtc # you can choose other formats of trajectory
run 2000000
write_data system.${a2}.data

unfix myus
undump mydump
next a
jump SELF loop
""")






def pre_metad():
    print("""
# ------- Set collective variables --------
# Activate MOLINFO functionalities
# Compute the backbone dihedral angle phi, defined by atoms C-N-CA-C
# you should use MOLINFO shortcuts
phi: TORSION ATOMS=5,7,9,15
# Compute the backbone dihedral angle psi, defined by atoms N-CA-C-N
# here also you should to use MOLINFO shortcuts
psi: TORSION ATOMS=7,9,15,17

# -------- Activate *well-tempered metadynamics* in phi and psi ----------
metad: METAD ARG=phi,psi ...
# Deposit a Gaussian every 500 time steps, with initial height
# equal to 1.2 kJ/mol and bias factor equal to 8.
# *Set BIASFACTOR means turn on well-tempered metadynamics*
  PACE=500 HEIGHT=1.2 BIASFACTOR=8
# Gaussian width (sigma) should be chosen based on the CV fluctuations
# in unbiased run. The unit of SIGMA is the same as the collective variable.
  SIGMA=0.3,0.3
  TEMP=300.0
# Gaussians will be written to HILLS file and also stored on grid
  FILE=HILLS GRID_MIN=-pi,-pi GRID_MAX=pi,pi
...

# Print collective variables and bias to COLVAR file every 100 steps
PRINT ARG=phi,psi,metad.bias FILE=COLVAR STRIDE=100
""")    


# todo : 晚上umbrella sampling的后处理，应该有一个生成metadata.dat的脚本，然后还有一个plot_overlap的脚本
def post_us():
    print("""
# we use Alan Grossfield's WHAM program, 
# see http://membrane.urmc.rochester.edu/?page_id=126
# wham <CV_lo> <CV_hi> <num_window> <tor> <T> <numpad> metadata.dat PMF.dat [num_MC_trials randSeed]
wham -2.2 2.2 51 1e-8 300.0 0 metadata.dat PMF.dat
# PMF.dat: the output PMF of CV

# read the generated 'PMF.dat'
df = pd.read_csv("PMF.dat",sep=r"\s+",comment="#",header=None)


# plot PMF
fig, ax = plt.subplots()
ax.plot(df.iloc[:,0],df.iloc[:,1])
ax.set_xlabel(r"$z$ ($\mathrm{\AA}$)")
ax.set_ylabel("PMF (kcal/mol)")
# fig.savefig("m12_pmf.svg")

# ========= Overlap between neighboring windows =========

# read CV files
df = pd.read_csv("metadata.dat",sep=r"\s+", header=None)
cv_files = df.iloc[:,0]
spring_eqms = df.iloc[:,1]
spring_ks = df.iloc[:,2]
histograms = []
for cv_file in cv_files:
    cv_data = np.loadtxt(cv_file)
    cvs = cv_data[:,1]
    hist, bin_edges = np.histogram(cvs, bins = 50)
    prob = hist/hist.sum() # normalize
    bin_mid = (bin_edges[1:] + bin_edges[:-1]) / 2 # use the middle of bin as x axis
    histograms.append([bin_mid,prob])

# plot overlap
fig,ax = plt.subplots(figsize=(8,6))
for cpt, (bin_mid,prob) in enumerate(histograms):
    if cpt % 2 == 0:
        ax.plot(bin_mid, prob, color="C0", linewidth=1)
    else:
        ax.plot(bin_mid, prob, color="C1", linewidth=1)
        
        
# calculate the overlap ratio between neighboring windows
# If the overlap ratio is too small (<0.10), the PMF results will be suspicious.
# https://stackoverflow.com/questions/72931022/how-to-calculate-histogram-intersection
def calc_windows_overlap(cv_file1, cv_file2, bins=100):
    cv_data1 = np.loadtxt(cv_file1)
    cv_data2 = np.loadtxt(cv_file2)
    cvs1 = cv_data1[:,1]
    cvs2 = cv_data2[:,1]
    rng = min(cvs1.min(), cvs2.min()), max(cvs1.max(), cvs2.max())
    hist1, bin_edges1 = np.histogram(cvs1, bins = bins, range=rng)
    hist2, bin_edges2 = np.histogram(cvs2, bins = bins, range=rng)
    overlap = np.minimum(hist1,hist2)
    area_overlap = overlap.sum()
    prob_overlap = 2*area_overlap/(hist1.sum()+hist2.sum())
    
    # prob_overlap1 = area_overlap/hist1.sum()
    # prob_overlap2 = area_overlap/hist2.sum()
    # bin_edges1_mid = (bin_edges1[1:]+bin_edges1[:-1])/2
    # bin_edges2_mid = (bin_edges2[1:]+bin_edges2[:-1])/2
    # plt.plot(bin_edges1_mid,hist1)
    # plt.plot(bin_edges2_mid,hist2)
    # return prob_overlap1, prob_overlap2
    return prob_overlap
    
num_cv_files = len(cv_files)
prob_overlaps=[]
for i in range(num_cv_files-1):
    prob_overlap = calc_windows_overlap(cv_files[i],cv_files[i+1], bins=100)
    prob_overlaps.append(prob_overlap)
    
prob_overlaps = np.array(prob_overlaps)
plt.plot(prob_overlaps)
plt.ylabel("Overlap ratio")
""")



def post_metad():
    print("""
plumed sum_hills --hills <HILLS_file> --stride <N> --mintozero
""")


def post_driver():
    print("""
# `plumed driver`: post-process the trajectory file <trajectory.xyz> according to 
# the <plumed.dat> file. `plumed driver` supports multiple formats of trajectory files,
# including XYZ, GRO, XTC, TRR, DCD, etc.
plumed driver --plumed <plumed.dat> --ixyz <trajectory.xyz> --trajectory-stride <STRIDE> --timestep <TIMESTEP>
# TIMESTEP: default in ps, 0.001 as 1 fs.
""")   

def handle_pre(option):
    """Handle pre-processing snippets"""
    handlers = {
        "us": pre_us,
        "metad": pre_metad,
        "well-tempered metad": pre_metad,
        "well-tempered metadynamics": pre_metad,
        "metadynamics": pre_metad
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
        "us": post_us,
        "metad": post_metad,
        "well-tempered metad": post_metad,
        "well-tempered metadynamics": post_metad,
        "metadynamics": post_metad,
        "driver": post_driver
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