import numpy as np
import pandas as pd
import scipy.constants as C
import warnings
from ..avetime import AveTime
from ..utils import mapdim2col
from .base import Results

# $$
# \varepsilon_\mathrm{bulk} = 1 + \frac{(\langle M^2\rangle-\langle M \rangle^2)}{3\varepsilon_0 V k_B T} \tag{1}
# $$
# lammps

class DielectricHomo(AveTime):
    # only for homogeneous system with total charge as 0. For the one in interface and confinement, see https://maicos-devel.gitlab.io/maicos/explanations/dielectric.html
    def __init__(self, path=None, timestep: float=1.0, temperature=None, cell=None, V=None):
        super().__init__(path, timestep)
        if temperature is None:
            raise ValueError("Please specify the temperature!")
        self.T = temperature #  temperature
        self.results=Results()


    def run(self,
            volumn=None,
            calc_dim="xyz",
            dim_map=None,
            unit_trans=1.0,
            ):
        # todo: 对于NVT来说，volumn是固定的，但是如果是NPT呢？volumn应该是平均值，还是在每一个时刻先除以V然后再平均呢？

        # compute mydipole all dipole 
        # fix output_dipole all ave/time 1 1 1 c_mydipole c_mydipole[*] file dipole.time
        # 准确来说，并不需要输出c_mydipole，只需要c_mydipole[*]即可
        # compute dipole实际使用的unwrapped coords计算的dipole!
        df = self.read_file()
        cols = mapdim2col(dim_map, calc_dim)
        M2_s = np.zeros(len(cols))
        M_s = np.zeros(len(cols))
        for i, col in enumerate(cols):
            data_M = df.iloc[:,col] if type(col)==int else df.loc[:,col]
            # data_M = data_M*unit_trans
            M2_s[i] = np.mean(data_M**2)
            M_s[i]= np.mean(data_M)
        

        fluct_s = M2_s-M_s**2

        kB = C.Boltzmann
        eps_0 = C.epsilon_0
        # unit_trans = e*1e-10
        eps_s = 1 + fluct_s*(unit_trans**2)/(eps_0*volumn*kB*self.T)
        
        self.results.M2 = M2_s
        self.results.M = M_s
        self.results.fluct = fluct_s
        self.results.eps = eps_s
        self.results.eps_mean = np.mean(eps_s)
        
        
        # self.results= {
        #     'M2': M2_s,
        #     'M': M_s,
        #     'fluct':fluct_s,
        #     'eps': eps_s,
        #     'eps_mean': np.mean(eps_s)
        # }
        print("The sequence of results is as the parameter 'calc_dim'. For 'yxz', the elements of array is y,x,z")
        print(self.results)


        # elif data_type == "dump":
        #     ## MDAnalysis
        #     # slow
        #     print("Make sure the trajectory file including the **unwrapped** unscaled coordinates!")
        #     import MDAnalysis as mda
        #     from MDAnalysis.analysis.dielectric import DielectricConstant
        #     u = mda.Universe("system.data", "dump.lammpstrj",topology_format="DATA", format="LAMMPSDUMP",dt=YOUR_DT, lammps_coordinate_convention="auto", unwrap_images=True)

        #     if charges is None:
        #         raise ValueError("Please provide the charges of the atoms!")
            
        #     # 2024-06-23 zss 从system.data中的读取电荷信息就不需要手动设置电荷了            
        #     # u.add_TopologyAttr('charges')
        #     # for key, value in charges.items():
        #     #     # 因为lammps的lammpstrj文件不带电荷，所以需要手动设置电荷。
        #     #     # 现在只支持单个选择 电荷均相同的设置。

        #     #     _select = u.select_atoms(key)
        #     #     _select.atom.charges =  value * np.ones(_select.atoms.n_atoms)
            
        #     # oxygen = u.select_atoms("type 1")
        #     # hydrogen = u.select_atoms("type 2")

        #     # oxygen.atoms.charges = -0.8476 * np.ones(oxygen.atoms.n_atoms)
        #     # hydrogen.atoms.charges = 0.4238 * np.ones(hydrogen.atoms.n_atoms)
            
        #     dielec = DielectricConstant(u.atoms,temperature=self.T, make_whole=False, verbose=True)
        #     dielec.run()
        #     # 需要注释掉dielectric.py文件中if not np.allclose(self.atomgroup.total_charge(compound='fragments')这个判断语e句，因为lammps的dump文件不包含fragments.
        #     self.results = dielec.results
        #     print(self.results)

        # else: 
        #     raise ValueError(f"the data type '{data_type}' is not supported! Please use 'dump' or 'lmp' instead.")