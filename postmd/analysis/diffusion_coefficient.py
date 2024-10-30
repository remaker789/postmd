from ..avetime import GreenKubo, AveTime
import numpy as np
import warnings
import matplotlib.pyplot as plt

# compute vacf有问题，需要使用dump读取数据
# class VACF(GreenKubo):
    
#     # ?vacf是用三个方向速度平均还是先求出总速度然后系综平均?
#     # todo
#     def __init__(self, path=None, timestep: float=1.0, dim=3):
#         super().__init__(path, timestep)
#         self.dim = dim   # dimension of diffusion coefficient
    
    
#     def run(self, 
#              data_type="lmp",
#              interval:int=None,
#              nlag:int=None,
#              vacf_type="xyz",
#              unit_trans=1.0,
#              int_method="trap"):
#         # todo：所有green-kubo都需要设置interval, 因为并不是所有人都会每1步输出原数据。
#         if interval is None:
#             raise ValueError("Please specify the interval frames between neighboring frames")
#         if data_type == "lmp":
#             # data_type="lmp"
#             # compute myvacf all vacf
#             # fix output_vacf all ave/time 1 1 1 c_myvacf[*] file vacf.time

#             self.dim=len(vacf_type)
#             self.read_file()
#             if self.data.shape[1] < 5:
#                 raise ValueError("You need to output 4 indices of `compute vacf` using the wildcard asterisk '*'")
            
#             # the output of `compute vacf` is v(t)*v(0), we need to divide all data by v(0) to get the true velocity
#             self.data.iloc[:,1:5] = self.data.iloc[:,1:5]/(self.data.iloc[0,1:5]**0.5)
#             if vacf_type == "x":
#                 self.calc_acf(col=1,nlag=nlag,unit_trans=unit_trans)
#             elif vacf_type == "y":
#                 self.calc_acf(col=2,nlag=nlag,unit_trans=unit_trans)
#             elif vacf_type == "z":
#                 self.calc_acf(col=3,nlag=nlag,unit_trans=unit_trans)
#             # elif vacf_type == "xy" or vacf_type == "yx":
#             #     self.calc_acf(col=[1,2],unit_trans=unit_trans)
#             # elif vacf_type == "xz" or vacf_type == "zx":
#             #     self.calc_acf(col=[1,3],unit_trans=unit_trans)
#             # elif vacf_type == "yz" or vacf_type == "zy":
#             #     self.calc_acf(col=[2,3],unit_trans=unit_trans)
#             elif vacf_type == "xyz":
#                 self.calc_acf(col=4,nlag=nlag,unit_trans=unit_trans)
#             else:
#                 raise ValueError(f"the vacf_type '{vacf_type}' is not supported! Please use 'x', 'y', 'z', 'xyz' instead.")
            
#             self.integrate_acf(method=int_method)
#         else:
#             raise ValueError(f"the data type '{data_type}' is not supported! Please use 'lmp' instead.")
        
#     def fit(self, range=None, plot=True):
        
#         # 使用最后1/4数据判断int_acf最后是否收敛
#         from ..utils import judge_plateau
#         if judge_plateau(self.int_acf[int(len(self.int_acf)*0.75):], threshold=0.2):
#             print("The integrated ACF seems converged.")
#         else:
#             print("The integrated ACF seems *not* converged.")
        
#         self.nframe = len(self.int_acf)
#         # 设定range
#         if range is None or range == "last half" or range == "lh":
#             range = [int((self.nframe)/2.), int(self.nframe)+1]
#         elif isinstance(range, (int, float)):
#             range = [int(range), int(self.nframe)+1]
#         else:
#             raise ValueError(f"the range '{range}' is not supported! Please use 'last half', 'middle half', 'all', or a list of two integers [start, stop] instead.")
        
#         print(f"total numbers of MSD lag: {self.nframe}")
#         print(f"Now start linear fitting using range: {range[0]} to {range[1]-1}")
        
#         self.D = self.int_acf[range[0]:range[1]].mean()
#         print(f"The diffusion coefficient in {self.dim} dimension from VACF is {self.D}")
        
#         if plot:
#             fig, axes = plt.subplots(2,1,figsize=(10,5))
#             axes[0].plot(self.nlag, self.acf)
#             axes[0].set_xlabel("$t$ (fs)")
#             axes[0].set_ylabel("VACF")
            
#             axes[1].plot(self.nlag, self.int_acf)
#             axes[1].set_xlabel("$t$ (fs)")
#             axes[1].set_ylabel("integration of VACF")
    



class MSD(AveTime):
    def __init__(self, path=None, timestep: float=1.0):
        super().__init__(path, timestep)
        print("Make sure the trajectory file including the **unwrapped** unscaled coordinates!")


    def run(self,
             data_type="dump",
             interval:int=None, 
             start:int=0, step:int=1, stop:int=None, 
             select="all", calc_dim="xyz",fft=True):
        self.dim=len(calc_dim) # dimension of diffusion coefficient
        self.msd_type = calc_dim
        if data_type == "dump":
            if interval is None:
                raise ValueError("Please specify the interval frames between neighboring frames")
            import MDAnalysis as mda
            import MDAnalysis.analysis.msd as msd
            
            if stop is None:
                warnings.warn("Default set stop time to 1ns. You can sepcify 'stop' by 'stop=xxx' instead")
                stop = int(1e6/interval) # 1ns的MSD
            
            interval_time = self.timestep*interval # actual time between neighboring frames, in fs
            
            print(f"timestep: {self.timestep} fs; interval frames: {interval}; interval time: {interval_time} fs")
            print(f"MSD analysis frames with start:{start}, step:{step}, stop:{stop}")
            # read trajectory file
            u = mda.Universe(self.path,format="LAMMPSDUMP",dt=interval_time*1.e-3, lammps_coordinate_convention="auto",unwrap_images=True) # dt in ps
            
            # run MSD analysis
            _msd = msd.EinsteinMSD(u, select=select, msd_type=calc_dim,fft=fft)
            
            _msd.run(start=start, step=step, stop=stop, verbose=True) # step means the interval of frames
            
            # MSD results
            self.nframe = _msd.n_frames
            self.nlag = np.arange(self.nframe)*interval_time # the lag-time axis, in fs
            self.msd =  _msd.results.timeseries    # in Ang^2
        elif data_type == "lmp":
            # 使用compute msd计算而得
            # compute mymsd all msd
            # fix output_msd all ave/time 1 1 1 c_mymsd[*] file msd.time
            
            # 好像没有必要判断interval，因为step数据已经包含在输出文件中了
            # if interval is None:
            #     raise ValueError("Please specify the interval steps between neighboring output")
            _msd = AveTime(self.path) # 是不是可以删去？
            _msd.read_file() # 是不是可以直接替换成self.read_file()?
            if _msd.data.shape[1] < 5:
                raise ValueError("You need to output 4 indices of `compute msd` using the wildcard asterisk '*'")
            
            # todo: 或许用mapping来区分x,y,z比较好 2024-06-22
            
            self.nframe = len(_msd.data)
            self.nlag = _msd.data["TimeStep"]*self.timestep # in fs
            if calc_dim == "x":
                self.msd = _msd.data.iloc[:,1] # in Ang^2
            elif calc_dim == "y":
                self.msd = _msd.data.iloc[:,2] # in Ang^2
            elif calc_dim == "z":
                self.msd = _msd.data.iloc[:,3] # in Ang^2
            elif calc_dim == "xy" or calc_dim == "yx":
                self.msd = _msd.data.iloc[:,1]+_msd.data.iloc[:,2] # in Ang^2
            elif calc_dim == "xz" or calc_dim == "zx":
                self.msd = _msd.data.iloc[:,1]+_msd.data.iloc[:,3] # in Ang^2
            elif calc_dim == "yz" or calc_dim == "zy":
                self.msd = _msd.data.iloc[:,2]+_msd.data.iloc[:,3] # in Ang^2
            elif calc_dim == "xyz":
                self.msd = _msd.data.iloc[:,4] # in Ang^2
            else:
                raise ValueError(f"the msd_type '{calc_dim}' is not supported! Please use 'x', 'y', 'z', 'xy', 'yx', 'xz', 'zx', 'yz', 'zy', 'xyz' instead.")
                
        
        else: 
            raise ValueError(f"the data type '{data_type}' is not supported! Please use 'dump' or 'lmp' instead.")
            
        
    
    def fit(self, range=None, unit_trans=1.0, plot=True):
        # MSD(t)拟合
        
        # 设定range
        if range is None or range == "last half" or range == "lh":
            range = [int((self.nframe)/2.), int(self.nframe)+1]
        elif range == "middle half" or range == "mh":
            range = [int((self.nframe)/4.), int((self.nframe)/4.*3.)+1]
        elif range == "all":
            range = [0, int(self.nframe)+1]
        elif isinstance(range, (int, float)):
            range = [int(range), int(self.nframe)+1]
        elif isinstance(range, (list, tuple, np.ndarray)):
            if len(range)==2:
                range = [int(range[0]), int(range[1])+1]
            else:
                raise ValueError(f"the range '{range}' is not supported! Please use a list of two integers [start, stop] instead.")

        else:
            raise ValueError(f"the range '{range}' is not supported! Please use 'last half', 'middle half', 'all', or a list of two integers [start, stop] instead.")
        
        print(f"total numbers of MSD lag: {self.nframe}")
        print(f"Now start linear fitting using range: {range[0]} to {range[1]-1}")
        slope, _intcpt= np.polyfit(self.nlag[range[0]:range[1]], self.msd[range[0]:range[1]], 1)
        self.D = slope/(2.0*self.dim)*unit_trans
        self.diffusion_coeff = self.D
        print(f"The diffusion coefficient in {self.dim} dimension from MSD is {self.D}")
        if plot:
            fig, ax = plt.subplots()
            ax.plot(self.nlag, self.msd, label="MSD")
            ax.plot(self.nlag[range[0]:range[1]], slope*self.nlag[range[0]:range[1]]+_intcpt,'r--', label="linear fit", lw=3.0)
            ax.set_xlabel("$t$ (fs)")
            ax.set_ylabel("$\mathrm{MSD}$"+f"$_{{{self.msd_type}}}$"+" ($\mathrm{\AA^2}$)")
            ax.legend()
            
            