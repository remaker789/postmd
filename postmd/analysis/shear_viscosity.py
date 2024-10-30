from ..avetime import GreenKubo, AveTime
import numpy as np
import pandas as pd
import warnings
from ..utils import mapdim2col
from .base import Results

class ShearViscosity(GreenKubo):
    """The shear viscosity is calculated by Green-Kubo formula,
    .. math::
    
        \eta = \frac { V } { k _ { B } T } \int _ { 0 } ^ { \infty } d t \langle \tau _ { \alpha \beta } ( t ) \tau _ { \alpha \beta } ( 0 ) \rangle _ { t _ { 0 } }
    
    where :math:`V` is the volume of the system, :math:`k _ { B }` is the Boltzmann constant, :math:`T` is the temperature, and :math:`\tau _ { \alpha \beta }` is the off-diagonal(or traceless) components of stress tensor.
    """    
    # 2024-06-22 zss 剪切粘度计算需要分别计算xy,yz,xz方向，然后求和平均得到整体平均值。
  

    def __init__(self, path=None, timestep: float = 1):
        super().__init__(path, timestep)
    
    
    def run(self, nlag=None, mapping=None, unit_trans=1.0,int_method="trap"):
        # now we only support the calculation of the total shear viscosity.
        # 一般只计算xy,yz,xz三个应力张量，有的还计算1/2(tau_xx-tau_yy)和1/2(tau_yy-tau_zz)
        self.read_file()
        self.results=Results()      
        if mapping is None:
            raise ValueError("The mapping of the stress tensor to indexing is not specified. You should specify it using column names like mapping={'xy':'v_pxy','yz':'v_pyz','xz':'v_pxz'} or using column numbers like mapping={'xy':1,'yz':2,'xz':3}")
        
        acf_s = []
        int_acf_s = []
        
        for comp, col in mapping.items():
            self.calc_acf(data_type="raw",col=col,nlag=nlag,unit_trans=unit_trans)
            self.integrate_acf(method=int_method)
            self.results[comp]=Results({'nlag':self.nlag, 'time':self.nlag*self.timestep, 'acf':self.acf,'int_acf':self.int_acf})
            acf_s.append(self.acf)
            int_acf_s.append(self.int_acf)
        
        mean_acf = np.array(acf_s).mean(axis=0)
        std_acf = np.array(acf_s).std(axis=0,ddof=1)
        mean_int_acf = np.array(int_acf_s).mean(axis=0)
        std_int_acf = np.array(int_acf_s).std(axis=0,ddof=1)
        self.results["mean"]=Results({'nlag':self.nlag, 'time':self.nlag*self.timestep, 'acf':mean_acf, 'std_acf':std_acf, 'int_acf':mean_int_acf, 'std_int_acf':std_int_acf})
        
    def write_results(self, filename="shear_viscosity.xlsx"):
        df = pd.DataFrame.from_dict({(i,j): self.results[i][j] 
                           for i in self.results.keys() 
                           for j in self.results[i].keys()},
                        orient='columns')
        if filename.endswith("xlsx"):
            df.to_excel(filename)
        else:
            df.to_csv(filename, index=False)
    
    