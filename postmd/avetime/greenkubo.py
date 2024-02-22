import numpy as np
import statsmodels.api as sm
from scipy.integrate import simpson, trapezoid

from .avetime import AveTime


# 一般来说GreenKubo是对平衡态下某些参数的自相关函数进行时间上的积分。
class GreenKubo(AveTime):
    """calculate the auto-correlation function (acf) and corresponding integration
    for Green-Kubo formula. But you need to pay attention to the true equation of
    Green-Kubo formula, because this script just do the integration of acf, not
    including other parts, like multiply or divide by some properties.
    """    
    def __init__(self, T=298.0, timestep=1.0):
        super().__init__(T, timestep)
        self.acf=None
        self.nlag=None
        self.int_acf=None

    # # after using fix ave/correlate ave running
    # def integrate(self, *, Nevery:int=None, Nrepeat:int=None, Nfreq:int=None, 
    #               overwrite:bool=True, col_time:int or str = 1, col_data:int or str=None):
    #     print('''The data must be genetated using "ave running" in fix ave/correlate commmand''')
        
    #     if overwrite:
    #         # ! 不确定fix ave/correlate中的names是第几行。
    #         df = self.read_file(self.path, names=None, num_line=2)
    #     else:
    #         pass     # todo: 不是overwrite的数据，只需要取最后一个block
    
    
    def cal_acf(self, data_type="raw", col:int=None, nlag:int=None, nlag_col:int=1, unit_trans=1.0) -> np.ndarray:
        """calculate the acf from data file and store in ``self.nlag`` and ``self.acf``.

        Args:
            data_type (str, optional): the data type of self.data (raw or acf). Defaults to ``"raw"``.
            col (int, optional): the column number(start from 0) of data to process in self.data. Defaults to ``None``.
            nlag (int, optional): Limit the number of autocovariances returned.  Size of returned array is nlag + 1. Defaults to ``None``.
            nlag_col (int, optional): the column number(start from ``0``) of nlag in data, usually named "TimeDelta". Defaults to ``1``.
            unit_trans (float, optional): transform the unit of property to SI unit. Defaults to ``1.0``.

        Returns:
            None
        """
        if col is None:
            raise ValueError("Please specify the column number of data you want to process!")
        
        # if the data are generated from "fix ave/correlate type auto overwriting"
        if data_type=="acf":
            print("Please ensure that your data are generated from 'fix ave/correlate type auto overwriting' command!")
            self.acf=np.array(self.data.iloc[:,col]*(unit_trans**2))
            self.nlag=np.array(self.data.iloc[:,nlag_col])
        
        # if the data is just the raw data, not processed by acf
        elif data_type=="raw":
            print("Please ensure your data are generate as raw data!")
            # in lammps, the data in step=0 of run is usually unstable
            print("The script will drop out first line due to unstable data in the initial step.")
            self.acf = sm.tsa.stattools.acovf(self.data.iloc[1:,col]*unit_trans, nlag=nlag)
            self.nlag=np.arange(0,len(self.acf))
        else:
            raise ValueError(f"the data type '{data_type}' is not supported! Please use 'acf' or 'raw' instead.")
            
        
            
    
    # 这里的integrate只是用来对acf进行积分
    def integrate_acf(self, nlag=None, acf=None, *, method="simpson"):
        """Integrate the acf data to the Green-Kubo formula and store in ``self.int_acf``.

        Args:
            nlag (np.ndarray, optional): the nlag data. Defaults to ``None``.
            acf (np.ndarray, optional): the acf data. Defaults to ``None``.
            method (str, optional): the integration method (``"simpson"`` or ``"trap"``). Defaults to ``"simpson"``.
        """        
        # 我想建一个输出所有method的方法
        print(f"-------- Integrate ACF using {method} method ---------")
        acf = acf if acf else self.acf
        nlag = nlag if nlag else self.nlag
        
        if acf is None:
            raise ValueError("Please input acf!")
        if nlag is None:
            raise ValueError("Please input nlag!")
        int_acf=[]        
        
        for i in range(len(acf)):
            if method == "simpson":    
                _int = simpson(y=acf[:i+1], x=nlag[:i+1])
            elif method == "trap":
                _int = trapezoid(y=acf[:i+1], x=nlag[:i+1])
            else:
                raise ValueError(f"method '{method}' is not supported yet!")
            int_acf.append(_int)
            
        self.int_acf=np.array(int_acf)*self.timestep*1e-15
    
    def write_result(self,filename="gk.csv"):
        """write the data of nlag, acf and integration of acf to file

        Args:
            filename (str, optional): the filename for output. Defaults to ``"gk.csv"``.
        """               
        np.savetxt(f"acf-{filename}", np.c_[self.nlag, self.acf, self.int_acf], delimiter=',', header="nlag,acf,int_acf")
        
        
    # 2023-08-28 校正有限尺寸的摩擦系数，因为green-kubo计算的摩擦系数会降低，所以我们需要校正
    # 目前是根据Haruki Oga_2019_JCP的文章校正的    
    #pending 不确定会不会做。。。
    def finite_size_correct_friction_coeff(self):
        pass
        
    

    