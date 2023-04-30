import numpy as np
import statsmodels.api as sm
from scipy.integrate import simpson, trapezoid

from .AveTime import *


# 一般来说GreenKubo是对平衡态下某些参数的自相关函数进行时间上的积分。
class GreenKubo(AveTime):
    """calculate the auto-correlation function (acf) and corresponding integration
    for Green-Kubo formula. But you need to pay attention to the true equation of
    Green-Kubo formula, because this script just do the integration of acf, not
    including other parts, like multiply or divide by some properties.

    Args:
        AveTime: parent class that was used to process data from "fix ave/time" command 
    """    
    def __init__(self, T=298, timestep=1):
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
    
    
    def cal_acf(self, data_type="raw", col:int=None, nlag:int=None, nlag_col:int=1) -> np.ndarray:
        """calculate the acf from data file.

        Args:
            data_type (str, optional): the data type of self.data. Defaults to "raw".
            col (int, optional): the column number(start from 0) of data to process in self.data. Defaults to None.
            nlag (int, optional): Limit the number of autocovariances returned.  Size of returned array is nlag + 1. Defaults to None.

        Raises:
            ValueError: _description_

        Returns:
            _type_: _description_
        """
        if col is None:
            raise ValueError("Please specify the column number of data you want to process!")    
        if data_type=="acf":
            print("Please ensure that your data are generated from 'fix ave/correlate type auto overwriting' command!")
            # if the data are generated from "fix ave/correlate type auto overwriting"
            self.acf=np.array(self.data.iloc[:,col])
            self.nlag=np.array(self.data.iloc[:,nlag_col])
            
        elif data_type=="raw":
            print("Please ensure your data are generate as raw data!")
            # if the data is just the raw data, not processed by acf
            self.acf = sm.tsa.stattools.acovf(self.data.iloc[:,col], nlag=nlag)
            self.nlag=np.arange(0,len(self.acf))
        else:
            raise ValueError(f"the data type '{data_type}' is not supported! Please use 'acf' or 'raw' instead.")
            
        
            
    
    # 这里的integrate只是用来对acf进行积分
    def integrate_acf(self, nlag=None, acf=None, *, method="simpson"):
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
            
        self.int_acf=np.array(int_acf)
    

    