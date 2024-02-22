import re
from pathlib import Path

import numpy as np
import pandas as pd

from ..system import System


class AveTime(System):
    """process the data generated from "fix ave/time" command in LAMMPS. 
    (see https://docs.lammps.org/fix_ave_time.html)
    """

    def __init__(self, T=298, timestep=1):
        """initial the class
        Args:
            T (float, optional): temperature of system. Defaults to ``298.0`` [K].
            timestep (float, optional): temperature you set in LAMMPS input file. Defaults to ``1.0`` [fs].
        """
        super().__init__(T, timestep)
        self.path=None

    def set_path(self, path=None):
        self.path = path
        self.judge_path()
        print(f"You are processing the file: '{self.path}'")

    def judge_path(self):
        """judge whether the path of object is accessable.
        """
        if (self.path is None):  # 判断path或者self.path是否为空
            raise ValueError("path can not be empty!")
        elif not Path(self.path).is_file():
            raise ValueError(f"'{self.path}' is not a file!")

    def get_names(self, name_line:int) -> list:
        """get the names used in dataframe when read files
        
        Args:
            name_line (int, optional): the line number (start from ``0``) including the column names of data. Defaults to ``2``.
        
        Returns:
            list: a list of the column names of data
        """
        with open(self.path) as f:
            for i, line in enumerate(f.readlines(), 1): # enumerate(iterator,1) means counting from 1
                if i == name_line:
                    return re.split('\s+', line.strip())[1:] # split the line with regex '\s+' to list, and discard the first element "#"

    def read_file(self, path:str=None,end=None, name_line:int=2, names:list=None, **kwargs) -> pd.DataFrame:
        """read file and return DataFrame
        
        There are some default setting in read files:

        - ``comment="#"`` due to the nature of LAMMPS
        - ``sep="\s+"`` to match all space between data
        - ``**kwargs`` received parameters of read_csv and read_excel, depending on your cases. 

        Args:
            path (str): path to the file. Defaults to ``None``.
            names (list, optional): list of columns names of data. 
                                    Defaults to ``None``, meaning using the content in name_line of file.
            name_line (int, optional): when names=None, the content in ``line=<name_line>``
                                      will be used as the column names. Defaults to ``2``.

        Returns:
            pd.DataFrame: DataFrame object read from the file.
        """
        path = path if path else self.path
        self.set_path(path)
        names = names if names else self.get_names(name_line)
        df = pd.read_csv(self.path, names=names, sep='\s+', comment="#", **kwargs)
        if end is None:
            self.data = df
        else:
            self.data = df[:end]
        return df

