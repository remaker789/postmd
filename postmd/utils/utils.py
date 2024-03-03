import os
import shutil
from pathlib import Path
import numpy as np
import scipy
import scipy.constants as C



def average_replicates(data_arrays):
    """average the data from replicates.

    Args:
        data_arrays (list or np.ndarray): a list of data

    Returns:
        np.ndarray: the averaged data

    Examples:
        >>> import numpy as np
        >>> from postmd.utils import average_replicates
        >>> data_arrays = [
        >>>    np.array([4.3, 5.6, 3.8, 5.1, 4.9]),  # First dataset
        >>>    np.array([3.2, 4.5, 4.1, 3.7, 4.3]),  # Second dataset
        >>>    np.array([5.5, 6.2, 5.9, 6.1, 5.8]),   # Third dataset
        >>>    ]
        >>> # Calculate the average of each dataset (i.e., each array).
        >>> average = average_replicates(data_arrays)
        >>> print(f"Averages of replicates: {average}")
        Averages of replicates: [4.33333333 5.43333333 4.6        4.96666667 5.        ]
    """    
    data_arrays = np.array(data_arrays)
    # Ensure that the input is a list of numpy arrays.
    if data_arrays.ndim<=1:
        raise TypeError("Input should be a list of numpy arrays.")

    # Use numpy's mean function to calculate the average of each dataset (i.e., each array).
    average = np.mean(data_arrays, axis=0)

    return average

def cal_box_length(num, density=1.0, NA=None):
    """calculate the length of a cubic water box.

    Warning: 
        The built-in Avogadro constant in LAMMPS (units real or metal) is 6.02214129e23, see lammps/src/update.cpp, they write "force->mv2d = 1.0 / 0.602214129" for units real and units metal. **However, we defalutly used the Avogadro constant in scipy.constants is 6.022140857e23, which is the international standard.**


    Args:
        num (int): the number of water molecules.
        density (float, optional): the density of water, in g/cm^3. Defaults to 1.0.

    Returns:
        float: the length of a cubic water box, in Angstrom
        
    Examples:
        >>> import postmd.utils as utils
        >>> utils.cal_box_length(1000, density=1.0)
        The length of a cubic water box for 1000 water molecules and 1.0 g/cm^3 is 31.043047 Angstrom
    """
    if NA is None:
        NA = C.Avogadro

    mass_O=15.9994
    mass_H=1.008
    mass_water = mass_O+2*mass_H
    box_length = (num*density/(1/mass_water*NA))**(1/3) * 1e8
    print(f"The length of a cubic water box for {num} water molecules and {density} g/cm^3 is {box_length:.6f} Angstrom")


def create_dir(path, backup=True):
    """Create a directory at the specified 'path'. If the directory already exists and 'backup' is True, rename the original directory by appending '.bkXXX'.

    Args:
        path (str): The path of the directory to be created.
        backup (bool, optional): Whether to back up an existing directory. Default is ``True``.

    Examples:
        >>> import os
        >>> import postmd.utils as utils
        >>> 
        >>> print(os.listdir())
        ['createdir.py']
        >>> utils.create_dir("test")
        >>> print(os.listdir())      # create a new "test" dir
        ['createdir.py', 'test']
        >>> utils.create_dir("test") 
        >>> print(os.listdir())      # move orgin "test" dir to "test.bk000" dir
        ['createdir.py', 'test', 'test.bk000']
        >>> utils.create_dir("test")
        >>> print(os.listdir())      # move orgin "test" dir to "test.bk001" dir
        ['createdir.py', 'test', 'test.bk000', 'test.bk001']

    """
    path += "/"
    if os.path.isdir(path) and backup:
        dirname = os.path.dirname(path)
        counter = 0
        while True:
            bk_dirname = dirname + ".bk%03d" % counter # formatting .bkxxx
            if not os.path.isdir(bk_dirname):
                shutil.move(dirname, bk_dirname)
                break
            counter += 1
    os.makedirs(path)
    
    
def cumave(data):
    """calculate the cumulative average.

    Args:
        data (1d list): the data need to do the cumulative average

    Returns:
        np.ndarray: the cumulative average
    
    Examples:
        >>> import postmd.utils as utils
        >>> import numpy as np
        >>> array = np.arange(9)
        >>> print(array)
        [0 1 2 3 4 5 6 7 8]
        >>> utils.cumave(array)
        array([0. , 0.5, 1. , 1.5, 2. , 2.5, 3. , 3.5, 4. ])
            
    """    
    data = np.array(data)
    return np.cumsum(data)/np.arange(1,len(data)+1)





def stat_bin(x,y, bins=10,range=None):
    """statistic the mean and standard deviation(ddof=0) of x and y in each bin.
    Here we used the [scipy.stats.binned_statistic](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.binned_statistic.html) function.

    Args:
        x: (N,) array_like. A sequence of values to be binned.
        y: (N,) array_like. The data on which the statistic will be computed. This must be the same shape as x, or a set of sequences - each the same shape as x. If values is a set of sequences, the statistic will be computed on each independently.
        bins (int or sequence of scalars, optional): If bins is an int, it defines the number of equal-width bins in the given range (10 by default). If bins is a sequence, it defines the bin edges, including the rightmost edge, allowing for non-uniform bin widths. Values in x that are smaller than lowest bin edge are assigned to bin number 0, values beyond the highest bin are assigned to bins[-1]. If the bin edges are specified, the number of bins will be, (nx = len(bins)-1). Defaults to 10.
        range ((float, float) or [(float, float)], optional): The lower and upper range of the bins. If not provided, range is simply (x.min(), x.max()). Values outside the range are ignored. Defaults to None.

    Returns:
        tuple: (x_mean, x_std, y_mean, y_std)
    """    
    x_mean, _, _ = scipy.stats.binned_statistic(x, x, statistic='mean', bins=bins, range=range)
    x_std, _, _ = scipy.stats.binned_statistic(x, x, statistic='std', bins=bins, range=range)
    y_mean, _, _ = scipy.stats.binned_statistic(x, y, statistic='mean', bins=bins, range=range)
    y_std, _, _ = scipy.stats.binned_statistic(x, y, statistic='std', bins=bins, range=range)
    return x_mean, x_std, y_mean, y_std


def judge_file(path):
    """judge whether the path is a file.
    """

    path = os.path.abspath(path)
    if not Path(path).is_file():
        raise ValueError(f"'{path}' is not a file!")  
    

def judge_dir(path):
    """judge whether the path is a folder.
    """
    path = os.path.abspath(path)
    if not Path(path).is_dir():
        raise ValueError(f"'{path}' is not a dir!") 
    

