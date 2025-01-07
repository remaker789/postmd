import os
import shutil
from pathlib import Path
import numpy as np
import scipy
import scipy.constants as C
from functools import wraps


fig_ext = ["png", "pdf", "jpg", "jpeg", "gif", "bmp", "svg", "eps", "ps", "tif", "tiff", "xbm", "xpm", "xwd", "png", "pdf", "jpg", "jpeg", "gif", "bmp", "svg", "eps", "ps", "tif", "tiff", "xbm", "xpm", "xwd"] # common extension of figures
video_ext = ["mp4", "avi", "mov", "mpeg", "mpg", "wmv", "mkv", "flv", "webm", "gif"] # common extenstion of videos


def mapdim2col(dim_map:dict, dim:str):
    # dim = "xyz"
    # dim_map = {'x': 'value1', 'y': 'value2', 'z': 'value3'}  # 假设这是dim_map字典
    # return mapped_col=['value1', 'value2', 'value3']
    
    # 将dim中的每个字符映射到dim_map中的值
    mapped_col = [dim_map[char] for char in dim]
    return mapped_col




def calc_replicas_mean_std(data_arrays,ddof=0):
    """average the data from replicates.

    Args:
        data_arrays (list or np.ndarray): a list of data

    Returns:
        np.ndarray: the averaged data

    Examples:
        >>> import numpy as np
        >>> from postmd.utils import calc_replicas_mean_std
        >>> data_arrays = [
        >>>    np.array([4.3, 5.6, 3.8, 5.1, 4.9]),  # First dataset
        >>>    np.array([3.2, 4.5, 4.1, 3.7, 4.3]),  # Second dataset
        >>>    np.array([5.5, 6.2, 5.9, 6.1, 5.8]),   # Third dataset
        >>>    ]
        >>> # Calculate mean and std of each dataset (i.e., each array).
        >>> mean, std = calc_replicas_mean_std(data_arrays)
        >>> print(f"replicas mean: {mean}")
        Averages of replicates: [4.33333333 5.43333333 4.6        4.96666667 5.        ]
    """    
    data_arrays = np.array(data_arrays)
    # Ensure that the input is a list of numpy arrays.
    if data_arrays.ndim<=1:
        raise TypeError("Input should be a list of numpy arrays.")

    _mean = np.mean(data_arrays, axis=0)
    _std = np.std(data_arrays, axis=0, ddof=ddof)

    return _mean, _std

def calc_box_length(num, density=1.0, NA=None):
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
        >>> utils.calc_box_length(1000, density=1.0)
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




def backup(func):
    # 用于作为写入文件备份的装饰器。目前目录备份还未实现。
    @wraps(func)
    def decorated(*args, **kwargs):
        if len(args)>0:
            print(1)
            path = args[0]
        else:
            print(kwargs)
            print(2)
            path = kwargs["path"]# if "path" in kwargs.keys else kwargs["filepath"]
        
        path = os.path.abspath(path)  # path = dname + fname
        fname = os.path.basename(path)
        dname = os.path.dirname(path)
        
        fext = os.path.splitext(fname)[1][1:]# file extension
        
        counter = 0
        if os.path.exists(path):
            while True:
                if fext in fig_ext:
                    bk_fname = os.path.splitext(fname)[0] + ".bk%03d" % counter + os.path.splitext(fname)[1]
                else:
                    bk_fname = fname + ".bk%03d" % counter # formatting .bkxxx
                if not os.path.exists(bk_fname):
                    shutil.move(path, os.path.join(dname,bk_fname))
                    print(f"'{fname}' is backup to '{bk_fname}'")
                    break
                counter += 1
        return func(*args, **kwargs)
    return decorated
    
def cummean(data):
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
        >>> utils.cummean(array)
        array([0. , 0.5, 1. , 1.5, 2. , 2.5, 3. , 3.5, 4. ])
            
    """    
    data = np.array(data)
    return np.cumsum(data)/np.arange(1,len(data)+1)





def stats_mean_std_bins(x,y, bins=10,range=None):
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
    


# 效率很低
# def calc_green_kubo(data, nlag=None):
#     """calculate the auto-correlation function in the Green-Kubo formula

#     Args:
#         data (_type_): _description_
#         nlag (_type_, optional): _description_. Defaults to None.

#     Returns:
#         _type_: _description_
#     """    
#     print("Begin calculate the auto-correlation function in the Green-Kubo formula")
#     data = np.array(data)
#     n = len(data)
#     acf=np.zeros(nlag)
#     acf[0]=(data**2).sum()
#     for i in range(1,nlag):
#         for j in range(n-i):
#             acf[i]+=data[j]*data[j+i]
#     acf = acf / (n - np.arange(nlag))
#     return acf
    
# copy and paste from statsmodels/compat/scipy.py
def _next_regular(target):
    """
    Find the next regular number greater than or equal to target.
    Regular numbers are composites of the prime factors 2, 3, and 5.
    Also known as 5-smooth numbers or Hamming numbers, these are the optimal
    size for inputs to FFTPACK.

    Target must be a positive integer.
    """
    if target <= 6:
        return target

    # Quickly check if it's already a power of 2
    if not (target & (target - 1)):
        return target

    match = float("inf")  # Anything found will be smaller
    p5 = 1
    while p5 < target:
        p35 = p5
        while p35 < target:
            # Ceiling integer division, avoiding conversion to float
            # (quotient = ceil(target / p35))
            quotient = -(-target // p35)
            # Quickly find next power of 2 >= quotient
            p2 = 2 ** ((quotient - 1).bit_length())

            N = p2 * p35
            if N == target:
                return N
            elif N < match:
                match = N
            p35 *= 3
            if p35 == target:
                return p35
        if p35 < match:
            match = p35
        p5 *= 5
        if p5 == target:
            return p5
    if p5 < match:
        match = p5
    return match

   
def calc_acf(data, nlag=None, fft=False):
    """calculate the auto-correlation function in the Green-Kubo formula

    Args:
        data (_type_): _description_
        nlag (_type_, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """    
    print("Begin calculate the auto-correlation function in the Green-Kubo formula")
    
    # forked from acovf function of statsmodels package
    data = np.array(data,dtype=np.float64)
    n = len(data)
    if fft is False:
        acf = np.empty(nlag + 1)
        acf[0] = data.dot(data)
        for i in range(nlag):
            acf[i + 1] = data[i + 1 :].dot(data[: -(i + 1)])
    else:
        # find the next regular number considering *nlag* zero-padding
        n_fft = _next_regular(n+nlag+1)
        F_data = np.fft.fft(data, n=n_fft)
        acf = np.fft.ifft(F_data * np.conjugate(F_data))[:nlag+1] # Wiener-Khinchin theorem
        acf = acf.real
        
    acf = acf / (n - np.arange(nlag+1)) # averaged by the number of samples
    return acf
    

def judge_plateau(data, threshold=0.2):
    # 判断一段数据是否会偏离水平线太多
    data = np.array(data)
    _max = np.max(data)
    _min = np.min(data)
    _mean = np.mean(data)
    if (_max-_mean)/_mean < threshold and (_min-_mean)/_mean < threshold:
        return True
    else:
        return False
