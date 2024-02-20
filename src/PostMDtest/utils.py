import os
import shutil
import numpy as np
import scipy.constants as C
import scipy

NA = C.Avogadro

def create_path(path, backup=False):
    """创建路径为path的文件夹。如果文件夹已存在，则将原文件夹备份，并添加.bkXXX后缀

    Args:
        path (_type_): _description_
    """    
    path += "/"
    if os.path.isdir(path) and backup:
        dirname = os.path.dirname(path)
        counter = 0
        while True:
            bk_dirname = dirname + ".bk%03d" % counter
            if not os.path.isdir(bk_dirname):
                shutil.move(dirname, bk_dirname)
                break
            counter += 1
    os.makedirs(path)
    
    
def cum_ave(data):
    """calculate the cumulative average

    Args:
        data (1d list): the data need to do the cumulative average

    Returns:
        np.ndarray: the cumulative average
    """    
    data = np.array(data)
    return np.cumsum(data)/np.arange(1,len(data)+1)


def cal_box_length(num, density=1.0):
    """calculate the length of a cubic water box

    Args:
        num (int): the number of water molecules.
        density (float, optional): the density of water, in g/cm^3. Defaults to 1.0.

    Returns:
        float: the length of a cubic water box, in Angstrom
    """
    mass_O=15.9994
    mass_H=1.008
    mass_water = mass_O+2*mass_H
    box_length = (num*density/(1/mass_water*NA))**(1/3) * 1e8
    print(f"The length of a cubic water box for {num} water molecules and {density} g/cm^3 is {box_length:.4f} Angstrom")
    return box_length


def stat_bin(x,y, bins=10,range=None):
    """statistic the mean and standard deviation(ddof=0) of x and y in each bin.
    Here we used the [scipy.stats.binned_statistic](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.binned_statistic.html) function.

    Args:
        x: (N,) array_like. A sequence of values to be binned.
        y: (N,) array_like. The data on which the statistic will be computed. This must be the same shape as x, or a set of sequences - each the same shape as x. If values is a set of sequences, the statistic will be computed on each independently.
        bins (int or sequence of scalars, optional): If bins is an int, it defines the number of equal-width bins in the given range (10 by default). If bins is a sequence, it defines the bin edges, including the rightmost edge, allowing for non-uniform bin widths. Values in x that are smaller than lowest bin edge are assigned to bin number 0, values beyond the highest bin are assigned to bins[-1]. If the bin edges are specified, the number of bins will be, (nx = len(bins)-1). Defaults to 10.
        range ((float, float) or [(float, float)], optional): The lower and upper range of the bins. If not provided, range is simply (x.min(), x.max()). Values outside the range are ignored. Defaults to None.

    Returns:
        _type_: _description_
    """    
    x_mean, _, _ = scipy.stats.binned_statistic(x, x, statistic='mean', bins=bins, range=range)
    x_std, _, _ = scipy.stats.binned_statistic(x, x, statistic='std', bins=bins, range=range)
    y_mean, _, _ = scipy.stats.binned_statistic(x, y, statistic='mean', bins=bins, range=range)
    y_std, _, _ = scipy.stats.binned_statistic(x, y, statistic='std', bins=bins, range=range)
    return x_mean, x_std, y_mean, y_std