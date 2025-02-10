import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mathtext import _mathtext as mathtext
# mathtext.FontConstantsBase = mathtext.ComputerModernFontConstants

os.chdir(os.path.dirname(__file__)) # 如果是ipynb文件请注释此行。




import postmd

os.chdir("output")


def model(filename):
    x = np.arange(1,11,1)
    y = x**2
    plt.figure()
    plt.plot(x,y, '-o', label="legend")
    plt.text(2,14, '5aKPS.')
    plt.text(2,50, r'x$^2$+y$^2$=1')
    plt.xlabel(r"$E$ [V/$\mathrm{\AA}$]")
    plt.ylabel(r"$I_c$ [A]")
    plt.title('title')
    plt.legend()
    plt.savefig(filename)
    

styles = [
    ["classic"],
    ["science"],
    ["science", "no-latex"],
    ["Times+stix"],
    ["Arial"],
    ["Times"],
]



for style in styles:
    with plt.style.context(style):
        if len(style) == 1:
            outname = f"{style[0]}.pdf"
        else:
            outname = f"{'+'.join(style)}.pdf"
        if "Times" or "Arial" in style:
            mathtext.FontConstantsBase = mathtext.ComputerModernFontConstants
        model(outname)
