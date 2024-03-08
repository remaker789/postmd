import numpy as np
from postmd.animate import CumAnimate

import os

os.chdir(os.path.dirname(__file__)) # 如果是ipynb文件请删除此行。
if os.path.basename(os.path.abspath('.')) == 'output':
    pass
elif not os.path.exists('./output'):
    os.makedirs('./output')
os.chdir('./output') # 请注意，现在工作路径是output！
print(f"Note: Current working directory is '{os.getcwd()}'")


if __name__ == '__main__':
    x_data = np.linspace(0, 10, 100)  # Corresponding x-axis data
    
    # a set of y data
    y_datas = []
    for i in range(1,101):
        y_data = i * np.sin(np.linspace(0, 2 * np.pi, 100))
        y_datas.append(y_data)


    # ============================= range_mode = auto =============================
    animator_auto_mean = CumAnimate(x_data, y_datas, fps=20, range_mode='auto', mode="mean")
    animator_auto_mean.ax.set_xlabel("Time")
    animator_auto_mean.ax.set_ylabel("Mean of Sine Waves")
    # animator_auto_mean.show()
    animator_auto_mean.save("CumAnimate_range_auto_mode_mean.gif")

    animator_auto_sum = CumAnimate(x_data, y_datas, fps=20, range_mode='auto', mode="sum")
    animator_auto_sum.ax.set_xlabel("Time")
    animator_auto_sum.ax.set_ylabel("Sum of Sine Waves")
    # animator_auto_sum.show()
    animator_auto_sum.save("CumAnimate_range_auto_mode_sum.gif")

    animator_auto = CumAnimate(x_data, y_datas, fps=20, range_mode='auto', mode="sequence")
    animator_auto.ax.set_xlabel("Time")
    animator_auto.ax.set_ylabel("Sine Waves")
    # animator_auto.show()
    animator_auto.save("CumAnimate_range_auto_mode_sequence.gif")


    # ============================ range_mode = fix =============================
    animator_fix_mean = CumAnimate(x_data, y_datas, fps=20, range_mode='fix', mode="mean")
    animator_fix_mean.ax.set_xlabel("Time")
    animator_fix_mean.ax.set_ylabel("Mean of Sine Waves")
    # animator_fix_mean.show()
    animator_fix_mean.save("CumAnimate_range_fix_mode_mean.gif")
    
    animator_fix_sum = CumAnimate(x_data, y_datas, fps=20, range_mode='fix', mode="sum")
    animator_fix_sum.ax.set_xlabel("Time")
    animator_fix_sum.ax.set_ylabel("Sum of Sine Waves")
    # animator_fix_sum.show()
    animator_fix_sum.save("CumAnimate_range_fix_mode_sum.gif")

    animator_fix = CumAnimate(x_data, y_datas, fps=20, range_mode='fix', mode="sequence")
    animator_fix.ax.set_xlabel("Time")
    animator_fix.ax.set_ylabel("Sine Waves")
    # animator_fix.show()
    animator_fix.save("CumAnimate_range_fix_mode_sequence.gif")


