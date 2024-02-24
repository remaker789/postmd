CumAnimate class
===================

.. autoclass:: postmd.animate.CumAnimate
    :members:
    :undoc-members:
    :show-inheritance:

Example:

.. code-block:: python

    from postmd.animate import CumAnimate
    import numpy as np

    if __name__ == '__main__':
        y_datas = []
        for i in range(100):
            y_data = np.sin(np.linspace(0, 2 * np.pi, 100)) + i/100
            y_datas.append(y_data)

        x_data = np.linspace(0, 10, 100)  # Corresponding x-axis data

        animator_auto = CumAnimate(y_datas, fps=20, range_mode='auto')
        # animator_auto.show()
        animator_auto.save("cum_test1.mp4")

        animator_auto_sum = CumAnimate(y_datas, fps=20, range_mode='auto', mode="sum")
        # animator_auto.show()
        animator_auto_sum.save("cum_sum_test1.mp4")

        animator_auto_mean = CumAnimate(y_datas, fps=20, range_mode='auto', mode="mean")
        # animator_auto.show()
        animator_auto_mean.save("cum_mean_test1.mp4")


        animator_fix = CumAnimate(y_datas, fps=20, range_mode='fix')
        # animator_fix.show()
        animator_fix.save("cum_test2.mp4")

        animator_fix_sum = CumAnimate(y_datas, fps=20, range_mode='fix', mode="sum")
        # animator_fix.show()
        animator_fix_sum.save("cum_sum_test2.mp4")

        animator_fix_mean = CumAnimate(y_datas, fps=20, range_mode='fix', mode="mean")
        # animator_fix.show()
        animator_fix_mean.save("cum_mean_test2.mp4")





