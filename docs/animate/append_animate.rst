AppendAnimate class
=====================

.. autoclass:: postmd.animate.AppendAnimate
    :members:
    :undoc-members:
    :show-inheritance:



.. _append_animate_example:

Example
-------------


.. code-block:: python

    import numpy as np
    from postmd.animate import AppendAnimate

    if __name__ == '__main__':
        y_data = np.sin(np.linspace(0, 2 * np.pi, 100))  # Assume this is the sine wave data to be animated
        x_data = np.linspace(0, 10, 100)  # Corresponding x-axis data

        animator_auto = AppendAnimate(x_data, y_data, fps=20, range_mode='auto')
        # add label of x and y axis.
        animator_auto.ax.set_xlabel("Time")
        animator_auto.ax.set_ylabel("Sine Wave")
        # animator_auto.show()
        animator_auto.save("animate_sine_wave_auto.gif")



        animator_fix = AppendAnimate(x_data, y_data, fps=20, range_mode='fix')
        animator_fix.ax.set_xlabel("Time")
        animator_fix.ax.set_ylabel("Sine Wave")
        # animator_fix.show()
        animator_fix.save("animate_sine_wave_fix.gif")


.. figure:: ../images/animate_sine_wave_auto.gif
    :width: 300px
    :align: left

    animate_sine_wave_auto.gif

.. figure:: ../images/animate_sine_wave_fix.gif
    :width: 300px
    :align: right

    animate_sine_wave_fix.gif