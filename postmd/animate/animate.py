import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

"""

===================
test
===================

"""


class Animate:
    def __init__(self, *args, fps=20, dpi=72):
        if len(args) == 1:  # Only y-data provided
            self.y = np.array(args[0])
            self.time_steps = len(self.y)
            self.x = np.arange(self.time_steps)
        elif len(args) == 2:  # Both x and y data provided
            self.x, self.y = map(np.array, args)
            self.time_steps = len(self.y)
            assert len(self.x) == self.time_steps, "The length of x must match the length of y"
        else:
            raise ValueError("Invalid input. Provide either y-axis data or both x and y-axis data.")

        self.fps = fps
        self.dpi = dpi
        
        self.fig, self.ax = plt.subplots()
        self.line, = self.ax.plot([], [], lw=2)



    def _init(self):
        """init the animation
        """        
        pass


    def _update_data(self, n):
        pass

    def _update_frame(self, n):
        _xdata, _ydata = self._update_data(n)
        self.line.set_data(_xdata, _ydata)

        if self.range_mode == 'auto':
            # Dynamically update the y-axis limits
            self.ax.relim()
            self.ax.autoscale_view()

        return self.line,


    def show(self):
        """show the animation of data.
        """        
        ani = FuncAnimation(self.fig, self._update_frame, frames=range(self.time_steps),
                            init_func=self._init, blit=False, interval=1000/self.fps, repeat=False)
        plt.show()

    def save(self, filename, output_format=None):
        """save the animation of data to file.

        Args:
            filename (str): file name to output
            output_format (str, optional): the output format. Now only support ``gif``, ``mp4``, and ``avi``. Defaults to ``None``, meaning the output format is inferred from the filename extension.

        """        
        self.output_format=output_format
        if self.output_format is None:
            self.output_format = os.path.splitext(filename)[-1][1:]  # file extension

        ani = FuncAnimation(self.fig, self._update_frame, frames=range(self.time_steps),
                            init_func=self._init, blit=True, interval=1000/self.fps, repeat=False)

        if self.output_format.lower() == 'gif':
            ani.save(filename, writer='pillow', fps=self.fps)
        elif self.output_format.lower() in ['mp4', 'avi']:
            from matplotlib.animation import FFMpegWriter
            writer = FFMpegWriter(fps=self.fps)
            ani.save(filename, writer=writer)
        else:
            raise ValueError(f"Unsupported output format: {self.output_format}")

        print(f"Animation saved to file: {filename}")


