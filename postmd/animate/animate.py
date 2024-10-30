import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class BaseAnimate:
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


class AppendAnimate(BaseAnimate): 
    def __init__(self, *args, fps=20, dpi=72, range_mode='auto'):
        """animate data by appending data to the figure.

        Args:
            *args: a set of y data ``(ydata1, ydata2, ...)``, or ``x, (ydata1, ydata2, ...)``. Each `ydata` should be a 1D array.
            fps (int, optional): Frames per second (FPS) setting for animation. Defaults to ``20``.
            dpi (int, optional): Dots per inch (DPI) resolution for animation. Defaults to ``72``.
            range_mode (str, optional): Determines the behavior of x and y limits in the animation. Defaults to ``'auto'``.
                - If set to ``'auto'``, the x and y limits will automatically adjust to the data range.
                - If set to ``'fix'``, the x and y limits will be fixed to a specific range. The fixed range is from the min to max of data plus 5% margin.
        """        
        super().__init__(*args, fps=fps, dpi=dpi)
        self.range_mode = range_mode
        if range_mode == 'fix':
            self._range_mode_fix()


    def _range_mode_fix(self):
        # Set fixed x-axis limits (assuming the provided x data already represents time steps appropriately)
        x_min, x_max = np.min(self.x), np.max(self.x)
        padding = (x_max - x_min) * 0.05  # Add 5% margins on both sides
        self.ax.set_xlim(x_min - padding, x_max + padding)
        # Calculate and set fixed y-axis limits
        y_min, y_max = np.min(self.y), np.max(self.y)
        padding = (y_max - y_min) * 0.05  # Add 5% margins on both sides
        self.ax.set_ylim(y_min - padding, y_max + padding)


    def _init(self):
        self.line.set_data([], [])
        return self.line,

    def _update_data(self, n):
        _xdata=self.x[:n+1]
        _ydata=self.y[:n+1]
        return _xdata, _ydata
    
    def show(self):
        super().show()

    def save(self, filename, output_format=None):
        super().save(filename, output_format)



class CumAnimate(BaseAnimate):
    def __init__(self, *args, fps:int=20, dpi:int=72, range_mode:str='auto', mode:str='mean'):
        """The cumulative animation of a set of data, which is usually for show the change of data with number of experiments or cumulative average.

        Args:
            *args: a set of y data ``(ydata1, ydata2, ...)``, or ``x, (ydata1, ydata2, ...)``. Each `ydata` should be a 1D array.
            fps (int, optional): Frames per second (FPS) setting for animation. Defaults to ``20``.
            dpi (int, optional): Dots per inch (DPI) resolution for animation. Defaults to ``72``.
            range_mode (str, optional): Determines the behavior of x and y limits in the animation. Defaults to ``'auto'``.
                - If set to ``'auto'``, the x and y limits will automatically adjust to the data range.
                - If set to ``'fix'``, the x and y limits will be fixed to a specific range. The fixed range is from the min to max of data plus 5% margin.
            mode (str, optional): ``"mean"``, ``"sum"`` or ``"mean"``. Defaults to `None`.
                - If ``mode="mean"``, the cumulative average of ``(ydata1, ydata2, ...)`` will display in order as animation. The n-th frame is ``(ydata1+ydata2+...+ydatan)/n``
                - If ``mode="sum"``, the cumulative summation of ``(ydata1, ydata2, ...)`` will display in order as animation. The n-th frame is ``ydata1+ydata2+...+ydatan``
                - If ``mode="sequence"``, ``ydata1``, ``ydata2``, ... will display in sequence as animation
                

        """
        super().__init__(*args, fps=fps, dpi=dpi)
        self.range_mode = range_mode
        self.mode = mode
        self._post_process()
        if range_mode == 'fix':
            self._range_mode_fix()


    def _post_process(self):
        if self.mode == "sequence":
            self.ydata = self.y
        elif self.mode == 'sum':
            self.ydata = np.cumsum(self.y, axis=0)
        elif self.mode == 'mean':
            self.ydata = np.cumsum(self.y, axis=0) / (np.arange(1, self.time_steps + 1).reshape(-1,1))


    def _range_mode_fix(self):
        x_min, x_max = np.min(self.x), np.max(self.x)
        padding = (x_max - x_min) * 0.05  # Add 5% margins on both sides of x axis
        self.ax.set_xlim(x_min - padding, x_max + padding)

        y_min, y_max = np.min(self.ydata), np.max(self.ydata)
        padding = (y_max - y_min) * 0.05  # Add 5% margins on both sides of y axis
        self.ax.set_ylim(y_min - padding, y_max + padding)
        
        
    def _init(self):
        self.line.set_data([], [])
        return self.line,



    def _update_data(self, n):
        _ydata = self.ydata[n,:]
        return self.x, _ydata
    

    def show(self):
        super().show()

    def save(self, filename, output_format=None):
        super().save(filename, output_format)

