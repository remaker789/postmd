import numpy as np
from .animate import Animate


class AppendAnimate(Animate): 
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
