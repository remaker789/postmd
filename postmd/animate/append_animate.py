import numpy as np
from .animate import Animate


class AppendAnimate(Animate): 
    def __init__(self, *args, fps=20, dpi=72, range_mode='auto'):
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
