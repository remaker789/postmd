import numpy as np
from .animate import Animate


class CumAnimate(Animate):
    def __init__(self, *args, mode: str = None, fps: int = 20, dpi: int = 72, range_mode: str = 'auto'):
        """The cumulative animation of a set of data, which is usually for show the change of data with number of experiments or cumulative average.

        Args:
            *args: a set of y data ``(ydata1, ydata2, ...)``, or ``x, (ydata1, ydata2, ...)``. Each `ydata` should be a 1D array.
            mode (str, optional): ``None``, ``"sum"`` or ``"mean"``. Defaults to `None`.

                - If ``mode=None``, ``ydata1``, ``ydata2``, ... will display in order as animation
                - If ``mode="sum"``, the cumulative summation of ``(ydata1, ydata2, ...)`` will display in order as animation. The n-th frame is ``ydata1+ydata2+...+ydatan``
                - If ``mode="mean"``, the cumulative average of ``(ydata1, ydata2, ...)`` will display in order as animation. The n-th frame is ``(ydata1+ydata2+...+ydatan)/n``

            fps (int, optional): Frames per second (FPS) setting for animation. Defaults to ``20``.
            dpi (int, optional): Dots per inch (DPI) resolution for animation. Defaults to ``72``.
            range_mode (str, optional): Determines the behavior of x and y limits in the animation. Defaults to ``'auto'``.

                - If set to ``'auto'``, the x and y limits will automatically adjust to the data range.
                - If set to ``'fix'``, the x and y limits will be fixed to a specific range. The fixed range is from the min to max of data plus 5% margin.

        """
        super().__init__(*args, fps=fps, dpi=dpi)
        self.range_mode = range_mode
        self.mode = mode
        self._post_process()
        if range_mode == 'fix':
            self._range_mode_fix()


    def _post_process(self):
        if self.mode is None:
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

