animate package
================


We built ``postmd.animate`` subpackage to help to animate the data. This subpackage 
is based on the ``matplotlib.animation.FuncAnimation`` class. 

We provided two methods to preview and save the animation.
- ``show()``: preview the animation.
- ``save()``: save the animation to a file(``git``, ``mp4``, or ``avi``).


We provide two derived classes: ``AppendAnimate`` and ``CumAnimate``.

.. toctree::
   :maxdepth: 1
   :titlesonly:
   
   cum_animate
   append_animate


.. note::
   We have initialized ``self.fig, self.ax = plt.subplots()`` in the base class ``Animate``.
   Therefore, you can use ``self.fig`` and ``self.ax`` to customize the figure of animation.
   See :ref:`append_animate_example` for more details.



Customize
--------------------------

The Base class is ``Animate`` class. The following method may be overrided in the derived class:

- ``_init(self)``: initialize the data in ``self.ax`` 
- ``_update_data(self,n)``: tell ``FuncAnimation`` how to update the data in ``self.ax``. ``n`` means the n-th frame in the animation. Should return ``(xdata,ydata)``.
- ``_update_frame(self,n)``: tell ``FuncAnimation`` how to update the frame in the animation. ``n`` means the n-th frame in the animation. **Usually, there is no need to override this method.**

