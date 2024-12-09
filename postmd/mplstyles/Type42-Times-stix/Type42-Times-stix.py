import os
import numpy as np
import matplotlib.pyplot as plt

plt.style.use("Type42-Times-stix\Type42-Times-stix.mplstyle")
x = np.arange(1,11,1)
y = x**2
plt.plot(x,y, '-o', label="legend")
plt.text(2,14, '5aKPS.')
plt.text(2,50, r'$x^2+y^2=1$')
plt.xlabel(r"$E$ [V/$\mathdefault{\AA}$]")
plt.ylabel(r"$I_c$ [A]")
plt.title('title')
plt.legend()
plt.savefig(os.path.join(os.path.dirname(__file__), 'Type42-Times-stix_preview.pdf'))
plt.show()
