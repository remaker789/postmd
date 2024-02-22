================
GreenKubo
================
The Green-Kubo formula usually have the following form including an auto-correlation function(ACF) of a certain property,

.. math::
    C \int_{0}^{\infty}\left\langle p(0) \cdot p(t) \right\rangle \mathrm{d} t

where :math:`C` is a prefactor parameter, :math:`\left\langle p(0) \cdot p(t) \right\rangle` is ACF of property :math:`p`.

For example, the Green-Kubo formula(vscf) for calculate the diffusion coefficient is,

.. math::
    D=\frac{1}{3} \int_{0}^{\infty}\left\langle\mathbf{v}^{c}(t) \cdot \mathbf{v}^{c}(0)\right\rangle \mathrm{d} t

**The calculation of ACF is central in the Green-Kubo formula**. Here we devide Green-Kubo data into two types: ``raw`` and ``acf``:

- ``raw``: just the property :math:`p`
- ``acf``: ACF :math:`\left\langle p(0) \cdot p(t) \right\rangle`, usually generated from command ``fix ave/correlate type auto overwriting``, see `fix ave/correlate <https://docs.lammps.org/fix_ave_correlate.html>`_ command

raw-type data
=============

.. code-block:: python
    
    ## ================ For raw data ===================
    import numpy as np
    import matplotlib.pyplot as plt 
    import scipy.constants as C
    from postmd.avetime import GreenKubo
    e=C.e

    unit_trans = e/1e-10  # eV/Ang to N
    gk=GreenKubo()
    gk.read_file(path=raw_data_path, name_line=2)
    gk.cal_acf(data_type="raw",col=3,nlag=2000,unit_trans=unit_trans)
    plt.figure()
    plt.plot(gk.nlag,gk.acf)
    plt.title("computed from raw data")
    plt.xlabel("nlag")
    plt.ylabel("acf")
    plt.legend()


    A = 2*np.pi*2.335e-10*200e-10
    kB = C.Boltzmann
    T = 298
    plt.figure()
    gk.integrate_acf(method="simpson")
    plt.plot(gk.nlag, gk.int_acf/(A*kB*T), label="simpson method")
    gk.integrate_acf(method="trap")
    plt.plot(gk.nlag, gk.int_acf/(A*kB*T), label="trap method")
    # print( (gk.int_acf/(A*kB*T)).max())
    plt.xlabel("nlag")
    plt.ylabel("integration of acf")
    plt.title("computed from raw data")
    plt.legend()

acf-type data
=============

.. code-block:: python

    ## ================ For autocorrelation data ==================
    import numpy as np
    import matplotlib.pyplot as plt 
    import scipy.constants as C
    from postmd.avetime import GreenKubo
    e=C.e

    unit_trans = e/1e-10  # friction force的单位是eV/Ang
    gk=GreenKubo()
    gk.read_file(path=acf_data_path, name_line=3, skiprows=4)
    gk.cal_acf(data_type="acf",col=5,nlag_col=1,unit_trans=unit_trans)
    plt.figure()
    plt.plot(gk.nlag,gk.acf)
    plt.title("computed from acf data")
    plt.xlabel("nlag")
    plt.ylabel("acf")
    plt.legend()


    A = 2*np.pi*2.335e-10*200e-10
    kB = C.Boltzmann
    T = 298
    plt.figure()
    gk.integrate_acf(method="simpson")
    plt.plot(gk.nlag, gk.int_acf/(A*kB*T), label="simpson method")
    gk.integrate_acf(method="trap")
    plt.plot(gk.nlag, gk.int_acf/(A*kB*T), label="trap method")
    # print( (gk.int_acf/(A*kB*T)).max())
    plt.xlabel("nlag")
    plt.ylabel("integration of acf")
    plt.title("computed from acf data")
    plt.legend()