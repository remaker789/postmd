AveTime class
================

The base class ``AveTime`` defined some method including read the data


Examples:

>>> from postmd.avetime import AveTime
>>> myavetime = AveTime()
---------------- System Properties --------------
Temperature:    298.0 K
Timestep:       1.0 fs
-------------------------------------------------
>>> myavetime.set_path(r"./tests/data/ave_time.dat") 
You are processing the file: './tests/data/ave_time.dat'
>>> df = myavetime.read_file()
You are processing the file: './tests/data/ave_time.dat'
>>> df.head()
TimeStep      c_jVz  c_jQzIon   c_jQzNa  c_jVzIon  c_jVzAvg  c_jVzAvg_spce
0     10000  17021.600  111.1920  111.1920  111.1920  71.21990        71.3517
1     20000   -527.274   14.2479   14.2479   14.2479  -2.20617        -2.2849
2     30000   3750.990   42.8537   42.8537   42.8537  15.69450        15.6462
3     40000   4519.150   95.6945   95.6945   95.6945  18.90860        18.6644
4     50000   9854.170   30.8659   30.8659   30.8659  41.23080        41.4485