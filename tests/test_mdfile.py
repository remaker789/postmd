import os 

os.chdir(os.path.dirname(__file__))

from postmd.system import MDFile

path = r"data/ave_chunk_bin1d.dat"
data_bin1d = MDFile(path)
data_bin1d._set_path(path)
prop = data_bin1d._get_header(2)
print(prop)
