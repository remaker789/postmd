# import sys,os
# from pathlib import Path
# sys.path.insert(0,Path(__file__).parent.parent.as_posix())

# print(sys.path)
from postmd.log import LogFile

path = r"D:\Onedrive\Github\MyRepos\PostMD\tests\data\log.lammps"
logfile=LogFile()
logfile.set_path(path)
logfile.get_lammps_version()
logfile.extract_thermodata()