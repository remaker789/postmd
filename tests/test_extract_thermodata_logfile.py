from PostMD.FileOperation import *

path = r"D:\Onedrive\Github\MyRepos\PostMD\tests\data\log.lammps"
logfile=LogFile()
logfile.set_path(path)
logfile.get_lammps_version()
logfile.extract_thermodata()