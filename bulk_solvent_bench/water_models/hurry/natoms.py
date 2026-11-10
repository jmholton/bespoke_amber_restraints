# atom count of each Amber NetCDF trajectory given: "file natoms" per line.  usage: phenix.python natoms.py a.nc ...
import sys
from scipy.io import netcdf_file
for f in sys.argv[1:]:
    print(f, netcdf_file(f, "r", mmap=False).dimensions["atom"])
