"""Final project functions"""
def test(x):
    return x

from netCDF4 import Dataset

class WindData:
    def __init__(self, data_file_path):
        self.data_file_path = data_file_path
    
    def get_rootgrp(self):
        rootgrp = Dataset(self.data_file_path, "r", format="NETCDF4")
        return rootgrp
    
    def get_latitude(self):
        rootgrp = self.get_rootgrp()
        latitude = rootgrp.variables['latitude']
        return latitude