"""Final project functions"""
import numpy as np
from netCDF4 import Dataset

class WindData:
    def __init__(self, data_file_path):
        self.data_file_path = data_file_path
    
    def get_rootgrp(self):
        """
        Opens a NetCDF file and returns the root group.
        This method uses the `Dataset` class from the `netCDF4` library to open
        the NetCDF file specified by `self.data_file_path` in read-only mode.
        Returns:
            netCDF4.Dataset: The root group of the NetCDF file.
        Raises:
            FileNotFoundError: If the file specified by `self.data_file_path` does not exist.
            OSError: If there is an issue opening the file.
        """
        
        rootgrp = Dataset(self.data_file_path, "r", format="NETCDF4")
        return rootgrp
    
    def get_latitude(self):
        """
        Retrieves the latitude values from the NetCDF dataset.
        This method accesses the 'latitude' variable in the NetCDF file and 
        extracts its values into a NumPy array.
        Returns:
            numpy.ndarray: A NumPy array containing the latitude values.
        """

        rootgrp = self.get_rootgrp()
        latitudes = rootgrp.variables['latitude']
        latitudes_out = np.zeros(len(latitudes))
        for i in range(len(latitudes)):
            latitudes_out[i] = latitudes[i]

        return latitudes_out