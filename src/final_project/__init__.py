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
    
    def get_longitudes(self):
        """
        Retrieves the latitude values from the NetCDF dataset.
        This method accesses the 'latitude' variable in the NetCDF file and 
        extracts its values into a NumPy array.
        Returns:
            numpy.ndarray: A NumPy array containing the latitude values.
        """

        rootgrp = self.get_rootgrp()
        longitudes = rootgrp.variables['longitude']
        longitudes_out = np.zeros(len(longitudes))
        for i in range(len(longitudes)):
            longitudes_out[i] = longitudes[i]

        return longitudes_out
    

    
    def get_components_of_wind(self, component_name):
        """
        Retrieves the wind components from the NetCDF dataset.
        This method accesses the specified wind component variables in the NetCDF file
        and extracts their values into a dictionary.
        The keys of the dictionary are the names of the wind components, and the values for each loacation

        Input:
            component_name (list): A list of strings representing the names of the wind components to retrieve.
        Returns:
            dict: A dictionary containing the wind components for each location.
                  The keys are the names of the wind components, and the values are NumPy arrays of shape (4, n).
        """

        wind_data = {}   # Dictionary to store wind data for each component

        latitudes = self.get_latitude()    #Extract latitudes from the NetCDF file
        longitudes = self.get_longitudes() #Extract longitudes from the NetCDF file

        rootgrp = self.get_rootgrp()
        
        wind_data['locations'] = np.array([[latitudes[0], longitudes[0]],
                                      [latitudes[1], longitudes[0]],
                                      [latitudes[0], longitudes[1]],
                                      [latitudes[1], longitudes[1]]])  #This is the location of the wind data
        
        for i in range(len(component_name)):   # Loop through each component name
            # Extract the data for the current component
            data = rootgrp.variables[component_name[i]]
            wind_data[component_name[i]] = np.array([data[:, 0, 0], data[:, 1, 0], data[:, 0, 1], data[:, 1, 1]])

        return wind_data

