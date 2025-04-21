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
    
    def compute_wind_speed_direction(self, location, component_name):
        
        wind_data = self.get_components_of_wind(component_name)
        location_index = np.where((wind_data['locations'] == location).all(axis=1))[0][0]

        speed = {}
        direction = {}
        speed_10m = np.zeros(len(wind_data['u10'][0]))
        speed_100m = np.zeros(len(wind_data['u100'][0]))
        direction_10m = np.zeros(len(wind_data['u10'][0]))
        direction_100m = np.zeros(len(wind_data['u100'][0]))
        for i in range(len(wind_data[component_name[0]][0])):
            x_value_10m = wind_data['u10'][location_index][i]
            y_value_10m = wind_data['v10'][location_index][i]
            x_value_100m = wind_data['u100'][location_index][i]
            y_value_100m = wind_data['v100'][location_index][i]

            complex_vector_10m = x_value_10m + y_value_10m*1j
            complex_vector_100m = x_value_100m + y_value_100m*1j

            speed_10m[i] = np.abs(complex_vector_10m)
            speed_100m[i] = np.abs(complex_vector_100m)

            direction_math_10m = np.rad2deg(np.angle(complex_vector_10m))
            direction_math_100m = np.rad2deg(np.angle(complex_vector_100m))

            direction_10m[i] = (90 - direction_math_10m) % 360
            direction_100m[i] = (90 - direction_math_100m) % 360
        
        speed['10m'] = speed_10m
        speed['100m'] = speed_100m
        direction['10m'] = direction_10m
        direction['100m'] = direction_100m
        return speed, direction

    def interpolate_at_loc(self, speed_locs, direction_locs, height, loc_lon, loc_lat):
        latitudes = self.get_latitude()    #Extract latitudes from the NetCDF file
        longitudes = self.get_longitudes() #Extract longitudes from the NetCDF file
       
        x1, x2 = longitudes[1], longitudes[0]
        y1, y2 = latitudes[0], latitudes[1]

        x, y = loc_lon, loc_lat

        # Bilinear interpolation estimates wind at Horns Rev by weighting values from the four surrounding points based on Horns Rev’s position within the grid
        wx1 = (x2 - x) / (x2 - x1)    # weight for longitude closer to left
        wx2 = (x - x1) / (x2 - x1)    # weight for longitude closer to right
        wy1 = (y2 - y) / (y2 - y1)    # weight for latitude closer to bottom
        wy2 = (y - y1) / (y2 - y1)    # weight for latitude closer to top

        interpolated_speed = (
        speed_locs[3][height] * wx1 * wy1 +  # L4
        speed_locs[0][height] * wx2 * wy1 +  # L1
        speed_locs[2][height] * wx1 * wy2 +  # L3
        speed_locs[1][height] * wx2 * wy2    # L2
         )

        interpolated_direction = (
        direction_locs[3][height] * wx1 * wy1 +
        direction_locs[0][height] * wx2 * wy1 +
        direction_locs[2][height] * wx1 * wy2 +
        direction_locs[1][height] * wx2 * wy2
         )

        return interpolated_speed, interpolated_direction