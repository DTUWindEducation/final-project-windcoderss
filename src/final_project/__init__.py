"""Final project functions"""
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from scipy.stats import weibull_min
from windrose import WindroseAxes
from scipy.integrate import quad
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

        rootgrp.close()
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
        rootgrp.close()
        return longitudes_out
    
    def get_time(self):
        """
        Get the time stamps of the data
        """
        rootgrp = self.get_rootgrp()
        # Extract the timestamps for the hours:
        time_dim = rootgrp.dimensions['time']
        num_hours = len(time_dim)
        hours = np.arange(num_hours)
        rootgrp.close()

        return hours

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
        rootgrp.close()
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

            direction_10m[i] = (270 - direction_math_10m) % 360
            direction_100m[i] = (270 - direction_math_100m) % 360
        
        speed['10m'] = speed_10m
        speed['100m'] = speed_100m
        direction['10m'] = direction_10m
        direction['100m'] = direction_100m
        return speed, direction

    def interpolate_at_loc(self, speed_locs, direction_locs, height, loc_lon, loc_lat):
        latitudes = self.get_latitude()    #Extract latitudes from the NetCDF file
        longitudes = self.get_longitudes() #Extract longitudes from the NetCDF file
       
        x1, x2 = longitudes[0], longitudes[1]
        y1, y2 = latitudes[1], latitudes[0]

        y, x = loc_lat, loc_lon

        # Bilinear interpolation estimates wind at Horns Rev by weighting values from the four surrounding points based on Horns Rev’s position within the grid
        wx1 = (x2 - x) / (x2 - x1)    # closer to lon x1 (west)
        wx2 = (x - x1) / (x2 - x1)    # closer to lon x2 (east)
        wy1 = (y2 - y) / (y2 - y1)    # closer to lat y1 (south)
        wy2 = (y - y1) / (y2 - y1)    # closer to lat y2 (north)

        interpolated_speed = (
        speed_locs[0][height] * wx1 * wy1 +  # L4
        speed_locs[1][height] * wx2 * wy1 +  # L1
        speed_locs[2][height] * wx1 * wy2 +  # L3
        speed_locs[3][height] * wx2 * wy2    # L2
         )

        interpolated_direction = (
        direction_locs[0][height] * wx1 * wy1 +
        direction_locs[1][height] * wx2 * wy1 +
        direction_locs[2][height] * wx1 * wy2 +
        direction_locs[3][height] * wx2 * wy2
         )

        return interpolated_speed, interpolated_direction
    

# Functions that is not part of the class
def windspeed_at_height(wind_speed, height, alpha):
    """
    Calculate wind speed at a given height using the power law.
    """
    if height > 55:
        u = wind_speed['100m']*(height/100)**alpha
    else:
        u = wind_speed['10m']*(height/10)**alpha

    return u

def fit_and_plot_weibull(wind_speeds):
    """
    Fits a Weibull distribution to wind speed data and plots the PDF.
    
    Parameters:
    wind_speeds (array-like): Array of wind speed values [m/s].
    
    Returns:
    shape (float): Weibull shape parameter (k).
    scale (float): Weibull scale parameter (A).
    """
    wind_speeds = np.array(wind_speeds)

    # Fit the Weibull distribution (fix location to 0)
    shape, loc, scale = weibull_min.fit(wind_speeds, floc=0)

    # Generate x values for plotting
    x = np.linspace(0, max(wind_speeds) + 2, 100)
    pdf = weibull_min.pdf(x, shape, loc, scale)

    # Plot
    plt.figure(figsize=(8, 5))
    plt.hist(wind_speeds, bins=100, density=True, alpha=0.6, color='skyblue', label='Histogram')
    plt.plot(x, pdf, 'r-', lw=2, label=f'Weibull PDF\nk={shape:.2f}, A={scale:.2f}')
    plt.xlabel('Wind Speed [m/s]')
    plt.ylabel('Probability Density')
    plt.title('Weibull Distribution Fit to Wind Speed Data')
    plt.legend()
    plt.grid(True)
    plt.show(block=False)

    return shape, scale

def plot_wind_rose(wind_speed, wind_direction):
    
    ax = WindroseAxes.from_ax()
    ax.bar(wind_direction, wind_speed, normed=True, opening=0.8, edgecolor="white")
    ax.set_legend()

    return
class WindTurbine:
    def __init__(self, windturbine_file_path):
        self.windturbine_file_path = windturbine_file_path

    def read_data(self):
        windturbine_data = np.loadtxt(self.windturbine_file_path, skiprows=1, delimiter=",")

        return windturbine_data
    
    def get_power(self, windspeed):
        windturbine_data = self.read_data()
        P = np.interp(windspeed, windturbine_data[:, 0], windturbine_data[:, 1])

        return P
    
    def get_AEP(self, Weibull_shape, Weibull_scale, eta):
        windturbine_data = self.read_data()
        u_in = windturbine_data[1, 0]
        u_out = windturbine_data[-1, 0]
        k = Weibull_shape
        A = Weibull_scale
        
        def f(u):
            return (k/A) *(u/A)**(k-1)*np.exp(-(u/A)**k)
        
        def integrand(u):
            return self.get_power(u) * f(u)
        
        integral, _ = quad(integrand, u_in, u_out)
        
        AEP = eta * 8760 * integral
        return AEP
    
#Extra function no 1
    def plot_power_output(self, windspeed, time):
        power_output = self.get_power(windspeed)
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(time, power_output, label='Power Output')
        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('Power Output (kW)')
        ax.set_title('Power Output of Wind Turbine')
        ax.legend()
        ax.grid(True)
        plt.show(block=False)
        return
    
 #Extra function no 2   
    def plot_power_duration_curve(self, windspeed, time):
        power_output = self.get_power(windspeed)
        # Sort the power output in descending order
        sorted_power = np.sort(power_output)[::-1]
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(time, sorted_power, label='Power Output Sorted')
        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('Power Output (kW)')
        ax.set_title('Duration Curve for Power Output of Wind Turbine')
        ax.legend()
        ax.grid(True)
        plt.show(block=False)
        return