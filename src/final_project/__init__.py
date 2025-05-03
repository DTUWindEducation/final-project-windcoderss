"""Final project functions"""
import numpy as np
from netCDF4 import Dataset, num2date
import matplotlib.pyplot as plt
from scipy.stats import weibull_min
from windrose import WindroseAxes
from scipy.integrate import quad


class WindData:
    """WindData is a class for handling and analyzing wind data
    from NetCDF files.
    It provides methods to extract latitude, longitude, time, wind components,
    compute wind speed and direction, and perform bilinear interpolation for
    specific locations.
    Attributes:
        data_file_path (str): Path to the NetCDF file containing wind data.
    Methods:
        get_rootgrp(): Opens the NetCDF file and returns the root group.
        get_latitude(): Retrieves latitude values from the dataset.
        get_longitudes(): Retrieves longitude values from the dataset.
        get_time(): Extracts time-related information from the dataset.
        get_components_of_wind(component_name):
        Retrieves specified wind components.
        compute_wind_speed_direction(location, component_name):
        Computes wind speed and direction at 10m and 100m heights for
        a given location.
        interpolate_at_loc(speed_locs, direction_locs, loc_lon, loc_lat):
        Performs bilinear interpolation to estimate wind speed
        and direction at a specific location.
    """
    def __init__(self, data_file_path):
        self.data_file_path = data_file_path

    def get_rootgrp(self):
        """
        Opens a NetCDF file and returns the root group.
        Returns:
            netCDF4.Dataset: The root group of
                             the NetCDF file opened in read mode.
        """
        # Open the NetCDF file in read mode
        rootgrp = Dataset(self.data_file_path, "r", format="NETCDF4")
        return rootgrp

    def get_latitude(self):
        """
        Retrieves the latitude values from a NetCDF dataset.
        Returns:
            numpy.ndarray: An array containing the latitude values.
        """
        rootgrp = self.get_rootgrp()
        # Extract the latitude variable from the NetCDF file
        latitudes = rootgrp.variables['latitude']
        # Initialize an array to store the latitude values
        latitudes_out = np.zeros(len(latitudes))
        # Loop through the latitude values
        for i in range(len(latitudes)):
            latitudes_out[i] = latitudes[i]
        rootgrp.close()
        return latitudes_out

    def get_longitudes(self):
        """
        Retrieves the longitude values from the dataset.
        Returns:
            np.ndarray: Array of longitude values extracted from the dataset.
        """
        rootgrp = self.get_rootgrp()
        # Extract the longitude variable from the NetCDF file
        longitudes = rootgrp.variables['longitude']
        longitudes_out = np.zeros(len(longitudes))
        # Loop through the longitude values
        for i in range(len(longitudes)):
            longitudes_out[i] = longitudes[i]
        rootgrp.close()
        return longitudes_out

    def get_time(self):
        """
        Extracts time-related information from a NetCDF dataset.
        Returns:
            tuple: A tuple containing:
            - hours (np.ndarray): Array of hour indices.
            - years (np.ndarray): Array of corresponding years for each time.
        """
        rootgrp = self.get_rootgrp()
        # Extract the timestamps for the hours:
        time_dim = rootgrp.dimensions['time']
        num_hours = len(time_dim)
        hours = np.arange(num_hours)
        # Extract the years
        time_var = rootgrp.variables['time']
        dates = num2date(time_var[:], units=time_var.units,
                         calendar=getattr(time_var, 'calendar', 'standard'))
        years = [date.year for date in dates]
        years = np.array(years)  # Convert the list of years to a NumPy array
        # Close the rootgrp
        rootgrp.close()

        return hours, years

    def get_components_of_wind(self, component_name):
        """
        Retrieves the wind components from the NetCDF dataset.
        This method accesses the specified wind component variables
        in the NetCDF file and extracts their values into a dictionary.
        The keys of the dictionary are the names of the wind components,
        and the values for each loacation

        Input:
            component_name (list): A list of strings representing the names
            of the wind components to retrieve.
        Returns:
            dict: Dictionary containing the wind components for each location.
                  The keys are the names of the wind components,
                    and the values are NumPy arrays of shape (4, n).
        """

        wind_data = {}   # Dictionary to store wind data for each component

        latitudes = self.get_latitude()  # Extract lat from the NetCDF file
        longitudes = self.get_longitudes()  # Extract long from the NetCDF file
        rootgrp = self.get_rootgrp()
        wind_data['locations'] = np.array([[latitudes[0], longitudes[0]],
                                          [latitudes[1], longitudes[0]],
                                          [latitudes[0], longitudes[1]],
                                          [latitudes[1], longitudes[1]]])
        for i in range(len(component_name)):   # Loop through each component
            # Extract the data for the current component
            data = rootgrp.variables[component_name[i]]
            wind_data[component_name[i]] = np.array([data[:, 0, 0],
                                                    data[:, 1, 0],
                                                    data[:, 0, 1],
                                                    data[:, 1, 1]])
        rootgrp.close()
        return wind_data

    def compute_wind_speed_direction(self, location, component_name):
        """
        Computes wind speed and direction at 10m and 100m for a given location.
        Args:
            location (array-like):
                Coordinates of the location to compute wind data for.
            component_name (list):
                List of wind components (e.g., ['u10', 'v10', 'u100', 'v100']).
        Returns:
            tuple: A tuple containing:
                - speed (dict):
                Wind speed at 10m and 100m heights, with keys '10m' and '100m'.
                - direction (dict):
                Wind direction at 10m and 100m, with keys '10m' and '100m'.
        """
        wind_data = self.get_components_of_wind(component_name)
        location_index = np.where((wind_data['locations']
                                   == location).all(axis=1))[0][0]

        speed = {}
        direction = {}
        speed_10m = np.zeros(len(wind_data['u10'][0]))
        speed_100m = np.zeros(len(wind_data['u100'][0]))
        direction_10m = np.zeros(len(wind_data['u10'][0]))
        direction_100m = np.zeros(len(wind_data['u100'][0]))
        # Loop through each time step
        for i in range(len(wind_data[component_name[0]][0])):
            x_value_10m = wind_data['u10'][location_index][i]
            y_value_10m = wind_data['v10'][location_index][i]
            x_value_100m = wind_data['u100'][location_index][i]
            y_value_100m = wind_data['v100'][location_index][i]

            # Create a complex vector for 10m and 100m wind data
            complex_vector_10m = x_value_10m + y_value_10m*1j
            complex_vector_100m = x_value_100m + y_value_100m*1j

            # Calculate the speed and direction
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

    def interpolate_at_loc(self, speed_locs, direction_locs, loc_lon, loc_lat):
        """
        Perform interpolation to estimate wind speed and direction
         at a specific location.
        Args:
            speed_locs (list): Wind speed values at the
              four surrounding grid points [L4, L1, L3, L2].
            direction_locs (list): Wind direction values at the
              four surrounding grid points [L4, L1, L3, L2].
            loc_lon (float): Longitude of the target location.
            loc_lat (float): Latitude of the target location.
        Returns:
            tuple: Interpolated wind speed (float) and wind direction (float)
              at the target location.
        """

        latitudes = self.get_latitude()  # Extract latitudes from NetCDF file
        longitudes = self.get_longitudes()  # Extract longitudes from NetCDF

        x1, x2 = longitudes[0], longitudes[1]
        y1, y2 = latitudes[1], latitudes[0]

        y, x = loc_lat, loc_lon  # Coordinates of chosen location fx. Horns Rev

        # Interpolation estimates wind at Horns Rev by weighting
        #  values from the four surrounding points
        #  based on Horns Rev’s position within the grid
        wx1 = (x2 - x) / (x2 - x1)    # closer to lon x1 (west)
        wx2 = (x - x1) / (x2 - x1)    # closer to lon x2 (east)
        wy1 = (y2 - y) / (y2 - y1)    # closer to lat y1 (south)
        wy2 = (y - y1) / (y2 - y1)    # closer to lat y2 (north)

        interpolated_speed = (
                                speed_locs[0] * wx1 * wy1 +  # L4
                                speed_locs[1] * wx2 * wy1 +  # L1
                                speed_locs[2] * wx1 * wy2 +  # L3
                                speed_locs[3] * wx2 * wy2    # L2
                                )

        interpolated_direction = (
                                direction_locs[0] * wx1 * wy1 +
                                direction_locs[1] * wx2 * wy1 +
                                direction_locs[2] * wx1 * wy2 +
                                direction_locs[3] * wx2 * wy2
                                )

        return interpolated_speed, interpolated_direction


class WindTurbine:
    """
    WindTurbine is a class for handling and analyzing wind turbine data.
    It provides methods to read wind turbine data from a CSV file,  calculate
    power output based on wind speed, and compute Annual Energy Production
    using the Weibull distribution.
    Attributes:
        windturbine_file_path (str): Path to the wind turbine data file.
    Methods:
        read_data(): Reads wind turbine data from a CSV file.
        get_power(windspeed): Calculates power output based on wind speed.
        get_AEP(Weibull_shape, Weibull_scale, eta):
        Computes AEP using Weibull distribution.
        """
    def __init__(self, windturbine_file_path):
        # Path to the wind turbine data file
        self.windturbine_file_path = windturbine_file_path

    def read_data(self):
        """
        Reads wind turbine data from a CSV file.
        Returns:
            numpy.ndarray:
            A 2D array containing the wind turbine data loaded from the file.
        """
        # Load the wind turbine data from the CSV file
        windturbine_data = np.loadtxt(self.windturbine_file_path,
                                      skiprows=1, delimiter=",")

        return windturbine_data

    def get_power(self, windspeed):
        """
        Calculate the power output of a wind turbine 
        based on the given windspeed.
        Args:
            windspeed (float): The wind speed in m/s.
        Returns:
            float: The interpolated power output in watts.
        """

        windturbine_data = self.read_data()
        # Interpolate the power output based on the wind speed
        P = np.interp(windspeed, windturbine_data[:, 0],
                      windturbine_data[:, 1])
        return P

    def get_AEP(self, Weibull_shape, Weibull_scale, eta):
        """
        Calculate the Annual Energy Production (AEP) of a wind turbine.
        Parameters:
            Weibull_shape (float): Shape parameter (k).
            Weibull_scale (float): Scale parameter (A).
            eta (float): Efficiency factor of the wind turbine.
        Returns:
            float: The calculated AEP in kilowatt-hours (kWh).
        """

        windturbine_data = self.read_data()
        # Minimum wind speed for power generation
        u_in = windturbine_data[1, 0]
        # Maximum wind speed for power generation
        u_out = windturbine_data[-1, 0]
        k = Weibull_shape  # Weibull shape parameter for the selected year
        A = Weibull_scale  # Weibull scale parameter for the selected year

        # Weibull probability density function
        def f(u):
            return (k/A) * (u/A)**(k-1)*np.exp(-(u/A)**k)

        # Integrand function for AEP calculation
        def integrand(u):
            return self.get_power(u) * f(u)

        # Perform numerical integration to calculate AEP
        integral, _ = quad(integrand, u_in, u_out, limit=200,
                           epsabs=1e-5, epsrel=1e-5)

        # Calculate AEP in kWh (assuming 8760 hours in a year)
        AEP = eta * 8760 * integral
        return AEP

# Extra function no 1
    def compare_AEP(self, windspeed, eta, years, all_years):
        """
        Compare the Annual Energy Production (AEP) 
        of a wind turbine across different years.
        Parameters:
            windspeed (array-like): Measured wind speed data.
            eta (float): Efficiency factor of the wind turbine.
            years (array-like):
            Array of years corresponding to the wind speed data.
            all_years (array-like): List of years to compare AEP.
        Returns:
            None: Displays a plot of AEP (in MWh) for the given years.
        """

        AEP_all = []  # List to store AEP values for each year
        for i in range(len(all_years)):  # Loop through each year in all_years
            # Fit Weibull distribution for the selected year
            weibull_shape, weibull_scale = fit_and_plot_weibull(windspeed,
                                                                years,
                                                                all_years[i],
                                                                False)
            # Calculate AEP for the selected year
            AEP = self.get_AEP(weibull_shape, weibull_scale, eta)
            AEP_all.append(AEP)  # Append the AEP value to the list
        AEP_all = np.array(AEP_all)  # Convert the list to a NumPy array
        fig, ax = plt.subplots(figsize=(10, 6))  # Create a figure
        ax.plot(all_years, AEP_all/1000, label='AEP')
        ax.set_xlabel('Year')
        ax.set_ylabel('AEP (MWh)')
        ax.set_title('AEP of windturbine for different years')
        ax.legend()
        ax.grid(True)
        plt.show(block=False)
        return

# Extra function no 2
    def plot_power_duration_curve(self, windspeed, time, years, selected_year):
        """
        Plots the power duration curve for a wind turbine for a specific year.
        Args:
            windspeed (list or array): Wind speed data corresponding to time.
            time (list or array): Time data in hours.
            years (list or array): Year data corresponding to each time point.
            selected_year (int):
            The year for which the power duration curve is plotted.
        Returns:
            None: Displays a plot of the power duration curve.
        """

        power_output = []  # List to store power output values
        time_plot = []  # List to store time values for the plot
        for i in range(len(time)):  # Loop through each time point
            if years[i] == selected_year:  # Check if the year matches
                power = self.get_power(windspeed[i])  # Calculate the power
                power_output.append(power)  # Append the power output
                time_plot.append(time[i])  # Append the time to the list
        time_plot = time_plot - time_plot[0]  # Normalize time to start from 0

        # Sort the power output in descending order
        sorted_power = np.sort(power_output)[::-1]
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(time_plot, sorted_power, label='Power Output Sorted')
        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('Power Output (kW)')
        ax.set_title('Duration Curve for Power Output of Wind Turbine in year '
                     + str(selected_year))
        ax.legend()
        ax.grid(True)
        plt.show(block=True)
        return


# Functions that is not part of any classes
def combine_wind_data(wind_data_paths, location, component_names):
    """
    Combines wind speed and direction data from multiple files
      for a specific location.
    Args:
        wind_data_paths (list of str): Paths to the wind data files.
        location (tuple): Coordinates (latitude, longitude) of the location.
        component_names (list of str): Names of the wind components to process.
    Returns:
        tuple:
            - speed (dict): Wind speed at 10m and 100m heights as numpy arrays.
            - direction (dict): Wind direction at 10m and 100m as numpy arrays.
    """
    speed = {}
    direction = {}
    speed_10m = []
    speed_100m = []
    direction_10m = []
    direction_100m = []
    for wind_path in wind_data_paths:  # Loop through each wind data file
        winddata = WindData(wind_path)
        # Compute wind speed and direction for the location
        speed_loc, direction_loc = winddata.compute_wind_speed_direction(
            location, component_names)
        # Append the wind speed and direction data to the lists
        speed_10m.extend(speed_loc['10m'])
        speed_100m.extend(speed_loc['100m'])
        direction_10m.extend(direction_loc['10m'])
        direction_100m.extend(direction_loc['100m'])
    speed['10m'] = np.array(speed_10m)
    speed['100m'] = np.array(speed_100m)
    direction['10m'] = np.array(direction_10m)
    direction['100m'] = np.array(direction_100m)
    return speed, direction


def combine_time(wind_data_paths):
    """
    Combines time data from multiple wind data files.
    Args:
        wind_data_paths (list of str): List of file paths to wind data files.
    Returns:
        tuple: A tuple containing:
            - hours_combined (np.ndarray): Combined array of hourly time data.
            - years_combined (np.ndarray): Combined array of yearly time data.
    """
    hours_combined = []
    years_combined = []
    for wind_path in wind_data_paths:  # Loop through each wind data file
        winddata = WindData(wind_path)
        hours, years = winddata.get_time()  # Get the time data from the file
        # If this is not the first file, adjust the hours to avoid overlap
        if len(hours_combined) > 1:
            hours_combined.extend(hours + hours_combined[-1] + 1)
        else:
            hours_combined.extend(hours)
        years_combined.extend(years)  # Append the years to the list
    hours_combined = np.array(hours_combined)
    years_combined = np.array(years_combined)
    return hours_combined, years_combined


def windspeed_at_height(wind_speed, height, alpha):
    """
    Calculate wind speed at a given height using the power law.
    """
    if height > 55:  # If height is greater than 55m, use 100m wind speed
        u = wind_speed['100m']*(height/100)**alpha
    else:  # Otherwise, use 10m wind speed
        u = wind_speed['10m']*(height/10)**alpha
    return u


def fit_and_plot_weibull(wind_speeds, years, selected_year, plot):
    """
    Fits a Weibull distribution to wind speed data and plots the PDF.

    Parameters:
    wind_speeds (array-like): Array of wind speed values [m/s].

    Returns:
    shape (float): Weibull shape parameter (k).
    scale (float): Weibull scale parameter (A).
    """
    windspeed_selected_year = []  # List to wind speeds for the selected year
    for i in range(len(wind_speeds)):  # Loop through each wind speed value
        if years[i] == selected_year:  # Check if the year matches
            windspeed_selected_year.append(wind_speeds[i])
    windspeed_selected_year = np.array(windspeed_selected_year)

    # Fit the Weibull distribution (fix location to 0)
    shape, loc, scale = weibull_min.fit(windspeed_selected_year, floc=0)

    # Generate x values for plotting
    x = np.linspace(0, max(windspeed_selected_year) + 2, 100)
    pdf = weibull_min.pdf(x, shape, loc, scale)
    if plot == True:  # If plot is True, create a histogram and PDF plot
        # Plot
        plt.figure(figsize=(8, 5))
        plt.hist(windspeed_selected_year, bins=100, density=True, alpha=0.6,
                 color='skyblue', label='Histogram')
        plt.plot(x, pdf, 'r-', lw=2,
                 label=f'Weibull PDF\nk={shape:.2f}, A={scale:.2f}')
        plt.xlabel('Wind Speed [m/s]')
        plt.ylabel('Probability Density')
        plt.title('Weibull Distribution Fit to Wind Speed Data for year '
                  + str(selected_year))
        plt.legend()
        plt.grid(True)
        plt.show(block=False)

    return shape, scale


def plot_wind_rose(wind_speed, wind_direction):
    """
    Plots a wind rose diagram based on wind speed and direction data.
    Parameters:
    wind_speed (array-like): Array of wind speed values.
    wind_direction (array-like): Array of wind direction values in degrees.
    Returns:
    None: The function creates a plot but does not return any value.
    """
    # Create a WindroseAxes object for the wind rose plot
    ax = WindroseAxes.from_ax()
    ax.bar(wind_direction, wind_speed, normed=True,
           opening=0.8, edgecolor="white")
    ax.set_title(
        "Wind Rose Diagram showing wind speed and direction for all years")
    ax.set_legend()
    plt.show(block=False)
    return
