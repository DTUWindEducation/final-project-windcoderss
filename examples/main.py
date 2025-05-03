import final_project
import numpy as np
import os
import matplotlib.pyplot as plt

# Change the NetCDF4 file here
winddata_file_path = "./inputs/2006-2008.nc"
# Change the windturbine data file here
windturbine_file_path = "./inputs/NREL_Reference_5MW_126.csv"

# Coordinates for the location we want to get the speed and direction
loc_lon = 7.9061
loc_lat = 55.5297

# Power law exponent value with a standard value of 0.143
alpha = 0.143
# The height we want to estimate the wind speed
# which could be the hub height of the 5MW windturbine
height_z = 90

# Is the availability of the windturbine (1 for this project)
eta = 1
# The components name in the NetCDF4 file
component_name = ['u10', 'v10', 'u100', 'v100']

# Year that is selected to plot and calculate AEP,
# must be between 1997 and 2008
selected_year = 1998

all_years = [1997, 1998, 1999, 2000, 2001, 2002, 
             2003, 2004, 2005, 2006, 2007, 2008]
plot = True

winddata = final_project.WindData(winddata_file_path)  # Load NetCDF4 file

# Get the latitude and longitudes from the winddata
latitudes = winddata.get_latitude()
longitudes = winddata.get_longitudes()

# The four provided locations:
location1 = np.array([latitudes[1], longitudes[1]])
location2 = np.array([latitudes[0], longitudes[1]])
location3 = np.array([latitudes[0], longitudes[0]])
location4 = np.array([latitudes[1], longitudes[0]])

nc_files = ['./inputs/' + f for f in os.listdir('./inputs/') if
            f.endswith('.nc')]

# Get the speed and direction of the wind components for all locations
print("Combining wind data for all locations")
speed_loc1, direction_loc1 = final_project.combine_wind_data(
    nc_files, location1, component_name)
speed_loc2, direction_loc2 = final_project.combine_wind_data(
    nc_files, location2, component_name)
speed_loc3, direction_loc3 = final_project.combine_wind_data(
    nc_files, location3, component_name)
speed_loc4, direction_loc4 = final_project.combine_wind_data(
    nc_files, location4, component_name)

# Get the time vector from the NetCDF4 files
print("Combining time data")
time, years = final_project.combine_time(nc_files)

# Using power law to calculate wind speed at height z for all locations
print("Calculating wind speed at height z")
wind_speed_height_z_loc1 = final_project.windspeed_at_height(
    speed_loc1, height_z, alpha)
wind_speed_height_z_loc2 = final_project.windspeed_at_height(
    speed_loc2, height_z, alpha)
wind_speed_height_z_loc3 = final_project.windspeed_at_height(
    speed_loc3, height_z, alpha)
wind_speed_height_z_loc4 = final_project.windspeed_at_height(
    speed_loc4, height_z, alpha)

# Gather the speed and direction for height z
speed_locs = [wind_speed_height_z_loc1, wind_speed_height_z_loc2,
              wind_speed_height_z_loc3, wind_speed_height_z_loc4]
direction_locs = [direction_loc1['100m'], direction_loc2['100m'],
                  direction_loc3['100m'], direction_loc4['100m']]

# Interpolating wind speed and direction
# at the specified location (Hornsrev) using the wind data
print("Interpolating wind speed and direction at the specified location")
Hornsrev_speed, Hornsrev_direction = winddata.interpolate_at_loc(
    speed_locs, direction_locs, loc_lon, loc_lat)

# Fit and plot the weibull distribution for wind speed
# at the selected location and height
print("Fitting and plot Weibull distribution")
Hornsrev_weibull_shape, Hornsrev_weibull_scale = final_project.fit_and_plot_weibull(Hornsrev_speed, years, selected_year, plot)

# Plot a wind rose diagram for the selected location and height for all years
print("Plotting wind rose diagram")
final_project.plot_wind_rose(Hornsrev_speed, Hornsrev_direction)

# Define the windturbine object
windturbine = final_project.WindTurbine(windturbine_file_path)

# Compute AEP for the selected windturbine at the selected location and height
print("Calculating AEP")
AEP = windturbine.compute_AEP(Hornsrev_weibull_shape, Hornsrev_weibull_scale, eta)
print(f"The Annual Energy Production (AEP) for the year {selected_year} is {AEP:.2f} kWh.")

# Plot power output of the wind turbine
print("Compareing AEP for all years")
windturbine.compare_AEP(Hornsrev_speed, eta, years, all_years)

# plot duration curve of power output
print("Plotting power duration curve for the selected year")
windturbine.plot_power_duration_curve(
    Hornsrev_speed, time, years, selected_year)

print("All done!")
