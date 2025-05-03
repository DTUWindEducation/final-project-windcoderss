import final_project
import numpy as np
import os
import matplotlib.pyplot as plt

winddata_file_path = "./inputs/2006-2008.nc" # change the NetCDF4 file here
windturbine_file_path = "./inputs/NREL_Reference_5MW_126.csv" # Change the windturbine data file here

loc_lon = 7.9061 # Coordinates for the location we want to get the speed and direction
loc_lat = 55.5297 
alpha = 0.143 # Power law exponent value with a standard value of 0.143
height_z = 90 # The height we want to estimate the wind speed which is the hub height of the 5MW windturbine
eta = 1 # Is the availability of the windturbine (1 for this project)
component_name = ['u10', 'v10', 'u100', 'v100'] # The components name in the NetCDF4 file
selected_year = 1998 # Year that is selected to plot and calculate AEP, must be between 1997 and 2008
all_years = [1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008]
plot = True

winddata = final_project.WindData(winddata_file_path) # Load NetCDF4 file

# Get the latitude and longitudes from the winddata
latitudes = winddata.get_latitude()
longitudes = winddata.get_longitudes()

# The four provided locations:
location1 = np.array([latitudes[1], longitudes[1]])
location2 = np.array([latitudes[0], longitudes[1]])
location3 = np.array([latitudes[0], longitudes[0]])
location4 = np.array([latitudes[1], longitudes[0]])

nc_files = ['./inputs/' + f for f in os.listdir('./inputs/') if f.endswith('.nc')]

# Get the speed and direction of the wind components for all locations
speed_loc1, direction_loc1 = final_project.combine_wind_data(nc_files, location1, component_name)
speed_loc2, direction_loc2 = final_project.combine_wind_data(nc_files, location2, component_name)
speed_loc3, direction_loc3 = final_project.combine_wind_data(nc_files, location3, component_name)
speed_loc4, direction_loc4 = final_project.combine_wind_data(nc_files, location4, component_name)

# Get the time vector from the NetCDF4 files
time, years = final_project.combine_time(nc_files)


# Using power law to calculate wind speed at height z for all locations
wind_speed_height_z_loc1 = final_project.windspeed_at_height(speed_loc1, height_z, alpha)
wind_speed_height_z_loc2 = final_project.windspeed_at_height(speed_loc2, height_z, alpha)
wind_speed_height_z_loc3 = final_project.windspeed_at_height(speed_loc3, height_z, alpha)
wind_speed_height_z_loc4 = final_project.windspeed_at_height(speed_loc4, height_z, alpha)

# Gather the speed and direction for height z
speed_locs = [wind_speed_height_z_loc1, wind_speed_height_z_loc2, wind_speed_height_z_loc3, wind_speed_height_z_loc4]
direction_locs = [direction_loc1['100m'], direction_loc2['100m'], direction_loc3['100m'], direction_loc4['100m']]

# Interpolating wind speed and direction at the specified location (Hornsrev) using the wind data
Hornsrev_speed, Hornsrev_direction = winddata.interpolate_at_loc(speed_locs, direction_locs, loc_lon, loc_lat)

# Fit and plot the weibull distribution for wind speed at the selected location and height
Hornsrev_weibull_shape, Hornsrev_weibull_scale = final_project.fit_and_plot_weibull(Hornsrev_speed, years, selected_year, plot)

# Plot a wind rose diagram for the selected location and height
final_project.plot_wind_rose(Hornsrev_speed, Hornsrev_direction)

# Define the windturbine object
windturbine = final_project.WindTurbine(windturbine_file_path)

# Compute AEP for the selected windturbine at the selected location and height
AEP = windturbine.get_AEP(Hornsrev_weibull_shape, Hornsrev_weibull_scale, eta)

print(AEP)

# Plot power output of the wind turbine
windturbine.compare_AEP(Hornsrev_speed, eta, years, all_years)

# plot duration curve of power output
windturbine.plot_power_duration_curve(Hornsrev_speed, time, years, selected_year)

# Extract wind speeds at 10m for all locations
speed_10m_loc1 = speed_loc1['10m']
speed_10m_loc2 = speed_loc2['10m']
speed_10m_loc3 = speed_loc3['10m']
speed_10m_loc4 = speed_loc4['10m']

# Extract wind directions at 10m for all locations
direction_10m_loc1 = direction_loc1['10m']
direction_10m_loc2 = direction_loc2['10m']
direction_10m_loc3 = direction_loc3['10m']
direction_10m_loc4 = direction_loc4['10m']

# Create a figure with two subplots
fig, axs = plt.subplots(2, 1, figsize=(10, 12))
# Plot wind speed at 10m for all locations and Hornsrev interpolated speed
axs[0].plot(time, speed_10m_loc1, label='Loc1 10m Wind Speed')
axs[0].plot(time, speed_10m_loc2, label='Loc2 10m Wind Speed')
axs[0].plot(time, speed_10m_loc3, label='Loc3 10m Wind Speed')
axs[0].plot(time, speed_10m_loc4, label='Loc4 10m Wind Speed')
axs[0].plot(time, Hornsrev_speed, label='Hornsrev 10m Wind Speed', linestyle='--')
axs[0].set_xlabel('Time (hours)')
axs[0].set_ylabel('Wind Speed (m/s)')
axs[0].set_title('Wind Speed at 10m for All Locations and Hornsrev')
axs[0].legend()

# Plot wind direction at 10m for all locations and Hornsrev interpolated direction
axs[1].plot(time, direction_10m_loc1, label='Loc1 10m Wind Direction')
axs[1].plot(time, direction_10m_loc2, label='Loc2 10m Wind Direction')
axs[1].plot(time, direction_10m_loc3, label='Loc3 10m Wind Direction')
axs[1].plot(time, direction_10m_loc4, label='Loc4 10m Wind Direction')
axs[1].plot(time, Hornsrev_direction, label='Hornsrev 10m Wind Direction', linestyle='--')
axs[1].set_xlabel('Time (hours)')
axs[1].set_ylabel('Wind Direction (degrees)')
axs[1].set_title('Wind Direction at 10m for All Locations and Hornsrev')
axs[1].legend()

# Adjust layout and show the plot
plt.tight_layout()
plt.show(block=True)
