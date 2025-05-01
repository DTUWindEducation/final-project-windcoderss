import final_project
import numpy as np
import matplotlib.pyplot as plt

winddata_file_path = "./inputs/2006-2008.nc"
windturbine_file_path = "./inputs/NREL_Reference_5MW_126.csv"

# Coordinates for the location we want to get the speed and direction
loc_lon = 7.9061
loc_lat = 55.5297
height = '10m'         # 10m or 100m
alpha = 0.143           # standard value
height_z = 56
eta = 1
component_name = ['u10', 'v10', 'u100', 'v100']

winddata1 = final_project.WindData(winddata_file_path)

latitudes = winddata1.get_latitude()
longitudes = winddata1.get_longitudes()

location3 = np.array([latitudes[0], longitudes[0]])
location4 = np.array([latitudes[1], longitudes[0]])
location2 = np.array([latitudes[0], longitudes[1]])
location1 = np.array([latitudes[1], longitudes[1]])

speed_loc1, direction_loc1 = winddata1.compute_wind_speed_direction(location1, component_name)
speed_loc2, direction_loc2 = winddata1.compute_wind_speed_direction(location2, component_name)
speed_loc3, direction_loc3 = winddata1.compute_wind_speed_direction(location3, component_name)
speed_loc4, direction_loc4 = winddata1.compute_wind_speed_direction(location4, component_name)

speed_locs = [speed_loc1, speed_loc2, speed_loc3, speed_loc4]
direction_locs = [direction_loc1, direction_loc2, direction_loc3, direction_loc4]

# Interpolating wind speed and direction at the specified location (Hornsrev) using the wind data
Hornsrev_10m_speed, Hornsrev_10m_direction = winddata1.interpolate_at_loc(speed_locs, direction_locs, height, loc_lon, loc_lat)

# Using power law to calculate wind speed at height z for all locations
wind_speed_height_z_loc1 = final_project.windspeed_at_height(speed_loc1, height_z, alpha)
wind_speed_height_z_loc2 = final_project.windspeed_at_height(speed_loc2, height_z, alpha)
wind_speed_height_z_loc3 = final_project.windspeed_at_height(speed_loc3, height_z, alpha)
wind_speed_height_z_loc4 = final_project.windspeed_at_height(speed_loc4, height_z, alpha)

Hornsrev_weibull_shape, Hornsrev_weibull_scale = final_project.fit_and_plot_weibull(Hornsrev_10m_speed)

final_project.plot_wind_rose(Hornsrev_10m_speed, Hornsrev_10m_direction)

windturbine1 = final_project.WindTurbine(windturbine_file_path)
AEP = windturbine1.get_AEP(Hornsrev_weibull_shape, Hornsrev_weibull_scale, eta)
print(AEP)

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
axs[0].plot(speed_10m_loc1, label='Loc1 10m Wind Speed')
axs[0].plot(speed_10m_loc2, label='Loc2 10m Wind Speed')
axs[0].plot(speed_10m_loc3, label='Loc3 10m Wind Speed')
axs[0].plot(speed_10m_loc4, label='Loc4 10m Wind Speed')
axs[0].plot(Hornsrev_10m_speed, label='Hornsrev 10m Wind Speed', linestyle='--')
axs[0].set_xlabel('Time Index')
axs[0].set_ylabel('Wind Speed (m/s)')
axs[0].set_title('Wind Speed at 10m for All Locations and Hornsrev')
axs[0].legend()

# Plot wind direction at 10m for all locations and Hornsrev interpolated direction
axs[1].plot(direction_10m_loc1, label='Loc1 10m Wind Direction')
axs[1].plot(direction_10m_loc2, label='Loc2 10m Wind Direction')
axs[1].plot(direction_10m_loc3, label='Loc3 10m Wind Direction')
axs[1].plot(direction_10m_loc4, label='Loc4 10m Wind Direction')
axs[1].plot(Hornsrev_10m_direction, label='Hornsrev 10m Wind Direction', linestyle='--')
axs[1].set_xlabel('Time Index')
axs[1].set_ylabel('Wind Direction (degrees)')
axs[1].set_title('Wind Direction at 10m for All Locations and Hornsrev')
axs[1].legend()

# Adjust layout and show the plot
plt.tight_layout()
plt.show()