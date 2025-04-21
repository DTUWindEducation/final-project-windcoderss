import final_project
import numpy as np
import matplotlib.pyplot as plt

winddata_file_path = "./inputs/1997-1999.nc"

# Coordinates for the location we want to get the speed and direction
loc_lon = 7.9061
loc_lat = 55.5297
height = ['10m']         # 10m or 100m

winddata1 = final_project.WindData(winddata_file_path)

latitudes = winddata1.get_latitude()
longitudes = winddata1.get_longitudes()

location1 = np.array([latitudes[0], longitudes[0]])
location2 = np.array([latitudes[1], longitudes[0]])
location3 = np.array([latitudes[0], longitudes[1]])
location4 = np.array([latitudes[1], longitudes[1]])

component_name = ['u10', 'v10', 'u100', 'v100']
wind_data1 = winddata1.get_components_of_wind(component_name)

speed_loc1, direction_loc1 = winddata1.compute_wind_speed_direction(location1, component_name)
speed_loc2, direction_loc2 = winddata1.compute_wind_speed_direction(location2, component_name)
speed_loc3, direction_loc3 = winddata1.compute_wind_speed_direction(location3, component_name)
speed_loc4, direction_loc4 = winddata1.compute_wind_speed_direction(location4, component_name)

speed_locs = [speed_loc1, speed_loc2, speed_loc3, speed_loc4]
direction_locs = [direction_loc1, direction_loc2, direction_loc3, direction_loc4]


# Extract wind speeds at 10m and 100m
speed_10m = speed_loc1['10m']
speed_100m = speed_loc1['100m']

# Extract wind directions at 10m and 100m
direction_10m = direction_loc1['10m']
direction_100m = direction_loc1['100m']

# Create a figure with two subplots
fig, axs = plt.subplots(2, 1, figsize=(10, 12))

# Plot wind speed on the first subplot
axs[0].plot(speed_10m, label='10m Wind Speed')
axs[0].plot(speed_100m, label='100m Wind Speed')
axs[0].set_xlabel('Time Index')
axs[0].set_ylabel('Wind Speed (m/s)')
axs[0].set_title('Wind Speed at 10m and 100m')
axs[0].legend()

# Plot wind direction on the second subplot
axs[1].plot(direction_10m, label='10m Wind Direction')
axs[1].plot(direction_100m, label='100m Wind Direction')
axs[1].set_xlabel('Time Index')
axs[1].set_ylabel('Wind Direction (degrees)')
axs[1].set_title('Wind Direction at 10m and 100m')
axs[1].legend()

# Adjust layout and show the plot
plt.tight_layout()
plt.show()

Horsrev_10m = wind_data1.interpolate_at_loc()
print(interpolated_speed, interpolated_direction)