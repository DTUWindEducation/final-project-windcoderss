import final_project
import numpy as np
import matplotlib.pyplot as plt

winddata_file_path = "./inputs/2000-2002.nc"
winddata1 = final_project.WindData(winddata_file_path)

latitudes = winddata1.get_latitude()
longitudes = winddata1.get_longitudes()

location1 = np.array([latitudes[0], longitudes[0]])
location2 = np.array([latitudes[1], longitudes[0]])
location3 = np.array([latitudes[0], longitudes[1]])
location4 = np.array([latitudes[1], longitudes[1]])

component_name = ['u10', 'v10', 'u100', 'v100']
wind_data1 = winddata1.get_components_of_wind(component_name)

speed, direction = winddata1.compute_wind_speed_direction(location1, component_name)

# Extract wind speeds at 10m and 100m
speed_10m = speed['10m']
speed_100m = speed['100m']

# Extract wind directions at 10m and 100m
direction_10m = direction['10m']
direction_100m = direction['100m']

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

