import final_project

winddata_file_path = "./inputs/1997-1999.nc"
winddata1 = final_project.WindData(winddata_file_path)

latitudes1 = winddata1.get_latitude()
longitudes1 = winddata1.get_longitudes()

component_name = ['u10', 'v10', 'u100', 'v100']
wind_data1 = winddata1.get_components_of_wind(component_name)
print("Wind Data:", wind_data1['u10'])

