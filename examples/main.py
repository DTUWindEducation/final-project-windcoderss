import final_project

winddata_file_path = "./inputs/1997-1999.nc"
winddata1 = final_project.WindData(winddata_file_path)

latitudes1 = winddata1.get_latitude()

print(latitudes1)