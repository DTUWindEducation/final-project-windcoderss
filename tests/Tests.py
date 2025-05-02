"""
Check that the functions and classes are working as intended
"""
import final_project
import numpy as np
DATA_DIR = "./inputs/1997-1999.nc"
WINDTURBINE_DIR = "./inputs/NREL_Reference_5MW_126.csv"

def test_WindData_class():
    """Check if class points at the desired input file"""
    # given
    path_resp_file = DATA_DIR
    # when
    winddata1 = final_project.WindData(path_resp_file)
    # then
    assert winddata1.data_file_path == path_resp_file  # check if winddata1 is the correct file

def test_get_rootgrop():
    """Check if the root grp is reached"""
    # given
    path_resp_file = DATA_DIR
    winddata1 = final_project.WindData(path_resp_file)
    # when
    rootgrp = winddata1.get_rootgrp()
    # then
    assert hasattr(rootgrp, 'variables')  # check if rootgrp.variables exists

def test_get_latitude():
    """ Check that the output is a np.array"""
     # given
    path_resp_file = DATA_DIR
    winddata1 = final_project.WindData(path_resp_file)
    # when
    rootgrp = winddata1.get_rootgrp()
    # then
    assert isinstance(rootgrp.variables['latitude'][:], np.ndarray)  # check if latitude is a numpy array
    assert np.any(np.isclose(rootgrp.variables['latitude'][:], 55.75))  # check if any value is close to 55.75
    assert np.any(np.isclose(rootgrp.variables['latitude'][:], 55.5))   # check if any value is close to 55.5

def test_get_longitude():
    """ Check that the output is a np.array"""
     # given
    path_resp_file = DATA_DIR
    winddata1 = final_project.WindData(path_resp_file)
    # when
    rootgrp = winddata1.get_rootgrp()
    # then
    assert len(rootgrp.variables['longitude'][:]) == 2  # check if longitude array has two values
    assert np.any(np.isclose(rootgrp.variables['longitude'][:], 7.75))  # check if any value is close to 7.75
    assert np.any(np.isclose(rootgrp.variables['longitude'][:], 8))   # check if any value is close to 8.0

def test_get_time():
    """Check that the function get time works"""
    # given
    path_resp_file = DATA_DIR
    winddata1 = final_project.WindData(path_resp_file)
    # when
    time = winddata1.get_time()
    assert isinstance(time, np.ndarray)

def test_get_components_of_wind():
    """Check if wind components is collected in the function"""
    # given
    component_name = ['u10', 'v10', 'u100', 'v100']
    path_resp_file = DATA_DIR
    winddata1 = final_project.WindData(path_resp_file)
    # when
    wind_data1 = winddata1.get_components_of_wind(component_name)
    # then
    assert isinstance(wind_data1, dict)  # check if wind_data1 is a dictionary
    assert 'locations' in wind_data1  # check if 'locations' key exists in the dictionary
    for component in component_name:
        assert component in wind_data1  # check if each component name exists in the dictionary
        assert len(wind_data1['u10']) == 4  # check if the component has size 4 which is the four locations

def test_compute_wind_speed_direction():
    """Check if wind speed and direction is calculated correctly"""
    # given
    component_name = ['u10', 'v10', 'u100', 'v100']
    path_resp_file = DATA_DIR
    winddata1 = final_project.WindData(path_resp_file)
    latitudes = winddata1.get_latitude()
    longitudes = winddata1.get_longitudes()
    location1 = np.array([latitudes[0], longitudes[0]])
    # when
    speed_loc1, direction_loc1 = winddata1.compute_wind_speed_direction(location1, component_name)
    # then
    assert isinstance(speed_loc1, dict)  # check if speed_loc1 is a dictionary
    assert isinstance(direction_loc1, dict)  # check if direction_loc1 is a dictionary
    assert len(speed_loc1) == 2  # check if speed_loc1 has two keys (10m and 100m)
    assert len(direction_loc1) == 2  # check if direction_loc1 has two keys (10m and 100m)

def test_interpolate_at_loc():
    """
    Test the interpolation of wind speed and direction at a specific location and height.
    """

    # given
    height = '10m'
    component_name = ['u10', 'v10', 'u100', 'v100']
 
    path_resp_file = DATA_DIR
    winddata1 = final_project.WindData(path_resp_file)

    latitudes = winddata1.get_latitude()
    longitudes = winddata1.get_longitudes()

    location2 = np.array([latitudes[0], longitudes[0]])
    location1 = np.array([latitudes[1], longitudes[0]])
    location4 = np.array([latitudes[0], longitudes[1]])
    location3 = np.array([latitudes[1], longitudes[1]])

    speed_loc1, direction_loc1 = winddata1.compute_wind_speed_direction(location1, component_name)
    speed_loc2, direction_loc2 = winddata1.compute_wind_speed_direction(location2, component_name)
    speed_loc3, direction_loc3 = winddata1.compute_wind_speed_direction(location3, component_name)
    speed_loc4, direction_loc4 = winddata1.compute_wind_speed_direction(location4, component_name)

    speed_locs = [speed_loc1[height], speed_loc2[height], speed_loc3[height], speed_loc4[height]]
    direction_locs = [direction_loc1[height], direction_loc2[height], direction_loc3[height], direction_loc4[height]]

    loc_lat = latitudes[0]
    loc_lon = longitudes[0]
    # when
    New_speed, New_direction = winddata1.interpolate_at_loc(speed_locs, direction_locs, loc_lat, loc_lon)
    # then
    assert isinstance(New_speed, np.ndarray)  # check if New_speed is a numpy array
    assert isinstance(New_direction, np.ndarray)  # check if New_direction is a numpy array
    assert len(New_speed) == len(speed_loc1['10m'])  # check if the length of New_speed matches the length of speed_locs
    assert len(New_direction) == len(direction_loc1['10m'])  # check if the length of New_direction matches the length of direction_locs

def test_windspeed_at_height():
    """
    Test the calculation of wind speed at a specific height using given wind data and parameters.
    """
    # Given
    path_resp_file = DATA_DIR
    component_name = ['u10', 'v10', 'u100', 'v100']
    winddata1 = final_project.WindData(path_resp_file)
    alpha = 0.143
    height_z = 15
    latitudes = winddata1.get_latitude()
    longitudes = winddata1.get_longitudes()
    location1 = np.array([latitudes[1], longitudes[1]])
    speed_loc1, direction_loc1 = winddata1.compute_wind_speed_direction(location1, component_name)

    # when
    wind_speed_height_z_loc1 = final_project.windspeed_at_height(speed_loc1, height_z, alpha)
    # Then
    assert len(wind_speed_height_z_loc1) == len(speed_loc1['10m'])  # check if lengths match
    assert np.all(wind_speed_height_z_loc1 > speed_loc1['10m'])  # check if all values are larger

def test_fit_and_plot_weibull():
    """
    Test the function that fits and plots the weibull distribution
    """
    # Given
    path_resp_file = DATA_DIR
    component_name = ['u10', 'v10', 'u100', 'v100']
    winddata1 = final_project.WindData(path_resp_file)
    latitudes = winddata1.get_latitude()
    longitudes = winddata1.get_longitudes()
    location1 = np.array([latitudes[1], longitudes[1]])
    speed_loc1, _ = winddata1.compute_wind_speed_direction(location1, component_name)
    wind_speeds = speed_loc1['10m']
    # when 
    shape, scale = final_project.fit_and_plot_weibull(wind_speeds)
    print(shape)
    # Then
    assert isinstance(shape, float) # Shape parameter should be a float
    assert isinstance(scale, float) # Scale parameter should be a float
    assert shape > 0 # Shape parameter should be positive
    assert scale > 0 # Scale parameter should be positive

def test_plot_wind_rose():
    """
    Test the function that fits and plots the weibull distribution
    """
    # Given
    path_resp_file = DATA_DIR
    component_name = ['u10', 'v10', 'u100', 'v100']
    winddata1 = final_project.WindData(path_resp_file)
    latitudes = winddata1.get_latitude()
    longitudes = winddata1.get_longitudes()
    location1 = np.array([latitudes[1], longitudes[1]])
    speed_loc1, direction_loc1 = winddata1.compute_wind_speed_direction(location1, component_name)
    # When / Then
    try:
        final_project.plot_wind_rose(speed_loc1['10m'], direction_loc1['10m'])
    except Exception as e:
        assert False, f"plot_wind_rose raised an exception: {e}"

def test_plot_power_output():
    """
    Test the function that plots the power output of the wind turbine
    """
    # Given
    path_resp_file = DATA_DIR
    component_name = ['u10', 'v10', 'u100', 'v100']
    winddata1 = final_project.WindData(path_resp_file)
    path_windturbine_file = WINDTURBINE_DIR
    wind_turbine = final_project.WindTurbine(path_windturbine_file)
    latitudes = winddata1.get_latitude()
    longitudes = winddata1.get_longitudes()
    location1 = np.array([latitudes[1], longitudes[1]])
    speed_loc1, _ = winddata1.compute_wind_speed_direction(location1, component_name)
    time = winddata1.get_time()
    # When / Then
    try:
        wind_turbine.plot_power_output(speed_loc1['10m'], time)
    except Exception as e:
        assert False, f"plot_power_output raised an exception: {e}"
    
    
def test_plot_power_duration_curve():
    """
    Test the function that plots the power output of the wind turbine
    """
    # Given
    path_resp_file = DATA_DIR
    component_name = ['u10', 'v10', 'u100', 'v100']
    winddata1 = final_project.WindData(path_resp_file)
    path_windturbine_file = WINDTURBINE_DIR
    wind_turbine = final_project.WindTurbine(path_windturbine_file)
    latitudes = winddata1.get_latitude()
    longitudes = winddata1.get_longitudes()
    location1 = np.array([latitudes[1], longitudes[1]])
    speed_loc1, _ = winddata1.compute_wind_speed_direction(location1, component_name)
    time = winddata1.get_time()
    # When / Then
    try:
        wind_turbine.plot_power_duration_curve(speed_loc1['10m'], time)
    except Exception as e:
        assert False, f"plot_power_duration_curve raised an exception: {e}"