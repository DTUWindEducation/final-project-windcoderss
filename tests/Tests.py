"""
Check that the functions and classes are working as intended
"""
import final_project
import numpy as np
DATA_DIR = "./inputs/1997-1999.nc"

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
    assert np.any(np.isclose(rootgrp.variables['latitude'][:], 7.75))  # check if any value is close to 7.75
    assert np.any(np.isclose(rootgrp.variables['latitude'][:], 8.0))   # check if any value is close to 8.0

def test_get_longitude():
    """ Check that the output is a np.array"""
     # given
    path_resp_file = DATA_DIR
    winddata1 = final_project.WindData(path_resp_file)
    # when
    rootgrp = winddata1.get_rootgrp()
    # then
    assert len(rootgrp.variables['longitude'][:]) == 2  # check if longitude array has two values
    assert np.any(np.isclose(rootgrp.variables['longitude'][:], 55.75))  # check if any value is close to 7.75
    assert np.any(np.isclose(rootgrp.variables['longitude'][:], 55.5))   # check if any value is close to 8.0

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