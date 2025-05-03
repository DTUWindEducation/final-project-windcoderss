[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/zjSXGKeR)
# Our Great Package

Team: Windcoderss

## Overview
The structure of our final project are as follows:<br />
FINAL-PROJECT-WINDCODERSS<br />
├── inputs/<br />
│   ├── 1997-1999.nc<br />
│   ├── 2000-2002.nc<br />
│   ├── 2003-2005.nc<br />
│   ├── 2006-2008.nc<br />
│   ├── NREL_Reference_5MW_126.csv<br />
│   └── NREL_Reference_15MW_126.csv <br />
├── outputs/<br />
├── src/final_project<br />
│   └── __init__.py<br />
├── tests/<br />
│   └── Tests.py<br />
├── examples/<br />
│   └── main.py<br />
├── .gitignore<br />
├── LICENSE<br />
├── README.md<br />
└── pyproject.toml<br />


## Quick-start guide
Before installing the following packages are required:<br />
* NetCDF4
* numpy
* matplotlib
* scipy
* windrose

Installation guide:<br />
In order to install the package do the following:<br />
1. Go to anaconda prompt.
2. Change directory to the package folder FINAL-PROJECT-WINDCODERSS.
3. Run: "pip install -e ."

This will install the package final_project

## Architecture

This package is designed to analyze and visualize wind turbine performance metrics, including power output, wind characteristiscs, and energy production over time using various statistical and plotting tools.<br />
The package consist of two main classes and a few standalone helper functions.
1. `WindData` Class<br />
    This class handles interactions with the NetCDF wind data files, including:
    * Opening datasets
    * Extracting latitude, longitude and time
    * Retrieveing wind component (`u`, `v`)
    * Calculating wind speed and direction
    * Bilinear interpolation for custom locations

2. `WindTurbine` Class <br />
    Handles wind turbine power data, including:
    * Reading turbine performance curves from .csv file
    * Interpolating power output for any wind speed
    * Computing AEP (Annual Energy Production)
    * Plotting power duration curves (Extra functions)
    * Compairing AEP over years (Extra functions)

## Peer review

[ADD TEXT HERE!]
