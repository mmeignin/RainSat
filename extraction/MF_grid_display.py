#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 12:53:11 2023

Code that co-localizes SEVIRI data with MF radar mosaic data at 1 km.
Requires: satpy, cartopy, netCDF4

@author: meignin
"""
import os
from satpy.scene import Scene
import datetime
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from netCDF4 import Dataset
import random

# Load radar mosaic grid data (MFGrid)
mf_grid_file = '/net/nfs/precipitations/nicoEUMETSAT/grid_france_1536-1km_200906110000-219912310000.nc'

# Open and read the grid file
nc_fid = Dataset(mf_grid_file, 'r')
radar_latitudes = nc_fid.variables['latitude'][:,:]
radar_longitudes = nc_fid.variables['longitude'][:,:]
nc_fid.close()
for i in range(1,80):
    # Define years, months, and days
    years =  ["2015", "2011", "2013", "2021", "2012", "2015", "2016", "2017", "2018", "2014", "2019", "2020", "2009"]
    months = list(range(1, 13))  # All months from 1 to 12
    days = list(range(1, 32))  # All days from 1 to 31

    # Randomly select a year, month, and day
    year = random.choice(years)
    year = '2019'
    month = random.choice(months)
    day = random.choice(days)

    # Format the year, month, and day for file naming
    year_str = str(year)
    month_str = f"{month:02d}"
    day_str = f"{day:02d}"

    print(f"Processing for {year_str}-{month_str}-{day_str}")

    # Reading the MF Mosaic data
    file_pattern = f'/net/nfs/precipitations/SEDOO/Mosaiques/LMDO_1536_1km_5min/{year_str}/{year_str}{month_str}/{year_str}{month_str}{day_str}/cumul_france_1536-1km-5min_{year_str}{month_str}{day_str}*.nc'
    file_list = glob(file_pattern)
    if len(file_list) <1:
        break
    random_file_per_time = random.choice(file_list)

    
    print(f"Processing file { 1}/{len(file_list)}")

    # Open the file and read the necessary variables
    nc_fid = Dataset(random_file_per_time, 'r')
    latitudes = nc_fid.variables['radar_image_latitude'][0, :]
    longitudes = nc_fid.variables['radar_image_longitude'][0, :]
    rainfall_rate = nc_fid.variables['ACRR'][:,:,:]  # Accumulated rainfall rate (in 1/100 mm)
    quality_flag = nc_fid.variables['QUALITY'][:,:,:]  # Quality flag
    rainfall_rate = np.transpose(rainfall_rate[0,:,:].astype(float))
    quality_flag = np.transpose(quality_flag[0,:,:].astype(int))

    # Preprocessing the rainfall data
    rainfall_rate[rainfall_rate > 65000] = -99.0
    rainfall_rate[quality_flag > 200] = -98.0
    rainfall_rate_full = np.copy(np.transpose(rainfall_rate))
    rainfall_rate[quality_flag < 0] = -97.0  # Apply quality threshold

    # Convert rainfall from 1/100 mm (over 5 minutes) to mm/hr
    rainfall_rate[rainfall_rate > 0.0] = rainfall_rate[rainfall_rate > 0.0] / 100.0 * 12.0
    rainfall_rate_full[rainfall_rate_full > 0] = rainfall_rate_full[rainfall_rate_full > 0] / 100.0 * 12.0

    # Select radar images with significant rainfall
    if np.sum((rainfall_rate >= 10.0) & (quality_flag >= 80)) < 10:
        print("Not enough significant rainfall detected. Skipping...")
        
    else:
        # Extract time information from the file
        time_var = nc_fid.variables['time']
        time_units = getattr(time_var, 'units')
        time_values = [int(time_units[14:18]), int(time_units[19:21]), int(time_units[22:24]),
                        int(time_units[25:27]), int(time_units[28:30]), int(time_units[31:33])]
        time_var = time_var[0]
        reference_time = datetime.datetime(time_values[0], time_values[1], time_values[2],
                                            time_values[3], time_values[4], time_values[5])

        # Calculate the actual timestamp for the mosaic data
        mosaic_time = reference_time + datetime.timedelta(seconds=int(time_var))
        print(f"Mosaic time: {mosaic_time}")

        # Display the rainfall rate (LRR) as an image
        plt.figure(figsize=(10, 8))
        plt.imshow(quality_flag, cmap='jet',
                    norm=LogNorm(vmin=0.01, vmax=50))
        plt.colorbar(label='Rainfall Rate (mm/hr)')
        plt.title(f"Rainfall Rate for {year_str}-{month_str}-{day_str} {mosaic_time.strftime('%H:%M:%S')}")
        #plt.show()
        plt.savefig(f"2019/rainfall_rate_{year_str}-{month_str}-{day_str}-{mosaic_time.strftime('%H%M%S')}.png")
        plt.close()
        # Close the netCDF file
        nc_fid.close()
        exit(0)