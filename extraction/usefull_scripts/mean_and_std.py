import os
import random
import sys
import numpy as np
from glob import glob
from netCDF4 import Dataset
import itertools

# Constants
MF_GRID_FILE = '/net/nfs/precipitations/nicoEUMETSAT/grid_france_1536-1km_200906110000-219912310000.nc'
GRID_LAT_LONG_VARS = ['latitude', 'longitude']
RAINFALL_VAR = 'ACRR'
QUALITY_VAR = 'QUALITY'
QUALITY_THRESHOLD = 200
RAINFALL_THRESHOLD = 10.0
MM_CONVERSION = 100.0
HOUR_CONVERSION = 12.0

# Function to load the radar grid data
def load_radar_grid():
    with Dataset(MF_GRID_FILE, 'r') as nc_fid:
        latitudes = nc_fid.variables[GRID_LAT_LONG_VARS[0]][:]
        longitudes = nc_fid.variables[GRID_LAT_LONG_VARS[1]][:]
    return latitudes, longitudes

# Function to gather matching files based on year, month, day
def get_matching_files(years, months, days):
    all_matching_files = []
    for year, month, day in itertools.product(years, months, days):
        file_pattern = f'/net/nfs/precipitations/SEDOO/Mosaiques/LMDO_1536_1km_5min/{year}/{year}{month:02d}/{year}{month:02d}{day:02d}/cumul_france_1536-1km-5min_{year}{month:02d}{day:02d}*.nc'
        all_matching_files.extend(glob(file_pattern))
    return all_matching_files

# Function to process each selected file
def process_file(file):
    try:
        with Dataset(file, 'r') as nc_fid:
            print(nc_fid)
            exit(0)
            latitudes = nc_fid.variables['radar_image_latitude'][0, :]
            longitudes = nc_fid.variables['radar_image_longitude'][0, :]
            rainfall_rate = nc_fid.variables[RAINFALL_VAR][:,:,:].astype(float)
            quality_flag = nc_fid.variables[QUALITY_VAR][:,:,:].astype(int)

        # Preprocess rainfall data
        rainfall_rate = np.transpose(rainfall_rate[0, :, :])
        quality_flag = np.transpose(quality_flag[0, :, :])
        rainfall_rate[rainfall_rate > 65000] = -99.0
        rainfall_rate[quality_flag > QUALITY_THRESHOLD] = -98.0
        rainfall_rate[quality_flag < 0] = -97.0

        # Convert rainfall from 1/100 mm to mm/hr
        rainfall_rate[rainfall_rate > 0] *= HOUR_CONVERSION / MM_CONVERSION

        return rainfall_rate, quality_flag

    except Exception as e:
        print(f"Error processing {file}: {e}", file=sys.stderr, flush=True)
        return None, None

# Update streaming statistics
def update_statistics(rainfall_rate, mean_rainfall, M2, n_data_points):
    valid_rainfall = rainfall_rate[rainfall_rate > 0]
    for rainfall in valid_rainfall.flatten():  # Flatten to iterate over all elements
        if rainfall > 0:
            n_data_points += 1
            delta = rainfall - mean_rainfall
            mean_rainfall += delta / n_data_points
            M2 += delta * (rainfall - mean_rainfall)
    return mean_rainfall, M2, n_data_points

# Main function to process all files
def process_files():
    # Load grid data
    latitudes, longitudes = load_radar_grid()

    # Define years, months, and days
    years = ["2015", "2011", "2013", "2021", "2012", "2015", "2016", "2017", "2018", "2014", "2019", "2020", "2009"]
    months = list(range(1, 13))
    days = list(range(1, 32))

    # Get matching files
    all_matching_files = get_matching_files(years, months, days)
    if not all_matching_files:
        print("No matching files found.", flush=True)
        sys.exit(1)

    print(f"Found {len(all_matching_files)} matching files.", flush=True)

    # Randomly select files
    random_files = random.sample(all_matching_files, min(200, len(all_matching_files)))

    # Initialize streaming statistics
    mean_rainfall = 0
    M2 = 0  # Sum of squared differences
    n_data_points = 0
    min_val = np.inf
    max_val = -np.inf

    # Process each selected file
    for file in random_files:
        print(f"Processing file {file}", file=sys.stderr, flush=True)

        # Process the file
        rainfall_rate, quality_flag = process_file(file)
        if rainfall_rate is None or quality_flag is None:
            continue

        # Select radar images with significant rainfall
        if np.sum((rainfall_rate >= RAINFALL_THRESHOLD) & (quality_flag >= 80)) < 10:
            print("Not enough significant rainfall detected. Skipping...")
            continue

        # Update min and max values
        rainfall_rate = rainfall_rate[rainfall_rate > 0]
        min_val = min(np.min(rainfall_rate), min_val)
        max_val = max(np.max(rainfall_rate), max_val)
        print(f"Current min: {min_val}, max: {max_val}")

        # Update streaming statistics
        mean_rainfall, M2, n_data_points = update_statistics(rainfall_rate, mean_rainfall, M2, n_data_points)

    # Calculate final statistics
    if n_data_points > 0:
        variance = M2 / n_data_points
        streaming_std_rainfall = np.sqrt(variance)
    else:
        streaming_std_rainfall = 0.0

    # Print final statistics
    print(f"Final min: {min_val}, max: {max_val}")
    print(f"Streaming mean: {mean_rainfall}")
    print(f"Streaming std: {streaming_std_rainfall}")

# Run the main function
if __name__ == "__main__":
    process_files()
