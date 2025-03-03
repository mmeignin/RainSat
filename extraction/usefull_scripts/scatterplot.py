import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import random
from glob import glob
import sys
import itertools
import glob


# Function to generate a scatter plot of Rainfall vs Quality Flag for a single file
def plot_rainfall_vs_quality_single(rainfall_rates, quality_flags, filename):
    plt.figure(figsize=(10, 6))
    plt.scatter(rainfall_rates, quality_flags, c='blue', label='Rainfall Rate vs Quality Flag', alpha=0.6, s=1)
    plt.xlabel('Rainfall Rate (mm/hr)')
    plt.ylabel('Quality Flag')
    plt.ylim(0, 100)
    plt.title('Rainfall Rate vs Quality Flag (Single File)')
    plt.grid(True)
    plt.legend()
    plt.savefig(filename)
    plt.close()



# Define years, months, and days to analyze
years = ["2015", "2011", "2013", "2021", "2012", "2015", "2016", "2017", "2018", "2014", "2019", "2020", "2009"]
months = range(1, 13)  # Months from 1 to 12
days = range(1, 32)  # All days from 1 to 31

# List to hold all matching files
all_matching_files = []

# Loop through each combination of year, month, and day
for year, month, day in itertools.product(years, months, days):
    # Format month and day as 2-digit strings
    month_str = f"{month:02d}"
    day_str = f"{day:02d}"

    # Construct the file pattern
    file_pattern = f'/net/nfs/precipitations/SEDOO/Mosaiques/LMDO_1536_1km_5min/{year}/{year}{month_str}/{year}{month_str}{day_str}/cumul_france_1536-1km-5min_{year}{month_str}{day_str}*.nc'

    # Use glob to find all files matching the pattern
    matching_files = glob.glob(file_pattern)

    # If matching files are found, add them to the list
    if matching_files:
        all_matching_files.extend(matching_files)

# Print the results
if all_matching_files:
    print(f"Found {len(all_matching_files)} matching files:")
else:
    print("No matching files found.")

# Initialize lists to store grid point rainfall rates and quality flags for one file
rainfall_rates = []
quality_flags = []

# Randomly pick one file from the list
num_files = min(100, len(all_matching_files))
random_files = random.sample(all_matching_files, num_files)

print(f"Processing file {random_files}")



for random_file in random_files:
    print(f"Processing file {random_file}",file=sys.stdout)
    # Open the file and read variables
    nc_fid = Dataset(random_file, 'r')
    rainfall_rate = nc_fid.variables['ACRR'][:,:,:]  # Accumulated rainfall rate (in 1/100 mm)
    quality_flag = nc_fid.variables['QUALITY'][:,:,:]  # Quality flag
    rainfall_rate = np.transpose(rainfall_rate[0,:,:].astype(float))
    quality_flag = np.transpose(quality_flag[0,:,:].astype(int))

    # Preprocessing the rainfall data
    rainfall_rate[rainfall_rate > 65000] = -99.0
    rainfall_rate[quality_flag > 200] = -98.0
    rainfall_rate[quality_flag < 0] = -97.0  # Apply quality threshold
    rainfall_rate[np.isnan(rainfall_rate)] = -96.0

    # Convert rainfall from 1/100 mm (over 5 minutes) to mm/hr
    rainfall_rate[rainfall_rate > 0.0] = rainfall_rate[rainfall_rate > 0.0] / 100.0 * 12.0

    # Check if there are enough significant rainfall points
    if np.sum((rainfall_rate >= 10.0) & (quality_flag >= 80)) < 10:
        print(f"Not enough significant rainfall detected for {year}-{month_str}-{day_str}. Skipping file...")
        continue    
    # Append all grid points' values to the lists (no averaging)
    rainfall_rates.append(rainfall_rate.flatten())  # Flatten the 2D array to 1D for scatter
    quality_flags.append(quality_flag.flatten())    # Flatten the 2D array to 1D for scatter

    # Close the netCDF file
    nc_fid.close()

# Convert lists to numpy arrays for analysis
rainfall_rates = np.concatenate(rainfall_rates)
quality_flags = np.concatenate(quality_flags)

# Scatter Plot: Rainfall Rate vs Quality Flag (single file)
scatter_filename = 'results/rain_quality_results/rainfall_quality_scatter.png'
plot_rainfall_vs_quality_single(rainfall_rates, quality_flags, scatter_filename)

print(f"Plots saved as {scatter_filename}")
