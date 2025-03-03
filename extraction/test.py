import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import random
from glob import glob
import sys
import itertools
from matplotlib.colors import LogNorm

# Function to generate image plots of Rainfall Rate and Quality Flag
def plot_rainfall_and_quality_maps(rainfall_rate, quality_flag, filename):
    """
    Generates 2D heatmaps of rainfall rates and quality flags.

    Arguments:
    rainfall_rate -- 2D numpy array of rainfall rates (mm/hr)
    quality_flag -- 2D numpy array of quality flags
    filename -- String, name of the output plot file
    """
    plt.close('all')
    fig, axes = plt.subplots(1, 3, figsize=(12, 6))

    # Rainfall Rate Heatmap
    im1 = axes[0].imshow(rainfall_rate, cmap="jet", origin="upper", norm=LogNorm (vmin=0.1, vmax=np.max(rainfall_rate)))
    axes[0].set_title("Rainfall Rate (mm/hr)")
    fig.colorbar(im1, ax=axes[0], orientation="vertical")

    # Quality Flag Heatmap
    im2 = axes[1].imshow(quality_flag, cmap="viridis", origin="upper", vmin=0, vmax=100)
    axes[1].set_title("Quality Flag")
    fig.colorbar(im2, ax=axes[1], orientation="vertical")
    
    # Quality Flag Heatmap
    im3 = axes[2].imshow(rainfall_rate[rainfall_rate > 200.], cmap="jet", origin="upper")
    axes[2].set_title("high rr")
    fig.colorbar(im3, ax=axes[2], orientation="vertical")


    # Save and close figure
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close(fig) 

# Define years, months, and days to analyze
years =  ["2015", "2011", "2013", "2021", "2012", "2015", "2016", "2017", "2018", "2014", "2019", "2020", "2009"]
months = range(1, 13)  # Months from 1 to 12
days = range(1, 32)  # All days from 1 to 31

# List to hold all matching files
all_matching_files = []

# Loop through each combination of year, month, and day
for year, month, day in itertools.product(years, months, days):
    month_str = f"{month:02d}"
    day_str = f"{day:02d}"

    # Construct the file pattern
    file_pattern = f'/net/nfs/precipitations/SEDOO/Mosaiques/LMDO_1536_1km_5min/{year}/{year}{month_str}/{year}{month_str}{day_str}/cumul_france_1536-1km-5min_{year}{month_str}{day_str}*.nc'

    # Use glob to find all files matching the pattern
    matching_files = glob(file_pattern)

    if matching_files:
        all_matching_files.extend(matching_files)

# Print the results
if all_matching_files:
    print(f"Found {len(all_matching_files)} matching files.", flush=True)
else:
    print("No matching files found.", flush=True)
    sys.exit(1)  # Exit script if no files are found

# Randomly pick up to 200 files from the list
num_files = min(len(all_matching_files), len(all_matching_files))
random_files = random.sample(all_matching_files, num_files)

print(f"Processing {num_files} random files...", flush=True)

# Process each selected file
for random_file in random_files:
    try:
        print(f"Processing file {random_file}", file=sys.stdout, flush=True)

        # Open the NetCDF file
        with Dataset(random_file, 'r') as nc_fid:

            # Read variables
            rainfall_rate = nc_fid.variables['ACRR'][:,:,:]  # Accumulated rainfall rate (in 1/100 mm)
            quality_flag = nc_fid.variables['QUALITY'][:,:,:]  # Quality flag

            # Extract and transpose data
            rainfall_rate = np.transpose(rainfall_rate[0,:,:].astype(float))
            quality_flag = np.transpose(quality_flag[0,:,:].astype(int))

            # Preprocessing rainfall data
            rainfall_rate[rainfall_rate > 65000] = np.nan  # Mask extreme invalid values
            rainfall_rate[quality_flag > 200] = np.nan  # Exclude bad quality flags
            rainfall_rate[quality_flag < 0] = np.nan  # Apply quality threshold
            rainfall_rate[np.isnan(rainfall_rate)] = np.nan  # Handle NaNs

            # Convert rainfall from 1/100 mm (over 5 minutes) to mm/hr
            rainfall_rate = rainfall_rate / 100.0 * 12.0

            if np.sum((rainfall_rate > 200) & (quality_flag > 60) )> 1:
                # Generate unique filename for each plot
                plot_filename = f"rainfall_quality_map_{random_file.split('/')[-1].replace('.nc', '.png')}"

                # Plot and save
                plot_rainfall_and_quality_maps(rainfall_rate, quality_flag, plot_filename)

                print(f"Plot saved as {plot_filename}", flush=True)


    except Exception as e:
        print(f"Error processing {random_file}: {e}", file=sys.stdout, flush=True)
        continue

print("Processing complete.", flush=True)
