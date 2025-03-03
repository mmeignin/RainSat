import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import cartopy.crs as ccrs
import glob
import os
import random
from tqdm import tqdm  # For tracking progress

# File path to NetCDF grid data
file_path = "/net/nfs/precipitations/nicoEUMETSAT/grid_france_1536-1km_200906110000-219912310000.nc"

# Open the NetCDF file to get grid dimensions
with Dataset(file_path) as dataset:
    latitude = dataset.variables["latitude"][:]  # (1536, 1536)
    longitude = dataset.variables["longitude"][:]  # (1536, 1536)

# Initialize a cumulative rain rate array
cumul = np.zeros_like(latitude, dtype=float)

# Directory containing imagette files
imagette_dir = "/net/nfs/precipitations/MatthieuM/second_dataset/*/*.npy"  # <-- Update this path!

# Retrieve all imagette files
all_images = glob.glob(imagette_dir)

# Sort the images by their directory and filename
all_images.sort(key=lambda x: (os.path.dirname(x), os.path.basename(x)))

# Optionally sample a subset of images
# For example, taking half of the images
all_images = random.sample(all_images,len(all_images))

if not all_images:
    print("No imagette files found.")
    exit()

# Determine the update step dynamically based on the total number of files
update_step = max(1, len(all_images) // 100)  # Ensure step is at least 1

# Set up the progress bar
with tqdm(total=len(all_images), desc="Processing Images", unit="file") as pbar:
    for i, path in enumerate(all_images):
        data = np.load(path, allow_pickle=True).item()

        # Extract rain rate and bounding box
        rain_rate = np.array(data['rain_rate'])
        bbox = np.array(data['bbox'])  # Bounding box (row_min, row_max, col_min, col_max)

        r_min, c_min, r_max, c_max = bbox  # Unpack the bounding box
        imagette_mask = np.zeros_like(cumul)

        # Apply rain rate values to the appropriate region
        imagette_mask[r_min:r_max, c_min:c_max] = rain_rate

        # Only accumulate rain rate values above 0.1 mm/hr
        cumul += np.where(imagette_mask > 0.1, imagette_mask, 0)

        # Update progress bar at dynamic intervals
        if (i + 1) % update_step == 0:
            pbar.update(update_step)

    # Ensure the progress bar reaches 100% at the end
    pbar.update(len(all_images) % update_step)


# Set up the figure and projection
plt.figure(figsize=(12, 10))
#projection = ccrs.PlateCarree()
projection = ccrs.Stereographic(central_longitude=0, central_latitude=90)
projection_str = str(projection)
import re
# Extract the projection type
proj_match = re.search(r"\+proj=(\w+)", projection_str)

proj_match = proj_match.group(1)  # Output: stere

ax = plt.axes(projection=projection)

T_cumul  = np.log1p(np.transpose(cumul))
# Plot cumulative precipitation
p = ax.pcolormesh(longitude, latitude, np.where(T_cumul >0,T_cumul,np.nan), cmap="jet", transform=ccrs.PlateCarree())

# Add coastlines and gridlines
ax.coastlines()

# Add title and colorbar
plt.title("Cumulative Rain Rate Over France")
plt.colorbar(p, ax=ax, label="Cumulative Rain Rate (mm)")

# Save and display the plot
plt.savefig(f"cumulative_rain_rate_{proj_match}.png", dpi=300)
plt.close()
