import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import cartopy.crs as ccrs
import glob
import itertools
import os
from matplotlib.colors import LogNorm  

# File path to base NetCDF data
file_path = '/net/nfs/precipitations/nicoEUMETSAT/grid_france_1536-1km_200906110000-219912310000.nc'

# Open the NetCDF file
with Dataset(file_path) as dataset:
    latitude = dataset.variables['latitude'][:]  
    longitude = dataset.variables['longitude'][:]  

# Define years, months, and days to analyze
years = ["2018"]
months = range(1, 13)
days = range(1, 32)

# List to hold all matching files
all_matching_files = []

# Find matching NetCDF files
for year, month, day in itertools.product(years, months, days):
    file_pattern = f"/net/nfs/precipitations/SEDOO/Mosaiques/LMDO_1536_1km_5min/{year}/{year}{month:02d}/{year}{month:02d}{day:02d}/cumul_france_1536-1km-5min_{year}{month:02d}{day:02d}*.nc"
    all_matching_files.extend(glob.glob(file_pattern))

all_matching_files.sort()

if not all_matching_files:
    print("No matching files found.")
    exit()

print(f"Found {len(all_matching_files)} matching files.")

# Use the first matching file
file_path = all_matching_files[0]

# Open and process the selected NetCDF file
with Dataset(file_path) as dataset:
    print(dataset)
    
    # Extract rain rate and quality data
    LRR = dataset['ACRR'][0, :, :].astype(float)
    QRR = dataset['QUALITY'][0, :, :].astype(int)
    
    # Apply data corrections
    LRR[(LRR > 65000) | (QRR > 200)] = np.nan  # Ignore invalid values
    LRR[QRR < 60] = np.nan
    LRR *= 12 / 100.0  # Convert to mm/hr

    # Mask low values
    rain_rate = np.where(LRR >= 0.1, LRR, np.nan)

# Create a figure
plt.figure(figsize=(10, 8))
projection = ccrs.PlateCarree()
ax = plt.axes(projection=projection)

# Plot the data with logarithmic normalization
p = ax.pcolormesh(longitude, latitude, rain_rate, cmap='jet', shading="auto",
                   transform=ccrs.PlateCarree(), norm=LogNorm(vmin=0.1, vmax=50))

# Add map features
ax.coastlines()

# ✅ **Proper gridlines instead of `ax.grid()`**
gl = ax.gridlines(draw_labels=True, linestyle="--", linewidth=0.5, color="gray")
gl.top_labels = False  # Remove top labels
gl.right_labels = False  # Remove right labels

plt.title('Rain Rate Over France (Log Scale)')
plt.colorbar(p, ax=ax, label='Rain Rate (mm/hr)')

# Save the plot
file_name = os.path.splitext(os.path.basename(file_path))[0]
plt.savefig(f"rain_rate_grid_{file_name}.png", dpi=300)

# Show the plot
plt.show()
