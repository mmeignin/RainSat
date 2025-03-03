import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# File path to your NetCDF data
file_path = '/net/nfs/precipitations/nicoEUMETSAT/grid_france_1536-1km_200906110000-219912310000.nc'

# Open the NetCDF file
dataset = Dataset(file_path)


# Extract latitude and longitude data (2D arrays)
latitude = dataset.variables['latitude'][:]  # Latitude (1536, 1536)
longitude = dataset.variables['longitude'][:]  # Longitude (1536, 1536)

# Check the shape of latitude and longitude to ensure they match (1536 x 1536 grid)
print(f"Latitude shape: {latitude.shape}")
print(f"Longitude shape: {longitude.shape}")

# Create a Stereographic projection using Cartopy's built-in Stereographic class
projection = ccrs.Stereographic(central_longitude=0, central_latitude=90)

# Plot the latitude and longitude grid with coastlines
plt.figure(figsize=(12, 10))

# Create a map projection (using Cartopy for better map support)
ax = plt.axes(projection=projection)

# Plot the grid using pcolormesh, overlay coastlines over it
p = ax.pcolormesh(longitude, latitude, np.ones_like(latitude), cmap='viridis', shading='auto', transform=ccrs.PlateCarree())

# Add coastlines on top of the grid, using the correct stereographic projection
ax.coastlines()

# Add title and labels
plt.title('Coastlines Over the Latitude and Longitude Grid (1536x1536) in Stereographic Projection')
plt.xlabel('Longitude')
plt.ylabel('Latitude')

# Add a colorbar to show the grid structure (it will just show the grid values since we are plotting 1s)
plt.colorbar(p, ax=ax, label='Grid Structure')
plt.grid()
# Save the plot to a file
plt.savefig("france_grid_with_coastlines_stereographic.png")

# Show the plot
plt.show()

# Close the dataset
dataset.close()
