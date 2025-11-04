import xarray as xr
import pandas as pd
import numpy as np
from shapely.geometry import Point
import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon


def get_removed_pixels(sat_file, unc_file):
    # Load satellite and uncertainty data for Uganda region
    sat_data = xr.open_dataset(sat_file).sel(lat=slice(10, 72), lon=slice(-180., -50))
    unc_data = xr.open_dataset(unc_file).sel(lat=slice(10, 72), lon=slice(-180., -50))
    
    # Calculate percent uncertainty
    percent = (unc_data['GWRPM25SIGMA'] / sat_data['GWRPM25']) * 100
    
    # Mask for removed pixels: percent uncertainty >= 90
    removed_mask = percent >= 90
    removed_percent = percent.where(removed_mask)
    
    return removed_percent

# Get removed pixels for each year
removed_20 = get_removed_pixels(
    '/data/acker/WashU_V5_GL/V5GL0502.HybridPM25.Global.202101-202112.nc',
    '/data/acker/WashU_GL_unc/V5GL0502.HybridPM25E.Global.202101-202112.nc'
)

removed_21 = get_removed_pixels(
    '/data/acker/WashU_V5_GL/V5GL0502.HybridPM25.Global.202201-202212.nc',
    '/data/acker/WashU_GL_unc/V5GL0502.HybridPM25E.Global.202201-202212.nc'
)

removed_22 = get_removed_pixels(
    '/data/acker/WashU_V5_GL/V5GL0502.HybridPM25.Global.202301-202312.nc',
    '/data/acker/WashU_GL_unc/V5GL0502.HybridPM25E.Global.202301-202312.nc'
)

# Stack along a new dimension 'year'
removed_combined = xr.concat([removed_20, removed_21, removed_22], dim='year')
removed_combined['year'] = [2021, 2022, 2023]  # assign actual years as coordinates

# Count how many years each pixel was removed
valid_counts = (~np.isnan(removed_combined)).sum(dim='year')

# Compute average percent uncertainty for pixels removed in 1+ years
mean_removed_percent = removed_combined.mean(dim='year', skipna=True)

# Mask to keep only pixels that were removed in at least one year
final_removed = mean_removed_percent.where(valid_counts > 0)

# Save result to NetCDF
final_removed.name = "Removed_Percent_Uncertainty"
final_removed.to_netcdf('/data/acker/ALA/paper2/GLv5_removed_pixels_2021-2023.nc')

print("Saved NetCDF of removed pixels with average percent uncertainty.")

var_name = "Removed_Percent_Uncertainty"           # Update with your variable name

ds = xr.open_dataset('/data/acker/ALA/paper2/GLv5_removed_pixels_2021-2023.nc')
data = ds[var_name]

# Select first time slice if needed
if "time" in data.dims:
    data = data.isel(time=0)

# Get lat/lon
lat = ds['lat'].values
lon = ds['lon'].values

# Handle 1D or 2D lat/lon
if lat.ndim == 1 and lon.ndim == 1:
    lat2d, lon2d = np.meshgrid(lat, lon, indexing='ij')
else:
    lat2d, lon2d = lat, lon

# Initialize geometry and attribute lists
geoms = []
values = []
lats = []
lons = []

nrows, ncols = data.shape
for i in range(nrows - 1):
    print(i)
    for j in range(ncols - 1):
        corners = [
            (lon2d[i, j],     lat2d[i, j]),
            (lon2d[i, j+1],   lat2d[i, j+1]),
            (lon2d[i+1, j+1], lat2d[i+1, j+1]),
            (lon2d[i+1, j],   lat2d[i+1, j])
        ]
        poly = Polygon(corners)
        val = data.values[i, j]

        if not np.isnan(val):  # Only keep valid cells
            geoms.append(poly)
            values.append(val)
            lats.append(lat2d[i, j])
            lons.append(lon2d[i, j])

# Create GeoDataFrame with lat/lon columns
gdf = gpd.GeoDataFrame({
    'removed': values,
    'lat': lats,
    'lon': lons,
    'geometry': geoms
}, crs="EPSG:4269")

# Save to shapefile
gdf.to_file("/data/acker/ALA/paper2/GL_removed_pixels_2021-2023_grids.shp")

