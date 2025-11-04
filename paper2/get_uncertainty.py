import xarray as xr
import pandas as pd
from shapely.geometry import Point
import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon


def process_nc_file(sat_file, unc_file):
    # Load satellite data; change lat and lon to the boundaries of Uganda
    sat_data = xr.open_dataset(sat_file).sel(lat=slice(10, 72), lon=slice(-180., -50))
    
    # Load uncertainty data; change lat and lon to the boundaries of Uganda
    unc_data = xr.open_dataset(unc_file).sel(lat=slice(10, 72), lon=slice(-180., -50))
    
    # Calculate percent uncertainty
    percent = (unc_data['GWRPM25SIGMA'] / sat_data['GWRPM25']) * 100

    # Combine into a new dataset
    filtered = xr.Dataset({
        'unc': percent,
        'lat': sat_data['lat'],
        'lon': sat_data['lon']
    })
    return filtered


# Process each year's files; right now this will create an average from 2021-2023
filtered_20 = process_nc_file(
    '/data/acker/WashU_V5_GL/V5GL0502.HybridPM25.Global.202101-202112.nc',
    '/data/acker/WashU_GL_unc/V5GL0502.HybridPM25E.Global.202101-202112.nc'
)

filtered_21 = process_nc_file(
    '/data/acker/WashU_V5_GL/V5GL0502.HybridPM25.Global.202201-202212.nc',
    '/data/acker/WashU_GL_unc/V5GL0502.HybridPM25E.Global.202201-202212.nc'
)

filtered_22 = process_nc_file(
    '/data/acker/WashU_V5_GL/V5GL0502.HybridPM25.Global.202301-202312.nc',
    '/data/acker/WashU_GL_unc/V5GL0502.HybridPM25E.Global.202301-202312.nc'
)

# Combine filtered DataFrames
combined = xr.concat([filtered_20, filtered_21, filtered_22], dim="year")
pm25_mean = combined['unc'].mean(dim="year")
print(pm25_mean)

pm25_mean.to_netcdf('/data/acker/ALA/paper2/GLv5_uncertainty_mean2021-2023.nc') #change this to your own data directory
print('netcdf done')

var_name = "unc"           # Update with your variable name

ds = xr.open_dataset('/data/acker/ALA/paper2/GLv5_uncertainty_mean2021-2023.nc')
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
    'unc': values,
    'lat': lats,
    'lon': lons,
    'geometry': geoms
}, crs="EPSG:4269")

# Save to shapefile
gdf.to_file("/data/acker/ALA/paper2/GL_uncertainty_2021-2023_grids.shp")

#print(mean_df)

#data_20_21 = pd.merge(filtered_df_20, filtered_df_21, left_index=True, right_index=True, how = 'inner')
#all_data = pd.merge(data_20_21, filtered_df_22, left_index=True, right_index=True, how = 'inner')
#all_data
