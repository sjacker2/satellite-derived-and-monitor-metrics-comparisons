import xarray as xr
import geopandas as gpd
from shapely.geometry import Polygon
import numpy as np

# --- Load NetCDF ---
nc_path = "/data/acker/ALA/paper2/GLv5_GWRPM25_mean2021-2023.nc"  # Update with your NetCDF file path
var_name = "GWRPM25"           # Update with your variable name

ds = xr.open_dataset(nc_path)
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
    'value': values,
    'lat': lats,
    'lon': lons,
    'geometry': geoms
}, crs="EPSG:4326")

# Save to shapefile
gdf.to_file("/data/acker/ALA/paper2/GL_2021-2023_grids.shp")
