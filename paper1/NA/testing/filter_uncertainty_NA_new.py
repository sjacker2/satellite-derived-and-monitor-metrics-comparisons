import xarray as xr
import pandas as pd
from shapely.geometry import Point
import geopandas as gpd

def process_nc_file(sat_file):
    # Load satellite data
    sat_data = xr.open_dataset(sat_file).to_dataframe()
    
    
    
    
    # Filter and return
    return sat_data
    
# Process each year's files
filtered_20 = process_nc_file(
    '/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2016001-2016364.nc'
    
)
filtered_21 = process_nc_file(
    '/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2017001-2017364.nc'
)
filtered_22 = process_nc_file(
    '/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2018001-2018364.nc'
    )

# Combine filtered DataFrames
combined_df = pd.concat([filtered_20, filtered_21, filtered_22])
print(combined_df)
combined_df.reset_index(inplace=True)
# Combine all filtered dataframes, group by location, and find mean
print('combined')
mean_df = combined_df.groupby(['lat', 'lon'])['GWRPM25'].mean().reset_index()
mean_df = mean_df.dropna()
print(mean_df)
print(mean_df['GWRPM25'].max())
# Use vectorized operations for geometry creation
print('converting to shapefile')
mean_df['geometry'] = gpd.points_from_xy(mean_df['lon'], mean_df['lat'])
mean_df.dropna(inplace=True)
# Convert to GeoDataFrame
df_gpd = gpd.GeoDataFrame(mean_df, geometry='geometry', crs="EPSG:4269")
df_gpd.to_file('/data/acker/ALA/NA/NAv5_GWRPM25_mean2016-2018.shp')
print('all done')
