import xarray as xr
import pandas as pd
from shapely.geometry import Point
import geopandas as gpd
import dask.dataframe as dd

def process_nc_file(sat_file):
    # Load satellite data
    sat_data = xr.open_dataset(sat_file).sel(lat=slice(10, 72), lon=slice(-180., -50)).to_dataframe()
        
    # Reset index and drop NaN values
    sat_data.reset_index(inplace=True)
    sat_data = sat_data.dropna(subset=['lat', 'lon', 'GWRPM25'])
    
    # Return relevant columns
    return sat_data[['lat', 'lon', 'GWRPM25']]

# Process NetCDF files into pandas DataFrames
filtered_20 = process_nc_file('/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2020001-2020364.nc')
filtered_21 = process_nc_file('/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2021001-2021364.nc')
filtered_22 = process_nc_file('/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2022001-2022364.nc')

# Convert pandas DataFrames to Dask DataFrames
filtered_20_dd = dd.from_pandas(filtered_20, npartitions=10)
filtered_21_dd = dd.from_pandas(filtered_21, npartitions=10)
filtered_22_dd = dd.from_pandas(filtered_22, npartitions=10)

# Combine all the datasets
combined_dd = dd.concat([filtered_20_dd, filtered_21_dd, filtered_22_dd])

# Perform the groupby operation
mean_dd = combined_dd.groupby(['lat', 'lon'])['GWRPM25'].mean()

# Compute the result and convert to a pandas DataFrame
mean_df = mean_dd.compute().reset_index()

# Convert to GeoDataFrame and save as a shapefile
print('Converting to shapefile...')
mean_df['geometry'] = mean_df.apply(lambda row: Point(row['lon'], row['lat']), axis=1)
df_gpd = gpd.GeoDataFrame(mean_df, geometry='geometry', crs="EPSG:4269")
df_gpd.to_file('/data/acker/ALA/NA/NAv5_GWRPM25_mean2020-2022.shp')
print('All done!')
