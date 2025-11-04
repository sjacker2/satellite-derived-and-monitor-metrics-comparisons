import xarray as xr
import pandas as pd
from shapely.geometry import Point
import geopandas as gpd

def process_nc_file(sat_file, unc_file):
    # Load satellite and uncertainty data
    sat_data = xr.open_dataset(sat_file).sel(lat=slice(10, 72), lon=slice(-180., -50))
    unc_data = xr.open_dataset(unc_file).sel(lat=slice(10, 72), lon=slice(-180., -50))
    
    # Calculate the 90th percentile of uncertainty
    percentile_90 = unc_data['GWRPM25SIGMAFULL'].quantile(0.9, dim=['lat', 'lon'])
    
    # Mask sat_data wherever unc_data exceeds the 90th percentile
    sat_data_filtered = sat_data.where(unc_data['GWRPM25SIGMAFULL'] <= percentile_90)
    
    # Return filtered satellite data
    return sat_data

# Process NetCDF files
filtered_20 = process_nc_file(
    '/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2020001-2020364.nc',
    '/data/acker/WashU_NA_unc/V5NA04.02.HybridPM25EFull.xNorthAmerica.2020001-2020364.nc'
)
'''filtered_21 = process_nc_file(
    '/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2021001-2021364.nc',
    '/data/acker/WashU_NA_unc/V5NA04.02.HybridPM25EFull.xNorthAmerica.2016001-2016364.nc'
)
filtered_22 = process_nc_file(
    '/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2022001-2022364.nc',
    '/data/acker/WashU_NA_unc/V5NA04.02.HybridPM25EFull.xNorthAmerica.2017001-2017364.nc'
)'''

# Stack the datasets into a single xarray object
#stacked = xr.concat([filtered_20, filtered_21, filtered_22], dim="dataset")

# Take the mean across the "dataset" dimension
#averaged_ds = stacked.mean(dim="dataset")

# Convert to pandas DataFrame
mean_df = filtered_20.to_dataframe().reset_index()
mean_df = mean_df.dropna()

# Convert to GeoDataFrame and save as a shapefile
print('Converting to shapefile...')
mean_df['geometry'] = mean_df.apply(lambda row: Point(row['lon'], row['lat']), axis=1)
df_gpd = gpd.GeoDataFrame(mean_df, geometry='geometry', crs="EPSG:4269")
df_gpd.to_file('/data/acker/ALA/NA/NAv5_GWRPM25_mean2020_test.shp')
print('All done!')
