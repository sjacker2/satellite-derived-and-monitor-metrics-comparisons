import xarray as xr
import pandas as pd
from shapely.geometry import Point
import geopandas as gpd

def process_nc_file(sat_file, unc_file, g_file):
    # Load satellite data
    sat_data = xr.open_dataset(sat_file).to_dataframe()
    
    # Load uncertainty data
    unc_data = xr.open_dataset(unc_file).sel(lat=slice(14, 68), lon=slice(-170., -52)).to_dataframe()
    gl_data = xr.open_dataset(g_file).sel(lat=slice(14, 68), lon=slice(-170., -52)).to_dataframe()

    sat_data = sat_data.reset_index()
    
    
    # Merge and calculate percent uncertainty
    merged = pd.merge(gl_data, unc_data, left_index=True, right_index=True, how='inner')
    merged.reset_index(inplace=True)
    
    # Calculate the 90th percentile globally
    merged['percent'] = (merged['GWRPM25SIGMA'] / merged['GWRPM25']) * 100
    #merged = merged.rename(columns={'lat_x':'lat','lon_x':'lon'})
    print(merged[merged['percent'] >= 90][['lat', 'lon', 'GWRPM25']])
    remove = merged[merged['percent'] >= 90][['lat', 'lon', 'GWRPM25']]
    filtered_n = sat_data.merge(remove[['lat', 'lon']], on=['lat', 'lon'], how='left', indicator=True)

    # Keep only rows that do not have a match in 'remove'
    filtered_n = filtered_n[filtered_n['_merge'] == 'left_only']

    # Drop the '_merge' column, as it is no longer needed
    filtered_n = filtered_n.drop(columns=['_merge'])    # Filter and return
    print(filtered_n)
    return filtered_n
    
# Process each year's files
filtered_20 = process_nc_file(
    '/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2015001-2015364.nc',
    '/data/acker/WashU_GL_unc/V5GL04.HybridPM25E.Global.201501-201512.nc',
    '/data/acker/WashU_V5_GL/V5GL04.HybridPM25.Global.201501-201512.nc'
)
filtered_21 = process_nc_file(
    '/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2016001-2016364.nc',
    '/data/acker/WashU_GL_unc/V5GL04.HybridPM25E.Global.201601-201612.nc',
    '/data/acker/WashU_V5_GL/V5GL04.HybridPM25.Global.201601-201612.nc'
)
filtered_22 = process_nc_file(
    '/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2017001-2017364.nc',
    '/data/acker/WashU_GL_unc/V5GL04.HybridPM25E.Global.201701-201712.nc',
    '/data/acker/WashU_V5_GL/V5GL04.HybridPM25.Global.201701-201712.nc'
)

# Combine filtered DataFrames
combined_df = pd.concat([filtered_20, filtered_21, filtered_22])

# Combine all filtered dataframes, group by location, and find mean
print('combined')
mean_df = combined_df.groupby(['lat', 'lon'])['GWRPM25'].mean().reset_index()
mean_df = mean_df.dropna()
print(mean_df)
print(mean_df['GWRPM25'].max())
# Use vectorized operations for geometry creation
print('converting to shapefile')
mean_df['geometry'] = gpd.points_from_xy(mean_df['lon'], mean_df['lat'])

# Convert to GeoDataFrame
df_gpd = gpd.GeoDataFrame(mean_df, geometry='geometry', crs="EPSG:4269")
df_gpd.to_file('/data/acker/ALA/NA/NAv5_GWRPM25_mean2015-2017_test_GL.shp')
print('all done')
