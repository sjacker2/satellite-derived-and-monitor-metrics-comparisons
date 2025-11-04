import xarray as xr
import pandas as pd
from shapely.geometry import Point
import geopandas as gpd

def process_nc_file(sat_file, unc_file):
    # Load satellite data
    sat_data = xr.open_dataset(sat_file).sel(lat=slice(10, 72), lon=slice(-180., -50)).to_dataframe()
    
    # Load uncertainty data
    unc_data = xr.open_dataset(unc_file).sel(lat=slice(10, 72), lon=slice(-180., -50)).to_dataframe()
    
    # Merge and calculate percent uncertainty
    merged = pd.merge(sat_data, unc_data, left_index=True, right_index=True, how='inner')
    merged.reset_index(inplace=True)
    merged['percent'] = (merged['GWRPM25SIGMA'] / merged['GWRPM25']) * 100
    
    # Filter and return
    return merged[merged['percent'] < 90][['lat', 'lon', 'GWRPM25']]

# Process each year's files
filtered_20 = process_nc_file(
    '/data/acker/WashU_V5_GL/V5GL0502.HybridPM25.Global.202201-202212.nc',
    '/data/acker/WashU_GL_unc/V5GL0502.HybridPM25E.Global.202201-202212.nc'
)


print(filtered_20)


#print(filtered_df_22)

# Combine all filtered dataframes, groupby location, find mean, and export as cs



'''print('coverting to shapefile')
mean_df['geometry'] = mean_df.apply(lambda row: Point(row['lon'], row['lat']), axis=1)
df_gpd = gpd.GeoDataFrame(mean_df, geometry='geometry', crs="EPSG:4269")
#mean_df.to_csv('/data/acker/ALA/GLv5_GWRPM25_mean2020-2022.csv', index=False)
df_gpd.to_file('/data/acker/ALA/GL/GLv5_GWRPM25_mean2021-2023.shp')
print('all done')
#print(mean_df)

#data_20_21 = pd.merge(filtered_df_20, filtered_df_21, left_index=True, right_index=True, how = 'inner')
#all_data = pd.merge(data_20_21, filtered_df_22, left_index=True, right_index=True, how = 'inner')
#all_data
'''