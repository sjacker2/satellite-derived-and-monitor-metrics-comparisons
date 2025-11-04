import xarray as xr
import pandas as pd
from shapely.geometry import Point
import geopandas as gpd

def process_nc_file(sat_file, unc_file):
    # Load satellite data; change lat and lon to the boundaries of Uganda
    sat_data = xr.open_dataset(sat_file).sel(lat=slice(10, 72), lon=slice(-180., -50))
    
    # Load uncertainty data; change lat and lon to the boundaries of Uganda
    unc_data = xr.open_dataset(unc_file).sel(lat=slice(10, 72), lon=slice(-180., -50))
    
    # Calculate percent uncertainty
    percent = (unc_data['GWRPM25SIGMA'] / sat_data['GWRPM25']) * 100

    # Mask values with percent uncertainty < 90
    mask = percent < 90
    filtered_pm25 = sat_data['GWRPM25'].where(mask)

    # Combine into a new dataset
    filtered = xr.Dataset({
        'GWRPM25': filtered_pm25,
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
pm25_mean = combined['GWRPM25'].mean(dim="year")
print(pm25_mean)

pm25_mean.to_netcdf('/data/acker/ALA/paper2/GLv5_GWRPM25_mean2021-2023.nc') #change this to your own data directory
print('all done')
#print(mean_df)

#data_20_21 = pd.merge(filtered_df_20, filtered_df_21, left_index=True, right_index=True, how = 'inner')
#all_data = pd.merge(data_20_21, filtered_df_22, left_index=True, right_index=True, how = 'inner')
#all_data
