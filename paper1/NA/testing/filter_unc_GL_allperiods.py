import xarray as xr
import pandas as pd
from shapely.geometry import Point
import geopandas as gpd

def process_nc_file(sat_file, unc_file, g_file, remove_pixels=None):
    # Load satellite data
    sat_data = xr.open_dataset(sat_file).to_dataframe()
    
    # Load uncertainty data
    unc_data = xr.open_dataset(unc_file).sel(lat=slice(14, 68), lon=slice(-170., -52)).to_dataframe()
    gl_data = xr.open_dataset(g_file).sel(lat=slice(14, 68), lon=slice(-170., -52)).to_dataframe()

    sat_data = sat_data.reset_index()
    
    # Merge and calculate percent uncertainty
    merged = pd.merge(gl_data, unc_data, left_index=True, right_index=True, how='inner')
    merged.reset_index(inplace=True)
    merged['percent'] = (merged['GWRPM25SIGMA'] / merged['GWRPM25']) * 100
    
    # Identify pixels above the 90% uncertainty threshold
    remove = merged[merged['percent'] >= 90][['lat', 'lon']]
    
    # Update remove_pixels tracking
    if remove_pixels is not None:
        remove_pixels.extend(remove[['lat', 'lon']].values.tolist())
    
    # Filter satellite data to exclude removed pixels
    filtered = sat_data.merge(remove, on=['lat', 'lon'], how='left', indicator=True)
    filtered = filtered[filtered['_merge'] == 'left_only'].drop(columns=['_merge'])
    return filtered

def process_time_period(time_period):
    label, *file_groups = time_period
    remove_pixels = []
    yearly_data = []
    #print(file_groups)
    # Process each year in the time period
    for year_files in file_groups:
        sat_file, unc_file, g_file = year_files
        filtered = process_nc_file(sat_file, unc_file, g_file, remove_pixels=remove_pixels)
        yearly_data.append(filtered)

    # Combine yearly data, calculate mean across years
    combined_df = pd.concat(yearly_data)
    averaged_df = combined_df.groupby(['lat', 'lon'])['GWRPM25'].mean().reset_index()
    
    return label, averaged_df, remove_pixels

# Define time periods and file paths
# Define all time periods with corresponding file paths
time_periods = [
    ("2015-2017", 
     [
         '/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2015001-2015364.nc',
         '/data/acker/WashU_GL_unc/V5GL04.HybridPM25E.Global.201501-201512.nc',
         '/data/acker/WashU_V5_GL/V5GL04.HybridPM25.Global.201501-201512.nc'
     ],
     [
         '/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2016001-2016364.nc',
         '/data/acker/WashU_GL_unc/V5GL04.HybridPM25E.Global.201601-201612.nc',
         '/data/acker/WashU_V5_GL/V5GL04.HybridPM25.Global.201601-201612.nc'
     ],
     [
         '/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2017001-2017364.nc',
         '/data/acker/WashU_GL_unc/V5GL04.HybridPM25E.Global.201701-201712.nc',
         '/data/acker/WashU_V5_GL/V5GL04.HybridPM25.Global.201701-201712.nc'
     ]),
    ("2016-2018", 
     [
         '/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2016001-2016364.nc',
         '/data/acker/WashU_GL_unc/V5GL04.HybridPM25E.Global.201601-201612.nc',
         '/data/acker/WashU_V5_GL/V5GL04.HybridPM25.Global.201601-201612.nc'
     ],
     [
         '/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2017001-2017364.nc',
         '/data/acker/WashU_GL_unc/V5GL04.HybridPM25E.Global.201701-201712.nc',
         '/data/acker/WashU_V5_GL/V5GL04.HybridPM25.Global.201701-201712.nc'
     ],
     [
         '/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2018001-2018364.nc',
         '/data/acker/WashU_GL_unc/V5GL04.HybridPM25E.Global.201801-201812.nc',
         '/data/acker/WashU_V5_GL/V5GL04.HybridPM25.Global.201801-201812.nc'
     ]),
    ("2017-2019", 
     [
         '/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2017001-2017364.nc',
         '/data/acker/WashU_GL_unc/V5GL04.HybridPM25E.Global.201701-201712.nc',
         '/data/acker/WashU_V5_GL/V5GL04.HybridPM25.Global.201701-201712.nc'
     ],
     [
         '/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2018001-2018364.nc',
         '/data/acker/WashU_GL_unc/V5GL04.HybridPM25E.Global.201801-201812.nc',
         '/data/acker/WashU_V5_GL/V5GL04.HybridPM25.Global.201801-201812.nc'
     ],
     [
         '/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2019001-2019364.nc',
         '/data/acker/WashU_GL_unc/V5GL04.HybridPM25E.Global.201901-201912.nc',
         '/data/acker/WashU_V5_GL/V5GL04.HybridPM25.Global.201901-201912.nc'
     ]),
    ("2018-2020", 
     [
         '/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2018001-2018364.nc',
         '/data/acker/WashU_GL_unc/V5GL04.HybridPM25E.Global.201801-201812.nc',
         '/data/acker/WashU_V5_GL/V5GL04.HybridPM25.Global.201801-201812.nc'
     ],
     [
         '/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2019001-2019364.nc',
         '/data/acker/WashU_GL_unc/V5GL04.HybridPM25E.Global.201901-201912.nc',
         '/data/acker/WashU_V5_GL/V5GL04.HybridPM25.Global.201901-201912.nc'
     ],
     [
         '/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2020001-2020364.nc',
         '/data/acker/WashU_GL_unc/V5GL04.HybridPM25E.Global.202001-202012.nc',
         '/data/acker/WashU_V5_GL/V5GL04.HybridPM25.Global.202001-202012.nc'
     ]),
    ("2019-2021", 
     [
         '/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2019001-2019364.nc',
         '/data/acker/WashU_GL_unc/V5GL04.HybridPM25E.Global.201901-201912.nc',
         '/data/acker/WashU_V5_GL/V5GL04.HybridPM25.Global.201901-201912.nc'
     ],
     [
         '/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2020001-2020364.nc',
         '/data/acker/WashU_GL_unc/V5GL04.HybridPM25E.Global.202001-202012.nc',
         '/data/acker/WashU_V5_GL/V5GL04.HybridPM25.Global.202001-202012.nc'
     ],
     [
         '/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2021001-2021364.nc',
         '/data/acker/WashU_GL_unc/V5GL04.HybridPM25E.Global.202101-202112.nc',
         '/data/acker/WashU_V5_GL/V5GL04.HybridPM25.Global.202101-202112.nc'
     ]),
    ("2020-2022", 
     [
         '/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2020001-2020364.nc',
         '/data/acker/WashU_GL_unc/V5GL04.HybridPM25E.Global.202001-202012.nc',
         '/data/acker/WashU_V5_GL/V5GL04.HybridPM25.Global.202001-202012.nc'
     ],
     [
         '/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2021001-2021364.nc',
         '/data/acker/WashU_GL_unc/V5GL04.HybridPM25E.Global.202101-202112.nc',
         '/data/acker/WashU_V5_GL/V5GL04.HybridPM25.Global.202101-202112.nc'
     ],
     [
         '/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2022001-2022364.nc',
         '/data/acker/WashU_GL_unc/V5GL04.HybridPM25E.Global.202201-202212.nc',
         '/data/acker/WashU_V5_GL/V5GL04.HybridPM25.Global.202201-202212.nc'
     ])
]


# Track removed pixels across all time periods
all_remove_pixels = []
all_outputs = []

# Process each time period
for time_period in time_periods:
    label, averaged_df, remove_pixels = process_time_period(time_period)
    print(f"Processed time period: {label} - {len(averaged_df)} rows after averaging")
    all_outputs.append((label, averaged_df))
    all_remove_pixels.extend(remove_pixels)  # Store all removed pixels

# Debug: How many pixels were flagged for removal?
print(f"Total removed pixels collected: {len(all_remove_pixels)}")

# Determine pixels removed > 25% of the time
remove_counts = pd.DataFrame(all_remove_pixels, columns=['lat', 'lon'])
remove_counts = remove_counts.value_counts().reset_index(name='count')

# Set threshold for frequent removals
threshold = len(time_periods) * 0.25
frequently_removed = remove_counts[remove_counts['count'] > threshold][['lat', 'lon']]

# Debug: Check frequently removed pixels
print(f"Number of pixels to be removed from all final outputs: {len(frequently_removed)}")
print(frequently_removed.head())

# Finalize outputs for each time period
final_outputs = []
for label, averaged_df in all_outputs:
    pre_filter_size = len(averaged_df)

    # Merge to identify rows to remove
    filtered_final = averaged_df.merge(frequently_removed, on=['lat', 'lon'], how='left', indicator=True)
    filtered_final = filtered_final[filtered_final['_merge'] == 'left_only'].drop(columns=['_merge'])

    post_filter_size = len(filtered_final)
    print(f"{label}: Before filtering: {pre_filter_size}, After filtering: {post_filter_size}")
    filtered_final.dropna(inplace=True)
    print('drop:',len(filtered_final))
    final_outputs.append((label, filtered_final))

# Save outputs as shapefiles
for label, filtered_final in final_outputs:
    filtered_final['geometry'] = gpd.points_from_xy(filtered_final['lon'], filtered_final['lat'])
    df_gpd = gpd.GeoDataFrame(filtered_final, geometry='geometry', crs="EPSG:4269")

    output_path = f'/data/acker/ALA/NA/NAv5_GWRPM25_mean_{label}_exact_GL.shp'
    df_gpd.to_file(output_path)
    print(f"Shapefile saved: {output_path}")
    break
print('done')
