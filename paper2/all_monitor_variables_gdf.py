import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
from shapely.geometry import Point

# Load county-level shapefile with PM2.5 and CDV
df = gpd.read_file('/data/acker/ALA/paper2/washu_2021-2023_merged.shp')
merged_df = df.rename(columns={'Design Val': 'Design Value'})

# Calculate differences and classification
merged_df['diff'] = merged_df['PM25_90th'] - merged_df['Design Value']
merged_df['abs_diff'] = abs(merged_df['diff'])

threshold = 9.0
def classify(row):
    cdv = row['Design Value']
    cdve = row['PM25_90th']
    if cdv > threshold and cdve > threshold:
        return 'TP'
    elif cdv > threshold and cdve <= threshold:
        return 'FN'
    elif cdv <= threshold and cdve > threshold:
        return 'FP'
    else:
        return 'TN'

#county area
merged_df['classification'] = merged_df.apply(classify, axis=1)
merged_df = merged_df.to_crs(epsg=5070)
merged_df['area_km2'] = merged_df.geometry.area / 1e6
merged_df = merged_df.to_crs(epsg=4269)
merged_df['size'] = merged_df['area_km2'].apply(lambda x: 'large' if x > 5000 else 'small')


# Load DV monitor data
dvs = pd.read_csv('/data/acker/EPA_DV/site_DVs_2021-2023.csv')
dvs.columns = dvs.columns.str.strip()  # Remove trailing whitespace from column names
valid = dvs.dropna(subset=['Valid DV'])
valid['geometry'] = gpd.points_from_xy(valid['Site Longitude'], valid['Site Latitude'])
valid = gpd.GeoDataFrame(valid, geometry='geometry', crs='EPSG:4269')

# Spatial join to assign monitors to counties
merged_df_2 = gpd.sjoin(valid, merged_df, how="left", predicate="within")

# Count monitors per county
monitor_counts = merged_df_2.groupby("GEOID").size().reset_index(name='monitor_count')
counties_with_monitor_counts = merged_df.merge(monitor_counts, on="GEOID", how="left")
counties_with_monitor_counts['monitor_count'] = counties_with_monitor_counts['monitor_count'].fillna(1).astype(int)

# Bin design values
bins = [0, 7, 10, float('inf')]
labels = ['<7', '7–10', '>10']
counties_with_monitor_counts['cdv_bin'] = pd.cut(counties_with_monitor_counts['Design Value'], bins=bins, labels=labels, right=False)

# Load grid-level data
grid = gpd.read_file('/data/acker/ALA/paper2/GL_2021-2023_grids.shp')
grid = grid.to_crs(epsg=4269)
joined = gpd.sjoin(grid, counties_with_monitor_counts, how='left', predicate='intersects')
grids = joined.dropna(subset=['value', 'GEOID'])
grids = grids.drop(columns=['index_right0'])

# Calculate monitor coverage
monitors_with_grid = gpd.sjoin(valid, grids, how="inner", predicate="within")
total_grids_per_county = grids.groupby('GEOID').size().reset_index(name='total_grids')
grids_with_monitor = monitors_with_grid.drop_duplicates(subset='geometry')
monitor_grids_per_county = grids_with_monitor.groupby('GEOID').size().reset_index(name='grids_with_monitor')

coverage_df = total_grids_per_county.merge(monitor_grids_per_county, on='GEOID', how='left')
coverage_df['grids_with_monitor'] = coverage_df['grids_with_monitor'].fillna(0)
coverage_df['monitor_coverage_pct'] = 100 * coverage_df['grids_with_monitor'] / coverage_df['total_grids']

# Merge all county-level variables
final_df = counties_with_monitor_counts.merge(coverage_df, on='GEOID', how='left')

# --- Calculate dist_km between max monitor and max grid cell per county ---

# Ensure all layers are in a projected CRS for distance
merged_df_2 = merged_df_2.to_crs(epsg=5070)
grids = grids.to_crs(epsg=5070)
final_df = final_df.to_crs(epsg=5070)

# Highest DV monitor location per county
highest_monitor = (
    merged_df_2.sort_values('Valid DV', ascending=False)
    .drop_duplicates('GEOID')
    .rename(columns={'geometry': 'monitor_geom'})
    [['GEOID', 'Valid DV', 'monitor_geom']]
)

# Highest PM2.5 grid cell per county
highest_grid = (
    grids.sort_values('PM25_90th', ascending=False)
    .drop_duplicates('GEOID')
    .rename(columns={'geometry': 'grid_geom'})
    [['GEOID', 'PM25_90th', 'grid_geom']]
)

# Merge and compute distance
comparison_df = highest_monitor.merge(highest_grid, on='GEOID', how='inner')
comparison_df['dist_km'] = comparison_df.apply(
    lambda row: row['monitor_geom'].distance(row['grid_geom']) / 1000, axis=1
)

# --- Final shapefile with everything merged ---
print(final_df.columns)
print(len(final_df))

# Merge in diff/classification/abs_diff
'''full_df = final_df.merge(
    merged_df[['GEOID', 'Design Value', 'classification', 'diff', 'abs_diff']],
    on='GEOID', how='left'
)
print(full_df.columns)
print(len(full_df))'''
# Merge in dist_km and geometries
full_df = final_df.merge(
    comparison_df[['GEOID', 'dist_km']],
    on='GEOID', how='left'
)
print(full_df.columns)

# Set geometry to monitor location and clean up
full_df = full_df.set_geometry('geometry')
full_df = gpd.GeoDataFrame(full_df, geometry='geometry', crs='EPSG:4269')
print(full_df.columns)

# Select columns to include in output
cols_to_keep = [
    'GEOID', 'PM25_90th', 'Design Value', 'classification',
    'diff', 'abs_diff', 'monitor_count', 'cdv_bin',
    'monitor_coverage_pct', 'dist_km', 'geometry', 'size'
]
output_df = full_df[cols_to_keep]

# Write to shapefile
output_df.to_file('/data/acker/ALA/paper2/explain_variables_1.shp')
