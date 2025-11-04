import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
from shapely.geometry import Point

df = gpd.read_file('/data/acker/ALA/washu_2021-2023_unmonitored.shp')
print(len(df))

#county area
df = df.to_crs(epsg=5070)
df['area_km2'] = df.geometry.area / 1e6
df['size'] = df['area_km2'].apply(lambda x: 'large' if x > 5000 else 'small')
df = df.to_crs(epsg=4269)

fires = gpd.read_file('/data/acker/shapefiles/mtbs_perims_DD.shp')
fires = fires[(fires['Ig_Date'] > '2021-01-01 00:00:00') & (fires['Ig_Date'] < '2024-01-01 00:00:00')]
merged_fires = gpd.sjoin(fires, df, how='inner', predicate='intersects')
# Create a column in mountain counties to indicate they are in the mountain region
merged_fires_regions =merged_fires[['GEOID']].copy()
merged_fires_regions['fire_region'] = 'Fire'
# Merge with all counties, marking non-mountain counties
all_counties_gdf = df.merge(merged_fires_regions, on='GEOID', how='left')
# Fill in missing values as 'Non-Mountain'
all_counties_gdf['fire_region'] = all_counties_gdf['fire_region'].fillna('Non-Fire')
all_counties_gdf = all_counties_gdf.drop_duplicates(subset='GEOID')
print(all_counties_gdf.head(5))


mountains = gpd.read_file('/data/acker/shapefiles/ne_10m_geography_regions_polys_ExportFeatures.shp')
mountains = mountains.to_crs(epsg=4269)
merged_mountains = gpd.sjoin(df, mountains, how='inner', predicate='intersects')
# Create a column in mountain counties to indicate they are in the mountain region
merged_mountains =merged_mountains[['GEOID']].copy()
merged_mountains['mountain_region'] = 'Mountain'
# Merge with all counties, marking non-mountain counties
all_counties_gdf = all_counties_gdf.merge(merged_mountains, on='GEOID', how='left')
# Fill in missing values as 'Non-Mountain'
all_counties_gdf['mountain_region'] = all_counties_gdf['mountain_region'].fillna('Non-Mountain')
all_counties_gdf = all_counties_gdf.drop_duplicates(subset=['GEOID'])
print(all_counties_gdf.head(5))


deserts = gpd.read_file('/data/acker/shapefiles/ne_10m_geography_regions_polys_ExportFeatures1.shp')
deserts = deserts.to_crs(epsg=4269)
merged_deserts = gpd.sjoin(df, deserts, how='inner', predicate='intersects')
# Create a column in mountain counties to indicate they are in the mountain region
merged_deserts =merged_deserts[['GEOID']].copy()
merged_deserts['desert_region'] = 'Desert'
# Merge with all counties, marking non-mountain counties
all_counties_gdf = all_counties_gdf.merge(merged_deserts, on='GEOID', how='left')
# Fill in missing values as 'Non-Mountain'
all_counties_gdf['desert_region'] = all_counties_gdf['desert_region'].fillna('Non-Desert')
all_counties_gdf = all_counties_gdf.drop_duplicates(subset=['GEOID'])
print(all_counties_gdf.head(5))


urban = gpd.read_file('/data/acker/shapefiles/tl_rd22_us_uac20.shp')
# Reproject both to equal-area projection before area calculations
all_counties_gdf = all_counties_gdf.to_crs(epsg=5070)
urban = urban.to_crs(epsg=5070)
# Perform geometric intersection: gives only urban areas inside each county
urban_intersection = gpd.overlay(all_counties_gdf, urban, how='intersection')
# Calculate area of the clipped (within-county) urban polygons
urban_intersection['urban_area'] = urban_intersection.geometry.area
# Calculate area of full counties
all_counties_gdf['county_area'] = all_counties_gdf.geometry.area
# Sum urban area per county
urban_area_by_county = (
    urban_intersection
    .groupby('GEOID')['urban_area']
    .sum()
    .reset_index()
)
# Merge back into full counties GeoDataFrame
all_counties_gdf = all_counties_gdf.merge(urban_area_by_county, on='GEOID', how='left')
# Fill missing values for counties with no urban area
all_counties_gdf['urban_area'] = all_counties_gdf['urban_area'].fillna(0)
# Calculate % of county that is urban
all_counties_gdf['urban_pct'] = all_counties_gdf['urban_area'] / all_counties_gdf['county_area']
# (Optional) Reproject back to original CRS if needed
all_counties_gdf = all_counties_gdf.to_crs(epsg=4269)
all_counties_gdf['urban_category'] = all_counties_gdf['urban_pct'].apply(lambda x: 'Urban (≥50%)' if x >= 0.50 else 'Non-Urban (<50%)')


print(all_counties_gdf.columns)
print(len(all_counties_gdf))

all_counties_gdf.to_file('/data/acker/ALA/paper2/all_variables_unmonitored.shp')


