import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from shapely.geometry import Point
df = pd.read_csv('/home/kysela/ALA/GLv5_meanGWRPM25(2020-2022).csv')


# change to Point geometry
df['geometry'] = df.apply(lambda row: Point(row['lon'], row['lat']), axis=1)
df_gpd = gpd.GeoDataFrame(df, geometry='geometry', crs="EPSG:4269")

#upload shapefile of U.S. counties
counties = gpd.read_file('/data/acker/shapefiles/tl_2020_us_county.shp')
counties = counties.drop(['STATEFP', 'COUNTYFP', "COUNTYNS", 'NAMELSAD', 'LSAD', 'CLASSFP', 'MTFCC', 'CSAFP', 'CBSAFP', 'METDIVFP', 'FUNCSTAT', 'ALAND', 'AWATER', 'INTPTLAT', 'INTPTLON'], axis=1)
df_gpd = df_gpd.to_crs(counties.crs)

import geopandas as gpd
from shapely.geometry import Point

# Initialize a new column for GEOID
df_gpd['GEOID'] = None

# Iterate through each point and check which county it falls into
for i, point in df_gpd.iterrows():
    for j, poly in counties.iterrows():
        # Check if the point is within the county polygon
        if point['geometry'].within(poly['geometry']):
            df_gpd.at[i, 'GEOID'] = poly['GEOID']
            print(poly['GEOID'])
            break  # Stop once we find the first matching county to avoid unnecessary checks

df_gpd.to_file('/data/acker/ALA/test.shp')