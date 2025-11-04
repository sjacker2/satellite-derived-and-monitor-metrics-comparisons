#!/home/harkey/miniconda3/envs/workenv/bin/python

# clip_write_WUSTLPM25_2shp : read in a Lizzy-made csv of QA'd WUSTL PM2.5 data 
# select data only for the area of interest, create polygons for each grid cell, 
# save as a shapefile for later use
# borrowing/learning from: 
# https://automating-gis-processes.github.io/CSC/notebooks/L2/geopandas-basics.html

import xarray as xr
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, Polygon
from fiona.crs import from_epsg

def process_nc_file(sat_file):
    # Load satellite and uncertainty data
    sat_data = xr.open_dataset(sat_file).sel(lat=slice(9, 73), lon=slice(-180., -49)).to_dataframe()
    # Reset index and drop NaN values
    sat_data.reset_index(inplace=True)
    sat_data = sat_data.dropna(subset=['lat', 'lon', 'GWRPM25'])
    
    # Return relevant columns
    return sat_data[['lat', 'lon', 'GWRPM25']]
    
    
   

# Process NetCDF files
filtered_20 = process_nc_file('/data/acker/WashU_NA_notfull/V5NA04.02.HybridPM25.NorthAmerica.2015001-2015364.nc')
filtered_21 = process_nc_file('/data/acker/WashU_NA_notfull/V5NA04.02.HybridPM25.NorthAmerica.2016001-2016364.nc')
filtered_22 = process_nc_file('/data/acker/WashU_NA_notfull/V5NA04.02.HybridPM25.NorthAmerica.2017001-2017364.nc')

# Combine all the datasets
combined_df = pd.concat([filtered_20, filtered_21, filtered_22])
print('combined')

# Perform the groupby operation
mean_df = combined_df.groupby(['lat', 'lon'])['GWRPM25'].mean().reset_index()
print('mean')

# read in file
df = mean_df
print(df)

# slice data for area of interest
'''minlat = 10
minlon = -179.9
maxlat = 72
maxlon = -50'''
minlat = 40.73
minlon = -112.6
maxlat = 41.19
maxlon = -111.5
#41.19, -112.6, 40.73, -111.5
df = df[(df.lat > minlat) & (df.lat < maxlat) & (df.lon < maxlon) & (df.lon > minlon)].reset_index().drop(columns=['index'])
print(df)

# calculate grid corners from center lat/lon, assuming a horizontal resolution of 0.01 degree
df['LL_lat'] = df.lat - 0.005
df['LL_lon'] = df.lon - 0.005
df['UL_lat'] = df.lat + 0.005
df['UL_lon'] = df.lon - 0.005
df['UR_lat'] = df.lat + 0.005
df['UR_lon'] = df.lon + 0.005
df['LR_lat'] = df.lat - 0.005
df['LR_lon'] = df.lon + 0.005
print(df)

# create geo data frame
newdata = gpd.GeoDataFrame()
newdata['geometry'] = None

index = 0
# loop over rows in dataframe to make polygons to put in the geometry column of newdata
print("looping over rows in input CSV, this will take a bit...")
for r in np.arange(0,len(df),1):
   print(r)
   poly = Polygon( [(df.LL_lon[r],df.LL_lat[r]), (df.UL_lon[r],df.UL_lat[r]), (df.UR_lon[r],df.UR_lat[r]), (df.LR_lon[r],df.LR_lat[r]), (df.LL_lon[r],df.LL_lat[r])] )
   newdata.loc[index, 'geometry'] = poly
   newdata.loc[index, 'GWRPM25'] = df.GWRPM25[r]
   index = index + 1
print("... done")
# add/create some projection info. --this is WGS84
newdata.crs = from_epsg(4326)

# see what we have; now have a shapefile
print(newdata)
#newdata.to_file('/data/acker/shapefiles/NAv5_mean2015-2017.shp')


#!/home/harkey/miniconda3/envs/workenv/bin/python
# ***********************************************
# this script will read in a one-timestep shapefile of 
# grid polygons+data (output from clip_WHIPSout.py + write_WHIPS2shp.py, or clip_write_WUSTLPM25_2shp.py) 
# and output values of that gridded data
# for each of the polygons in the "newpolys" file (~ census tracts or counties, whatever) 
# area weighting is applied to account for multiple grids overlapping the new polygon

import os 
import gc
import glob
import geopandas as gpd     # version 0.9.0
import xarray as xr         # version 0.20.1; make sure dask is installed too (v. 2021.10.0)
import numpy as np          # version 1.21.2
import pandas as pd         # version 1.4.1
#########

outfilename = "/data/acker/ALA/NA/NAv5_GWRPM25_mean2015-2017_49011.shp" # make sure this has a .shp extension

# shapefiles for the grid and new thing
# if using output from write_WHIPS2shp.py make sure it was run with include_data = y
#gridpolys = gpd.read_file("/home/harkey/for_Lizzy/GLv5_meanGWRPM25_20_22-California.shp")
gridpolys = newdata
newpolys = gpd.read_file('/data/acker/shapefiles/cb_2020_us_county_500k.shp')
# List of state abbreviations for CONUS, Alaska (AK), and Hawaii (HI)
states_to_include = [
    'UT'
]

# Filter counties to only include rows where STUSPS is in the specified list
newpolys = newpolys[newpolys['STUSPS'].isin(states_to_include)]
newpolys = newpolys[newpolys['COUNTYFP'] == '011']

# **********************************************
# **********************************************

# make sure everything's in same coordinate reference system--a projection kind (meters), not a geographic kind (lat/lon)
gridpolys = gridpolys.to_crs("EPSG:5070")
newpolys = newpolys.to_crs("EPSG:5070")
#print(gridpolys)
#print(newpolys)
#print(stophere)

# each row in the new polys shapefile should have a unique GEOID (or TRACTCE if looking at census tracts). 
# do a check
if (len(newpolys) != newpolys.GEOID.nunique()):
   print(' !!!!! this will not work')
   print(stophere)

# loop over indices is like looping over GEOIDs
# https://www.geeksforgeeks.org/different-ways-to-iterate-over-rows-in-pandas-dataframe/
for i in newpolys.index:
    # Select the current polygon
    get_this = newpolys.loc[[i]]

    # Clip the grid polygons to the current polygon
    new_clip = gpd.clip(gridpolys, get_this.geometry, keep_geom_type=True)

    # Check if there are overlapping grid cells
    if not new_clip.empty:
        try:  # For WUSTL PM2.5 data
            newpolys.at[i, 'maxGWRPM25'] = new_clip['GWRPM25'].max()
            newpolys.at[i, 'meanGWRPM25'] = new_clip['GWRPM25'].mean()
            newpolys.at[i, 'p90GWRPM25'] = new_clip['GWRPM25'].quantile(0.9)

        except:  # For WHIPS NO2 data
            newpolys.at[i, 'maxNO2col'] = new_clip['avgNO2col'].max()
    else:
        # If no overlap, assign NaN
        newpolys.at[i, 'maxGWRPM25'] = np.nan
        newpolys.at[i, 'maxNO2col'] = np.nan

    # Clean up memory
    del new_clip
    gc.collect()

# Write the output to a shapefile
newpolys.to_file(outfilename)
