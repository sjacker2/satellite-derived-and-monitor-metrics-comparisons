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
filtered_20 = process_nc_file('/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2015001-2015364.nc')
filtered_21 = process_nc_file('/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2016001-2016364.nc')
filtered_22 = process_nc_file('/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.2017001-2017364.nc')

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
minlat = 10
minlon = -179.9
maxlat = 72
maxlon = -50
#10, 72), lon=slice(-180., -50
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
newdata.to_file('/data/acker/shapefiles/NAv5_mean2015-2017.shp')


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

outfilename = "/data/acker/ALA/NA/NAv5_GWRPM25_mean2015-2017_counties.shp" # make sure this has a .shp extension

# shapefiles for the grid and new thing
# if using output from write_WHIPS2shp.py make sure it was run with include_data = y
#gridpolys = gpd.read_file("/home/harkey/for_Lizzy/GLv5_meanGWRPM25_20_22-California.shp")
gridpolys = newdata
newpolys = gpd.read_file('/data/acker/shapefiles/cb_2020_us_county_500k.shp')
# List of state abbreviations for CONUS, Alaska (AK), and Hawaii (HI)
states_to_include = [
    'AL', 'AK', 'AZ', 'AR', 'CA', 'CO', 'CT', 'DE', 'FL', 'GA', 'HI', 'ID', 'DC',
    'IL', 'IN', 'IA', 'KS', 'KY', 'LA', 'ME', 'MD', 'MA', 'MI', 'MN', 'MS',
    'MO', 'MT', 'NE', 'NV', 'NH', 'NJ', 'NM', 'NY', 'NC', 'ND', 'OH', 'OK',
    'OR', 'PA', 'RI', 'SC', 'SD', 'TN', 'TX', 'UT', 'VT', 'VA', 'WA', 'WV',
    'WI', 'WY'
]

# Filter counties to only include rows where STUSPS is in the specified list
newpolys = newpolys[newpolys['STUSPS'].isin(states_to_include)]

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
   # next line of code from https://stackoverflow.com/questions/17071871/how-do-i-select-rows-from-a-dataframe-based-on-column-values
   get_this = newpolys.loc[ (newpolys['GEOID'] == newpolys['GEOID'][i]) ]
   # now use that row (slice) of newpolys as the thing to clip:
   new_clip = gpd.clip(newpolys, get_this.geometry, keep_geom_type=True) # keep geom. otherwise

   # now do a spatial join and get the weight (%) of each grid polygon
   # this section adapted from https://stackoverflow.com/questions/70051743/calculating-spatial-averages-for-each-country-after-spatial-join
   # more about lambda functions here: https://www.geeksforgeeks.org/applying-lambda-functions-to-pandas-dataframe
   grid_new_join = (                    # whole series of stacked functions/operations here to fill "grid_new_join"
       gpd.sjoin(new_clip, gridpolys)   # spatial join, default is left outer join, so will keep only things that intersect with new_clip
                                          # this gets rid of the geometry info. of the grid data--so next,
       .merge(                             # merge the grid geometry to the join output--and rename the grid geom. column to keep track of what's what
           gridpolys.loc[:, "geometry"],  
           left_on="index_right",
           right_index=True,
           suffixes=("", "_grid") 
       )
       .assign(                            # assign = add new column to dataframe (which is the output of the merge)
           overlap=lambda d: (             # lambda function, supposedly it's the fast way to do it instead of a formal/separate definition 
           d["geometry"]                   # intersect the two geometries. "d" is what's being passed to the lambda function (?)
           .intersection(gpd.GeoSeries(d["geometry_grid"], crs="EPSG:5070"))
           .area                           # get the fractional area
            / (d["geometry"].area)         # divide by total new polygon area to get a percent
           )#.round(6)                      # round to 6 decimal places
       )
   )
   # just to reiterate: "overlap" is the percent of a new polygon overlapped by a grid box. So the total overlap should be ~ 100%.
   #print("after lambda function",grid_new_join.columns)
   # remove extra columns just because
   try: # this is for newpolys = counties 
      grid_new_join = grid_new_join.drop(['COUNTYNS','AFFGEOID', 'NAMELSAD', 'STUSPS', 'LSAD', 'index_right', 'ALAND','AWATER'], axis=1)
   except: # this is for newpolys = census tracts
      grid_new_join = grid_new_join.drop(['NAMELSAD', 'MTFCC', 'FUNCSTAT', 'ALAND', 'AWATER'], axis=1)
#   print(grid_new_join)

   # clean things up to save memory
   del(new_clip)
   del(get_this)
   gc.collect()

   print("there are ",len(grid_new_join)," grid cells overlapping the ( #",i,"/",len(newpolys)," ) polygon with GEOID = ",newpolys['GEOID'][i])
#   print("**** the total overlap percent = ",sum(grid_new_join['overlap'].values))
   # https://stackoverflow.com/questions/56409042/is-there-a-way-to-add-a-column-to-a-geopandas-dataframe-using-a-single-value-geo
   # apply weights (overlap), sum, and assign value to new column   
   try: # for WUSTL PM2.5 output from clip_write_WUSTLPM25_2shp.py 
      newpolys.at[i,'awGWRPM25'] = sum(grid_new_join['GWRPM25'].values * grid_new_join['overlap'].values)
      newpolys.at[i,'maxGWRPM25'] = max(grid_new_join['GWRPM25'].values)
      newpolys.at[i, 'p90GWRPM25'] = grid_new_join['GWRPM25'].quantile(0.9)
      newpolys.at[i, 'meanGWRPM25'] = grid_new_join['GWRPM25'].mean()
   except: # for WHIPS NO2 output from clip_WHIPSout.py and write_WHIPS2shp.py
      newpolys.at[i,'awNO2col'] = sum(grid_new_join['avgNO2col'].values * grid_new_join['overlap'].values)
#   print(newpolys)
 
   del(grid_new_join)
   gc.collect()

# write the output!
newpolys.to_file(outfilename)

