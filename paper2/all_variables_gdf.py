import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
from shapely.geometry import Point

df = gpd.read_file('/data/acker/ALA/paper2/explain_variables_1.shp')

df = df.rename(columns={'monitor__1':'monitor_coverage_pct', 'classifica':'classification', 'Design Val':'Design Value', 'monitor_co':'monitor_count'})

fires = gpd.read_file('/data/acker/ALA/paper2/fire_counties.shp')

fires = fires.rename(columns={'fire_regio':'fire_region'})

mountains = gpd.read_file('/data/acker/ALA/paper2/mountain_counties.shp')

mountains = mountains.rename(columns={'mountain_r':'mountain_region'})

deserts = gpd.read_file('/data/acker/ALA/paper2/desert_counties.shp')

deserts = deserts.rename(columns={'desert_reg':'desert_region'})

urban = gpd.read_file('/data/acker/ALA/paper2/urban_counties.shp')

urban = urban.rename(columns={'urban_cate':'urban_category'})

final_df = df.merge(fires[['GEOID', 'fire_region']], on='GEOID', how='left')

final_df = final_df.merge(mountains[['GEOID', 'mountain_region']], on='GEOID', how='left')

final_df = final_df.merge(deserts[['GEOID', 'desert_region']], on='GEOID', how='left')

final_df = final_df.merge(urban[['GEOID', 'urban_category']], on='GEOID', how='left')

print(final_df.columns)

final_df.to_file('/data/acker/ALA/paper2/all_variables.shp')


