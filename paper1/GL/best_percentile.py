import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
from scipy.stats import rankdata

# TIME PERIODS AND SHAPEFILES
time_periods = {
    "2016-2018": "/data/acker/ALA/GL/GLv5_GWRPM25_mean2016-2018.shp",
    "2017-2019": "/data/acker/ALA/GL/GLv5_GWRPM25_mean2017-2019.shp",
    "2018-2020": "/data/acker/ALA/GL/GLv5_GWRPM25_mean2018-2020.shp",
    "2019-2021": "/data/acker/ALA/GL/GLv5_GWRPM25_mean2019-2021.shp",
    "2020-2022": "/data/acker/ALA/GL/GLv5_GWRPM25_mean2020-2022.shp",
    "2021-2023": "/data/acker/ALA/GL/GLv5_GWRPM25_mean2021-2023.shp"
}

percentiles_to_test = np.arange(0.01, 1.00, 0.01)
results_summary = []

# LOAD COUNTIES ONCE
counties = gpd.read_file('/data/acker/shapefiles/cb_2020_us_county_500k.shp')
states_to_include = ['AL', 'AK', 'AZ', 'AR', 'CA', 'CO', 'CT', 'DE', 'FL', 'GA', 'HI', 'ID', 'DC',
                     'IL', 'IN', 'IA', 'KS', 'KY', 'LA', 'ME', 'MD', 'MA', 'MI', 'MN', 'MS',
                     'MO', 'MT', 'NE', 'NV', 'NH', 'NJ', 'NM', 'NY', 'NC', 'ND', 'OH', 'OK',
                     'OR', 'PA', 'RI', 'SC', 'SD', 'TN', 'TX', 'UT', 'VT', 'VA', 'WA', 'WV',
                     'WI', 'WY']
counties_conus = counties[counties['STUSPS'].isin(states_to_include)]
counties_conus = counties_conus.drop(["COUNTYNS", 'NAMELSAD', 'LSAD', 'ALAND', 'AWATER', 'AFFGEOID'], axis=1)

def load_and_join_satellite_to_counties(shapefile_path):
    df = gpd.read_file(shapefile_path)
    counties_reproj = counties_conus.to_crs(df.crs)
    results = []
    chunk_size = 10000
    for i in range(0, len(df), chunk_size):
        chunk = df.iloc[i:i+chunk_size]
        chunk_result = gpd.sjoin(chunk, counties_reproj[['GEOID', 'geometry','STATEFP','COUNTYFP']],
                                 how="left", predicate="within")
        chunk_result.dropna(subset=['GEOID'], inplace=True)
        results.append(chunk_result)
    return pd.concat(results, ignore_index=True), df.crs

for period, shapefile_path in time_periods.items():
    print(f"\nProcessing {period}")
    
    df_counties, crs = load_and_join_satellite_to_counties(shapefile_path)

    # EPA Data
    gdf_epa = pd.read_csv(f'/data/acker/EPA_DV/{period}.csv')
    gdf_epa.rename(columns={period: 'Design Value'}, inplace=True)
    gdf_epa['State FIPS'] = gdf_epa['State FIPS'].apply(lambda x: str(x).zfill(2))
    gdf_epa['County FIPS'] = gdf_epa['County FIPS'].apply(lambda x: str(x).zfill(3))
    gdf_epa['ID'] = gdf_epa['State FIPS'] + "_" + gdf_epa['County FIPS']
    gdf_epa.loc[gdf_epa['Design Value'] == 0, 'Design Value'] = np.nan
    gdf_epa.dropna(subset=['Design Value'], inplace=True)

    best_pearson = {'percentile': None, 'value': -1}
    best_spearman = {'percentile': None, 'value': -1}

    for p in percentiles_to_test:
        county_percentile = df_counties.groupby('GEOID')['GWRPM25'].quantile(p).reset_index()
        county_percentile = county_percentile.rename(columns={'GWRPM25': 'PM25_90th'})
        county_percentile = county_percentile.merge(counties[['GEOID', 'STATEFP', 'COUNTYFP']], on='GEOID', how='left')
        county_percentile['ID'] = county_percentile['STATEFP'] + "_" + county_percentile['COUNTYFP']

        # Drop counties not in EPA
        county_percentile.loc[~county_percentile['ID'].isin(gdf_epa['ID']), 'PM25_90th'] = np.nan
        county_percentile.dropna(inplace=True)

        # Sort and rank
        county_ranked_epa = gdf_epa.sort_values(by='Design Value', ascending=False).reset_index(drop=True)
        county_ranked_epa.reset_index(inplace=True)

        county_ranked = county_percentile.sort_values(by='PM25_90th', ascending=False).reset_index(drop=True)
        county_ranked.reset_index(inplace=True)

        # Merge
        merged_df = county_ranked.merge(county_ranked_epa[['ID', 'Design Value']], on='ID', how='inner')
        merged_df.dropna(subset=['PM25_90th', 'Design Value'], inplace=True)

        if len(merged_df) < 10:
            continue

        # CDV and CDVE arrays
        cdv = merged_df['Design Value'].values
        cdve = merged_df['PM25_90th'].values

        # Pearson
        cdv_mean = np.mean(cdv)
        cdve_mean = np.mean(cdve)
        numerator = np.sum((cdv - cdv_mean) * (cdve - cdve_mean))
        denominator = np.sqrt(np.sum((cdv - cdv_mean)**2) * np.sum((cdve - cdve_mean)**2))
        pearson_r = numerator / denominator

        # Spearman
        cdv_ranks = rankdata(cdv, method='min')
        cdve_ranks = rankdata(cdve, method='min')
        d_i = cdv_ranks - cdve_ranks
        n = len(cdv)
        spearman_numerator = 6 * np.sum(d_i**2)
        spearman_denominator = n * (n**2 - 1)
        spearman_r = 1 - (spearman_numerator / spearman_denominator)

        if pearson_r > best_pearson['value']:
            best_pearson = {'percentile': p, 'value': pearson_r}
        if spearman_r > best_spearman['value']:
            best_spearman = {'percentile': p, 'value': spearman_r}

    print(f"  Best Pearson @ {best_pearson['percentile']:.2f} => r = {best_pearson['value']:.3f}")
    print(f"  Best Spearman @ {best_spearman['percentile']:.2f} => r_s = {best_spearman['value']:.3f}")

    results_summary.append({
        'time_period': period,
        'best_pearson_percentile': best_pearson['percentile'],
        'best_pearson_r': best_pearson['value'],
        'best_spearman_percentile': best_spearman['percentile'],
        'best_spearman_r_s': best_spearman['value']
    })

# Print final summary
summary_df = pd.DataFrame(results_summary)
print("\n=== Summary Table ===")
print(summary_df)

# Best overall percentiles
overall_best_pearson = summary_df.groupby('best_pearson_percentile')['best_pearson_r'].mean()
overall_best_spearman = summary_df.groupby('best_spearman_percentile')['best_spearman_r_s'].mean()
print(f"\nBest Overall Pearson Percentile: {overall_best_pearson}")
print(f"Best Overall Spearman Percentile: {overall_best_spearman}")
