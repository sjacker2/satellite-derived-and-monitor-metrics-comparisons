import xarray as xr

import xarray as xr

def process_nc_file(sat_file, unc_file, g_file, remove_pixels=None):
    """
    Processes a single year of data from NetCDF files without converting to DataFrame.
    - sat_file: Satellite PM2.5 dataset
    - unc_file: Uncertainty dataset
    - g_file: Global PM2.5 dataset
    - remove_pixels: List to store pixels with high uncertainty (optional)
    Returns filtered xarray dataset.
    """
    # Load datasets directly as xarray
    sat_data = xr.open_dataset(sat_file)
    unc_data = xr.open_dataset(unc_file).sel(lat=slice(14, 68), lon=slice(-170., -52))
    gl_data = xr.open_dataset(g_file).sel(lat=slice(14, 68), lon=slice(-170., -52))

    # Merge datasets on lat/lon coordinates
    merged = xr.merge([gl_data, unc_data])

    # **Fix 1: Prevent division by zero / NaN issues**
    valid_mask = (merged["GWRPM25"].notnull()) & (merged["GWRPM25SIGMA"].notnull()) & (merged["GWRPM25"] > 0)
    merged = merged.where(valid_mask, drop=True)  # Keep only valid pixels

    # **Calculate percent uncertainty (only for valid data)**
    merged["percent"] = (merged["GWRPM25SIGMA"] / merged["GWRPM25"]) * 100

    # **Fix 2: Properly filter only where uncertainty is high**
    remove_mask = (merged["percent"] >= 90) & valid_mask  # Ensure we only remove valid pixels

    # Extract lat/lon values correctly from the 2D mask
    remove_lat, remove_lon = xr.broadcast(merged["lat"], merged["lon"])
    remove_lat = remove_lat.where(remove_mask, drop=True).values.flatten()
    remove_lon = remove_lon.where(remove_mask, drop=True).values.flatten()

    print(f"Pixels flagged for removal: {len(remove_lat)}")  # Debugging print

    # Store removed pixels (optional tracking)
    if remove_pixels is not None:
        lat_lon_pairs = list(zip(remove_lat, remove_lon))
        remove_pixels.extend(lat_lon_pairs)

    # **Fix 3: Ensure mask is correctly applied to satellite data**
    remove_mask = remove_mask.fillna(False).astype(bool)  # Ensure boolean type
    sat_filtered = sat_data.where(~remove_mask, drop=True)

    return sat_filtered


def process_time_period(time_period):
    """
    Processes all years in a given time period.
    Returns the averaged dataset across years and a list of removed pixels.
    """
    label, *file_groups = time_period
    remove_pixels = []
    yearly_data = []

    # Process each year in the time period
    for year_files in file_groups:
        sat_file, unc_file, g_file = year_files
        filtered = process_nc_file(sat_file, unc_file, g_file, remove_pixels=remove_pixels)
        yearly_data.append(filtered)

    # Compute multi-year mean (over time) without DataFrame conversion
    combined_ds = xr.concat(yearly_data, dim="time").mean(dim="time")

    return label, combined_ds, remove_pixels


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
# Process all time periods
all_remove_pixels = []
all_outputs = []

for time_period in time_periods:
    label, averaged_ds, remove_pixels = process_time_period(time_period)
    print(f"Processed time period: {label} - Shape: {averaged_ds.sizes}")
    all_outputs.append((label, averaged_ds))
    all_remove_pixels.extend(remove_pixels)

# **Fix: Convert removed pixel list to dataset for proper filtering**
remove_counts = xr.Dataset(
    {
        "count": ("pixel", [all_remove_pixels.count(pix) for pix in all_remove_pixels])
    },
    coords={
        "lat": ("pixel", [pix[0] for pix in all_remove_pixels]),
        "lon": ("pixel", [pix[1] for pix in all_remove_pixels])
    }
)

# Set threshold for frequent removals
threshold = len(time_periods) * 0.25

# Find frequently removed pixels
frequently_removed = remove_counts.where(remove_counts["count"] > threshold, drop=True)

# **Fix: Convert frequently_removed to boolean mask**
frequently_removed_mask = xr.zeros_like(averaged_ds["GWRPM25"], dtype=bool)

for i in range(len(frequently_removed["lat"])):
    lat = frequently_removed["lat"].values[i]
    lon = frequently_removed["lon"].values[i]
    frequently_removed_mask = frequently_removed_mask | ((averaged_ds.lat == lat) & (averaged_ds.lon == lon))

# Final filtering for each time period
final_outputs = []
for label, averaged_ds in all_outputs:
    pre_filter_size = averaged_ds.sizes

    # **Fix: Ensure correct mask application**
    filtered_final = averaged_ds.where(~frequently_removed_mask, drop=True)

    post_filter_size = filtered_final.sizes
    print(f"{label}: Before filtering: {pre_filter_size}, After filtering: {post_filter_size}")

    # Drop NaN values explicitly
    filtered_final = filtered_final.dropna(dim="time", how="all")

    final_outputs.append((label, filtered_final))

# **Save results as NetCDF files**
for label, filtered_final in final_outputs:
    output_path = f'/data/acker/ALA/NA/NAv5_GWRPM25_mean_{label}_exact_GL.nc'
    filtered_final.to_netcdf(output_path)
    print(f"NetCDF file saved: {output_path}")

print("Processing complete.")
