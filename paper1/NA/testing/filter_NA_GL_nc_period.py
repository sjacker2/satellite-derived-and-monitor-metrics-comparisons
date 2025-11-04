import xarray as xr
import numpy as np
import json

def process_nc_file(sat_file, unc_file, g_file):
    """
    Processes a single year of data from NetCDF files.
    Returns the filtered dataset and removed pixel locations.
    """
    # Load datasets
    sat_data = xr.open_dataset(sat_file)
    unc_data = xr.open_dataset(unc_file).sel(lat=slice(14, 68), lon=slice(-170., -52))
    gl_data = xr.open_dataset(g_file).sel(lat=slice(14, 68), lon=slice(-170., -52))

    # Merge datasets
    merged = xr.merge([gl_data, unc_data])

    # Ensure valid data
    valid_mask = (merged["GWRPM25"].notnull()) & (merged["GWRPM25SIGMA"].notnull()) & (merged["GWRPM25"] > 0)
    merged = merged.where(valid_mask, drop=True)

    # Compute percent uncertainty
    merged["percent"] = (merged["GWRPM25SIGMA"] / merged["GWRPM25"]) * 100
    

    # Identify pixels with high uncertainty
    remove_mask = (merged["percent"] >= 90) & valid_mask
    # Count the number of pixels where percent uncertainty >= 90
    num_high_uncertainty_pixels = remove_mask.sum().item()  # Convert to integer
    print(f"Number of pixels with uncertainty >= 90%: {num_high_uncertainty_pixels}")
    remove_lat, remove_lon = xr.broadcast(merged["lat"], merged["lon"])
    remove_lat = remove_lat.where(remove_mask, drop=True).values.flatten()
    remove_lon = remove_lon.where(remove_mask, drop=True).values.flatten()

    # Store removed pixels
    removed_pixels = list(zip(remove_lat, remove_lon))
    # Filter satellite data
    remove_mask = remove_mask.fillna(False).astype(bool).copy()
    # Ensure remove_mask has the same shape as sat_data before applying it
   
    # Convert boolean mask to float for interpolation
    remove_mask_float = remove_mask.astype(float)

    # Perform interpolation to align lat/lon coordinates
    #remove_mask_interp = remove_mask_float.interp(lat=sat_data.lat, lon=sat_data.lon, method="nearest")
    #remove_mask = remove_mask.sel(lat=sat_data.lat, lon=sat_data.lon, method="nearest")
    remove_mask_interp = remove_mask_float.reindex(lat=sat_data.lat, lon=sat_data.lon, method="nearest")

    # Debug: Check `remove_mask` before applying
    print(f"Remove mask sample:\n{remove_mask.isel(lat=slice(0, 5), lon=slice(0, 5))}")

    # Ensure mask is strictly boolean
    remove_mask = remove_mask_interp.astype(bool)
    print(f"Remove mask type: {remove_mask.dtype}")


    # Apply mask explicitly
    sat_filtered = sat_data.copy()  # Ensure we don’t modify the original
    # Debug: Count how many NaNs were introduced
    print(f"NaN count before filtering: {sat_filtered['GWRPM25'].isnull().sum().item()}")
    sat_filtered["GWRPM25"] = sat_filtered["GWRPM25"].where(~remove_mask, np.nan)

    # Debug: Count how many NaNs were introduced
    print(f"NaN count after filtering: {sat_filtered['GWRPM25'].isnull().sum().item()}")

    # Drop all NaNs aggressively
    #filtered_data = sat_filtered.where(~sat_filtered["GWRPM25"].isnull(), drop=True)
    # Print final dataset shape
    #print(f"After filtering {sat_file}, dataset shape: {filtered_data['GWRPM25'].shape}")


    return sat_filtered, removed_pixels

def process_time_period(label, file_groups, save_path):
    """
    Processes a single time period and saves the results.
    """
    removed_pixels_all_years = []
    yearly_data = []

    for year_files in file_groups:
        sat_file, unc_file, g_file = year_files
        filtered, removed_pixels = process_nc_file(sat_file, unc_file, g_file)
        yearly_data.append(filtered)
        print(filtered)
        #print(removed_pixels)
        removed_pixels_all_years.extend(removed_pixels)

    # Compute multi-year mean
    combined_ds = xr.concat(yearly_data, dim="time").mean(dim="time")
    # Debug: Count how many NaNs were introduced
    print(f"NaN count after filtering and averaging: {combined_ds['GWRPM25'].isnull().sum().item()}")
    # Save dataset
    output_path = f"{save_path}/NAv5_GWRPM25_mean_{label}_exact_GL.nc"
    combined_ds.to_netcdf(output_path)
    print(f"NetCDF file saved: {output_path}")

    # Save removed pixels as JSON
    remove_pixels_path = f"{save_path}/removed_pixels_{label}.json"
    with open(remove_pixels_path, 'w') as f:
        json.dump(removed_pixels_all_years, f)

    print(f"Removed pixels saved: {remove_pixels_path}")

# Example usage
year1 = '2015'
year2 = '2016'
year3 = '2017'
time_periods = [
    ("2015-2017", [
        (f'/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.{year1}001-{year1}364.nc',
         f'/data/acker/WashU_GL_unc/V5GL04.HybridPM25E.Global.{year1}01-{year1}12.nc',
         f'/data/acker/WashU_V5_GL/V5GL04.HybridPM25.Global.{year1}01-{year1}12.nc'),
        (f'/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.{year2}001-{year2}364.nc',
         f'/data/acker/WashU_GL_unc/V5GL04.HybridPM25E.Global.{year2}01-{year2}12.nc',
         f'/data/acker/WashU_V5_GL/V5GL04.HybridPM25.Global.{year2}01-{year2}12.nc'),
        (f'/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.{year3}001-{year3}364.nc',
         f'/data/acker/WashU_GL_unc/V5GL04.HybridPM25E.Global.{year3}01-{year3}12.nc',
         f'/data/acker/WashU_V5_GL/V5GL04.HybridPM25.Global.{year3}01-{year3}12.nc')
    ])
]

#you need to switch it to do the averaging after you remove the 25%
# Run processing for a single period (modify save path as needed)
process_time_period("2015-2017", time_periods[0][1], "/data/acker/ALA/NA/netcdfs")
