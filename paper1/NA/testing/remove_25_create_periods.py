import xarray as xr
import numpy as np
import json
import glob
import os

def load_removed_pixels(json_files):
    """
    Loads removed pixels from all JSON files and returns a frequency count.
    """
    pixel_counts = {}

    for file in json_files:
        with open(file, 'r') as f:
            removed_pixels = json.load(f)
            for pixel in removed_pixels:
                pixel_tuple = tuple(pixel)
                if pixel_tuple in pixel_counts:
                    pixel_counts[pixel_tuple] += 1
                else:
                    pixel_counts[pixel_tuple] = 1

    return pixel_counts

def apply_final_mask(dataset, frequent_pixels):
    """
    Removes the most frequently removed pixels from the dataset.
    """
    
    # Replace frequent removal pixels with NaNs
    print(f"NaN count before filtering: {dataset['GWRPM25'].isnull().sum().item()}")
    
    # Convert lat/lon pairs to NumPy arrays
    frequent_lats, frequent_lons = np.array(frequent_pixels).T

    # Apply faster masking using `isin()`
    mask_lat = np.isin(dataset["lat"], frequent_lats)
    mask_lon = np.isin(dataset["lon"], frequent_lons)

    # Create a final mask using broadcasting
    mask = mask_lat[:, None] & mask_lon[None, :]
    dataset["GWRPM25"] = dataset["GWRPM25"].where(~mask, np.nan)

    print(f"NaN count after filtering: {dataset['GWRPM25'].isnull().sum().item()}")

    return dataset

def process_final_files(nc_files, json_files, save_path, threshold_factor=0.25):
    """
    Processes all years by removing the most frequently removed pixels.
    Then averages years together for each period.
    """
    # Load removed pixels
    pixel_counts = load_removed_pixels(json_files)

    # Determine threshold
    threshold = len(json_files) * threshold_factor
    frequent_pixels = [pixel for pixel, count in pixel_counts.items() if count > threshold]

    # Process each dataset
    processed_years = []
    for file in nc_files:
        dataset = xr.open_dataset(file)
        print(f"Processing {file}...")
        filtered_dataset = apply_final_mask(dataset, frequent_pixels)

        # Save filtered dataset per year
        output_path = file.replace(".nc", "_final.nc")
        filtered_dataset.to_netcdf(output_path)
        print(f"Filtered dataset saved: {output_path}")

        processed_years.append((file, filtered_dataset))

    # Now average years into periods
    periods = {
        "2015-2017": ["Filtered_GWRPM25_2015_final.nc", "Filtered_GWRPM25_2016_final.nc", "Filtered_GWRPM25_2017_final.nc"],
        "2016-2018": ["Filtered_GWRPM25_2016_final.nc", "Filtered_GWRPM25_2017_final.nc", "Filtered_GWRPM25_2018_final.nc"],
        "2017-2019": ["Filtered_GWRPM25_2017_final.nc", "Filtered_GWRPM25_2018_final.nc", "Filtered_GWRPM25_2019_final.nc"],
        "2018-2020": ["Filtered_GWRPM25_2018_final.nc", "Filtered_GWRPM25_2019_final.nc", "Filtered_GWRPM25_2020_final.nc"],
        "2019-2021": ["Filtered_GWRPM25_2019_final.nc", "Filtered_GWRPM25_2020_final.nc", "Filtered_GWRPM25_2021_final.nc"],
        "2020-2022": ["Filtered_GWRPM25_2020_final.nc", "Filtered_GWRPM25_2021_final.nc", "Filtered_GWRPM25_2022_final.nc"]
    }

    for period, years in periods.items():
        files = [xr.open_dataset(os.path.join(save_path, f)) for f in years if os.path.exists(os.path.join(save_path, f))]
        if len(files) == 3:
            combined_ds = xr.concat(files, dim="time").mean(dim="time")

            # Save final period NetCDF
            period_output_path = os.path.join(save_path, f"NAv5_GWRPM25_mean_{period}.nc")
            combined_ds.to_netcdf(period_output_path)
            print(f"Averaged period dataset saved: {period_output_path}")

# Load files
save_path = "/data/acker/ALA/NA/netcdfs"
nc_files = glob.glob(os.path.join(save_path, "Filtered_GWRPM25_*.nc"))
json_files = glob.glob(os.path.join(save_path, "Removed_Pixels_*.json"))

# Run final processing
process_final_files(nc_files, json_files, save_path)
