import xarray as xr
import numpy as np
import json
import glob

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
    Removes the most frequently removed pixels from the dataset. Need to change so that does adding nans like in first script.
    """
    mask = xr.zeros_like(dataset["GWRPM25"], dtype=bool)
    
    for lat, lon in frequent_pixels:
        mask |= (dataset.lat == lat) & (dataset.lon == lon)
    
    return dataset.where(~mask, drop=True)

def process_final_files(nc_files, json_files, save_path, threshold_factor=0.25):
    """
    Processes all time periods by removing the most frequently removed pixels.
    """
    # Load removed pixels
    pixel_counts = load_removed_pixels(json_files)

    # Determine threshold
    threshold = len(json_files) * threshold_factor
    frequent_pixels = [pixel for pixel, count in pixel_counts.items() if count > threshold]

    # Process each dataset
    for file in nc_files:
        dataset = xr.open_dataset(file)
        print(dataset)
        filtered_dataset = apply_final_mask(dataset, frequent_pixels)
        print(filtered_dataset)
        # Save filtered dataset
        output_path = file.replace(".nc", "_filtered.nc")
        filtered_dataset.to_netcdf(output_path)
        print(f"Filtered dataset saved: {output_path}")

# Load files
nc_files = glob.glob("/data/acker/ALA/NA/netcdfs/NAv5_GWRPM25_mean_*.nc")
json_files = glob.glob("/data/acker/ALA/NA/netcfs/removed_pixels_*.json")

# Run final processing
process_final_files(nc_files, json_files, "/data/acker/ALA/NA/netcdfs")
