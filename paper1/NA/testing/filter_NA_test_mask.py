import xarray as xr
import numpy as np
import json
import os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib import colors

def process_nc_file(sat_file, unc_file, g_file, year, save_path):
    """
    Processes a single year of data from NetCDF files and plots the remove_mask.
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
    num_high_uncertainty_pixels = remove_mask.sum().item()
    print(f"Year {year}: Number of pixels with uncertainty >= 90%: {num_high_uncertainty_pixels}")

    # Store removed pixels
    remove_lat, remove_lon = xr.broadcast(merged["lat"], merged["lon"])
    remove_lat = remove_lat.where(remove_mask, drop=True).values.flatten()
    remove_lon = remove_lon.where(remove_mask, drop=True).values.flatten()
    removed_pixels = list(zip(remove_lat, remove_lon))
    
    # Ensure remove_mask is strictly boolean
    remove_mask = remove_mask.fillna(False).astype(bool).copy()

    # Convert boolean mask to float for interpolation
    remove_mask_float = remove_mask.astype(float)

    remove_mask_interp = remove_mask_float.interp(lat=sat_data.lat, lon=sat_data.lon, method="nearest")

    # Ensure mask is strictly boolean
    remove_mask = remove_mask_interp.astype(bool)
    
    print(f"Shape of remove_mask: {remove_mask.shape}")

    # **Plot the remove_mask on a map**
    plot_remove_mask(remove_mask, year)

    # Apply the mask and introduce NaNs for removed pixels
    sat_filtered = sat_data.copy()
    
    # Debug: Count how many NaNs were introduced
    print(f"NaN count before filtering: {sat_filtered['GWRPM25'].isnull().sum().item()}")
    
    sat_filtered["GWRPM25"] = sat_filtered["GWRPM25"].where(~remove_mask, np.nan)
    
    # Debug: Count how many NaNs were introduced
    print(f"NaN count after filtering: {sat_filtered['GWRPM25'].isnull().sum().item()}")

    # Save filtered dataset
    output_nc = os.path.join(save_path, f"Filtered_GWRPM25_{year}.nc")
    sat_filtered.to_netcdf(output_nc)
    print(f"Filtered NetCDF saved: {output_nc}")

    # Save removed pixels as JSON
    output_json = os.path.join(save_path, f"Removed_Pixels_{year}.json")
    with open(output_json, 'w') as f:
        json.dump(removed_pixels, f)
    print(f"Removed pixels JSON saved: {output_json}")

    return sat_filtered


def plot_remove_mask(remove_mask, year):
    """
    Plots the boolean remove_mask on a map to visualize removed pixels.
    """
    fig = plt.figure(figsize=(12, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())

    # Define extent (North America)
    ax.set_extent([-115.8, -108.5, 38, 44.1], crs=ccrs.PlateCarree())

    # Add map features
    ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor="black")
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5, edgecolor="black")
    ax.add_feature(cfeature.STATES, linewidth=0.2, edgecolor="black")

    # Create 2D grid of lat/lon
    lon2d, lat2d = np.meshgrid(remove_mask.lon, remove_mask.lat)

    # Plot the mask using binary colors
    cmap = colors.ListedColormap(["blue", "red"])  # Blue = kept, Red = removed
    bounds = [0, 0.5, 1]  # Binary mask boundaries
    norm = colors.BoundaryNorm(bounds, cmap.N)

    con = ax.pcolormesh(lon2d, lat2d, remove_mask.values, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())

    # Add colorbar
    cbar = plt.colorbar(con, ax=ax, orientation="vertical", shrink=0.8, pad=0.05)
    cbar.set_label("Mask Values", fontsize=12)
    cbar.set_ticks([0, 1])
    cbar.set_ticklabels(["Kept Pixels", "Removed Pixels"])

    # Set title
    plt.title(f"Removed Pixels for Year {year}", fontsize=14)

    # Show the plot
    plt.savefig(f'mask_{year}.jpg')
    plt.show()


# Example years to process
years = ['2016']#["2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022"]
save_path = "/data/acker/ALA/NA/netcdfs"

# File patterns
for year in years:
    sat_file = f"/data/acker/WashU_V5_NA/V5NA04.02.HybridPM25.xNorthAmerica.{year}001-{year}364.nc"
    unc_file = f"/data/acker/WashU_GL_unc/V5GL04.HybridPM25E.Global.{year}01-{year}12.nc"
    g_file = f"/data/acker/WashU_V5_GL/V5GL04.HybridPM25.Global.{year}01-{year}12.nc"

    process_nc_file(sat_file, unc_file, g_file, year, save_path)
