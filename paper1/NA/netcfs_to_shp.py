import xarray as xr
import geopandas as gpd
import pandas as pd
import glob
from shapely.geometry import Point

def netcdf_to_shapefile(nc_file, save_path):
    """
    Converts a filtered NetCDF file to a shapefile using .to_dataframe().
    
    - nc_file: Path to the filtered NetCDF file
    - save_path: Directory to save the shapefile
    """
    # Load dataset
    ds = xr.open_dataset(nc_file)

    # Convert to DataFrame
    df = ds["GWRPM25"].to_dataframe().reset_index()  # Converts lat/lon grid into tabular format
    print('here')
    # Drop NaN values
    df = df.dropna(subset=["GWRPM25"])
    print(df)

    # Create geometry column using lat/lon
    df["geometry"] = df.apply(lambda row: Point(row["lon"], row["lat"]), axis=1)
    print('here')
    # Convert to GeoDataFrame
    gdf = gpd.GeoDataFrame(df, geometry="geometry", crs="EPSG:4326")  # WGS84 CRS

    # Define output shapefile path
    shapefile_name = nc_file.split("/")[-1].replace(".nc", ".shp")
    shapefile_path = f"{save_path}/{shapefile_name}"

    # Save as shapefile
    gdf.to_file(shapefile_path)
    print(f"Shapefile saved: {shapefile_path}")

def process_all_netcdfs(nc_dir, shapefile_output_dir):
    """
    Processes all filtered NetCDF files in a directory and converts them to shapefiles.
    
    - nc_dir: Directory containing filtered NetCDF files
    - shapefile_output_dir: Directory to save shapefiles
    """
    nc_files = glob.glob(f"{nc_dir}/NAv5_GWRPM25_mean*.nc")

    for nc_file in nc_files:
        try:
            netcdf_to_shapefile(nc_file, shapefile_output_dir)
        except Exception as e:
            print(f"Error processing {nc_file}: {e}")

# Example usage
nc_directory = "/data/acker/ALA/NA/netcdfs"  # Directory where filtered NetCDF files are stored
shapefile_output_directory = "/data/acker/ALA/NA/final"  # Output directory for shapefiles

process_all_netcdfs(nc_directory, shapefile_output_directory)
