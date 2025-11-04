import numpy as np
import xarray as xr
import os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Define full temporal range (years from 2015 to 2022)
AllYrs = np.arange(2015, 2023)

# Define source file directories for different datasets
idir = '/data/acker/WashU_V5_NA/'
gdir = '/data/acker/WashU_V5_GL/'
edir = '/data/acker/WashU_GL_unc/'

init = True  # Flag to initialize variables
for Yri in range(len(AllYrs)):
    # Construct file names for the given year
    ncfile = os.path.join(idir, f'V5NA04.02.HybridPM25.xNorthAmerica.{AllYrs[Yri]*1000+1}-{AllYrs[Yri]*1000+364}.nc')
    gfile = os.path.join(gdir, f'V5GL0502.HybridPM25.Global.{AllYrs[Yri]*100+1}-{AllYrs[Yri]*100+12}.nc')
    efile = os.path.join(edir, f'V5GL0502.HybridPM25E.Global.{AllYrs[Yri]*100+1}-{AllYrs[Yri]*100+12}.nc')
    
    if init:
        # Load regional coordinate grid definition
        SATLAT = xr.open_dataset(ncfile)['lat'].values
        SATLON = xr.open_dataset(ncfile)['lon'].values

        # Load global coordinate grid definition
        SATLAT2 = xr.open_dataset(efile)['lat'].values
        SATLON2 = xr.open_dataset(efile)['lon'].values

        # Find alignment indices for mapping global data to regional grid
        tlatspot2 = slice(np.where(np.abs(SATLAT[0] - SATLAT2) < 1e-4)[0][0], np.where(np.abs(SATLAT[-1] - SATLAT2) < 1e-4)[0][0] + 1)
        tlonspot2 = slice(np.where(np.abs(SATLON[0] - SATLON2) < 1e-4)[0][0], np.where(np.abs(SATLON[-1] - SATLON2) < 1e-4)[0][0] + 1)
        
        # Initialize data storage arrays for PM2.5 values
        PM25 = np.zeros((len(SATLAT), len(SATLON), len(AllYrs)), dtype=np.float32)
        gPM25 = np.copy(PM25)  # Global PM2.5 data
        gPM25E = np.copy(PM25)  # Global PM2.5 error data
        init = False  # Set flag to False after initialization
    
    # Read and store regional PM2.5 values
    PM25[:, :, Yri] = xr.open_dataset(ncfile)['GWRPM25'].values
    
    # Read and map global PM2.5 values to regional grid
    tPM25 = xr.open_dataset(gfile)['GWRPM25'].values
    gPM25[:, :, Yri] = tPM25[tlatspot2, tlonspot2]
    
    # Read and map global PM2.5 error values to regional grid
    tPM25 = xr.open_dataset(efile)['GWRPM25SIGMA'].values
    gPM25E[:, :, Yri] = tPM25[tlatspot2, tlonspot2]
    
    print(f'{AllYrs[Yri]} loaded')

# Process and save filtered data for each year
for Yri in range(len(AllYrs)):
    filtered_PM25 = np.copy(PM25[:, :, Yri])
    filtered_PM25[gPM25E[:, :, Yri] / gPM25[:, :, Yri] > 0.9] = np.nan  # Apply uncertainty filter
    filtered_PM25[np.sum((gPM25E / gPM25) > 0.9, axis=2) / PM25.shape[2] > 0.25] = np.nan  # Apply expanded filter

    # Save filtered PM2.5 data for each year to a separate NetCDF file
    filtered_ds = xr.Dataset(
        {
            "GWRPM25": (('lat', 'lon'), filtered_PM25)
        },
        coords={
            "lat": SATLAT,
            "lon": SATLON
        }
    )
    output_filename = f"/data/acker/ALA/NA/netcdfs/filtered_PM25_{AllYrs[Yri]}.nc"
    filtered_ds.to_netcdf(output_filename)
    print(f"Filtered data for {AllYrs[Yri]} saved to {output_filename}")
    
import numpy as np
import xarray as xr
import os

# Define year ranges for multi-year averages
year_ranges = [(2015, 2017), (2016, 2018), (2017, 2019), (2018, 2020), (2019, 2021), (2020, 2022)]

save_path = "/data/acker/ALA/NA/netcdfs"

# Process and save multi-year averages
for period in year_ranges:
    years = [f"filtered_PM25_{year}.nc" for year in range(period[0], period[1] + 1)]
    files = [xr.open_dataset(os.path.join(save_path, f)) for f in years if os.path.exists(os.path.join(save_path, f))]
    
    if len(files) == 3:
        combined_ds = xr.concat(files, dim="time").mean(dim="time")

        # Save final period NetCDF
        period_output_path = os.path.join(save_path, f"NAv5_GWRPM25_mean_{period[0]}-{period[1]}.nc")
        combined_ds.to_netcdf(period_output_path)
        print(f"Averaged data for {period[0]}-{period[1]} saved to {period_output_path}")

print("Finished computing multi-year averages.")


print('Finished.')