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
    gfile = os.path.join(gdir, f'V5GL0502.HybridPM25.Global.{AllYrs[Yri]*100+1}-{AllYrs[Yri]*100+12}.nc') #problem maybe bc using 04
    efile = os.path.join(edir, f'V5GL0502.HybridPM25E.Global.{AllYrs[Yri]*100+1}-{AllYrs[Yri]*100+12}.nc') #he did 02
    
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

# Define region of interest for visualization
Region = {
    'LAT': [35, 44],  # Latitude range
    'LON': [-116, -106],  # Longitude range
}

# Select year to plot (2016)
Yri = np.where(AllYrs == 2016)[0][0]

for i in range(1, 4):
    # Create map plot
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Miller()})
    ax.set_extent([Region['LON'][0], Region['LON'][1], Region['LAT'][0], Region['LAT'][1]], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.COASTLINE)
    
    # Apply different filtering criteria based on iteration
    if i == 1:
        print(PM25.shape)
        tdata = PM25[:, :, Yri]  # Unfiltered data
        print(tdata.shape)
        title = 'Unfiltered PM2.5'
    elif i == 2:
        tdata = PM25[:, :, Yri]
        tdata[gPM25E[:, :, Yri] / gPM25[:, :, Yri] > 0.9] = np.nan  # Apply uncertainty filter
        #tdata = np.nanmean(tdata, axis=2)
        title = 'Filtered PM2.5'
    elif i == 3:
        tdata = PM25[:, :, Yri]
        tdata[gPM25E[:, :, Yri] / gPM25[:, :, Yri] > 0.9] = np.nan  # Apply uncertainty filter
        tdata[np.sum((gPM25E / gPM25) > 0.9, axis=2) / PM25.shape[2] > 0.25] = np.nan  # Apply expanded filter
        #tdata = np.nanmean(tdata, axis=2)
        title = 'Expanded Filtered PM2.5'
    
    # Plot PM2.5 data on the map
    im = ax.pcolormesh(SATLON, SATLAT, tdata, transform=ccrs.PlateCarree(), cmap='jet')
    plt.colorbar(im, orientation='horizontal', label='PM2.5')
    ax.set_title(f'{title}: {AllYrs[Yri]}')
    print('here')
    plt.savefig(f'aaron_2016_{title}.png')
    plt.show()

print('Finished.')
