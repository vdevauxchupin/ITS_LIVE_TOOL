import xarray as xr
from rasterio.enums import Resampling
import numpy as np


# Fetch the path for essential python scripts

amount = 5
TYPE = 'D'
# Open Dataset and grab only one variable
ds = xr.open_dataset('Cubes/temp/MalaspinaGlacierCube_32607.nc').v

# Get the unique indices otherwise we can not resample
#_, index = np.unique(ds.mid_date, return_index=True)
#t = ds.mid_date
ds = ds.sortby('mid_date').resample(mid_date=f'{amount}{TYPE}').median(dim='mid_date', skipna=True)
ds = ds.rename(mid_date = 'time')

# Import pyproj to get the CRS
import pyproj

# Define the CRS
crs_from = pyproj.CRS('EPSG:32607')
crs_to = pyproj.CRS('EPSG:3413')

# Attribute a projection to the dataset (currently missing its projection)
ds = ds.rio.write_crs("EPSG:32607", inplace=True)

# Reproject the dataset
ds = ds.rio.reproject("EPSG:3413",resampling=Resampling.bilinear)

# Save it as a netcdf
ds.to_netcdf('Cubes/temp/Reprojected_Cube.nc',compute=True)