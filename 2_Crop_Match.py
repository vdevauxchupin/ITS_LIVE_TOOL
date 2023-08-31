import xarray as xr
from rasterio.enums import Resampling

amount = 5
TYPE = 'D'

# Open Dataset and grab only one variable size 18206x334x333
ds1 = xr.open_dataset('Reprojected_Cube.nc')
ds1 = ds1.v.sortby(ds1.time).resample(time=f'{amount}{TYPE}').median(dim="time", skipna = True)
ds2 = xr.open_dataset('Concat_2018_2022.nc')
ds2 = ds2.v.sortby(ds2.mid_date).resample(mid_date=f'{amount}{TYPE}').median(dim="mid_date", skipna = True)
ds2 = ds2.rename(mid_date = 'time')

# Check if ds1 includes ds2 spatially
if (ds2.x.min() >= ds1.x.min()) and (ds2.x.max() <= ds1.x.max()) and (ds2.y.min() >= ds1.y.min()) and (ds2.y.max() <= ds1.y.max()):
    # Select the spatial extent of the second dataset
    x_min, x_max = ds2.x.min().values, ds2.x.max().values
    y_min, y_max = ds2.y.min().values, ds2.y.max().values

    # Crop the first dataset to the spatial extent of the second dataset
    ds1 = ds1.sel(x=slice(x_min, x_max), y=slice(y_max, y_min))
    # Set the resolution of ds1_cropped to match ds2
    ds1 = ds1.interp_like(ds2, method='linear') 

# Check if ds2 includes ds1 spatially
elif (ds1.x.min() >= ds2.x.min()) and (ds1.x.max() <= ds2.x.max()) and (ds1.y.min() >= ds2.y.min()) and (ds1.y.max() <= ds2.y.max()):
    # Select the spatial extent of the second dataset
    x_min, x_max = ds1.x.min().values, ds1.x.max().values
    y_min, y_max = ds1.y.min().values, ds1.y.max().values

    # Crop the first dataset to the spatial extent of the second dataset
    ds2 = ds2.sel(x=slice(x_min, x_max), y=slice(y_max, y_min))
    # Set the resolution of ds1_cropped to match ds2
    ds2 = ds2.interp_like(ds2, method='linear') 

# Otherwise, the datasets do not overlap spatially
else:
    print('The datasets do not overlap spatially')

ds1 = xr.concat((ds1,ds2),dim='time')

ds1.to_netcdf('Cubes/temp/1984_2022.nc')


