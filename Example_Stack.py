import xarray as xr
import dask
from dask.diagnostics import ProgressBar
import pandas as pd

# Define boundaries to trim down the datacubes spatially, and avoid downloading their entire spatial footprint
xmin = -3339099.0323194754
xmax = -3251520.3157346323
ymax = 355078.5266675556
ymin = 273609.2968890464


cubes = ['http://its-live-data.s3.amazonaws.com/datacubes/v02/N50W130/ITS_LIVE_vel_EPSG3413_G0120_X-3350000_Y250000.zarr',
 'http://its-live-data.s3.amazonaws.com/datacubes/v02/N50W140/ITS_LIVE_vel_EPSG3413_G0120_X-3350000_Y350000.zarr',
 'http://its-live-data.s3.amazonaws.com/datacubes/v02/N60W140/ITS_LIVE_vel_EPSG3413_G0120_X-3250000_Y350000.zarr',
 'http://its-live-data.s3.amazonaws.com/datacubes/v02/N60W130/ITS_LIVE_vel_EPSG3413_G0120_X-3250000_Y250000.zarr']


# For example, here I want only the velocity components
variables_to_keep = ['vx','vy', 'v', 'date_dt', 'v','acquisition_date_img1', 'acquisition_date_img2']

# Update the list of variables to keep (otherwise it drops the dimensions)
variables_to_keep += ['mid_date', 'x', 'y']

# List of variables to drop for the download (we drop everything but the variables written below)
variables_drop = [ele for ele in list(
        xr.open_dataset(cubes[0], engine='zarr').variables
        ) if ele not in variables_to_keep
]

# Dummy operation replicating a function going through each 2D footprint pixel of the input datacube

# Create a sample 3D xarray dataset with multiple variables
data_var1 = np.random.rand(10, 100, 100)  # Shape: (mid_date, y, x)
data_var2 = np.random.rand(10, 100, 100)  # Shape: (mid_date, y, x)

coords = {"mid_date": np.arange(10), "y": np.arange(100), "x": np.arange(100)}
ds = xr.Dataset(
    {"var1": (["mid_date", "y", "x"], data_var1), "var2": (["mid_date", "y", "x"], data_var2)},
    coords=coords,
)

def custom_function(pixel):
    # Calculate the monthly average for each "x" and "y" pixel
    monthly_average = pixel.resample(mid_date="1M").mean(dim="mid_date")

    # Create a new time dimension representing the monthly means
    new_mid_date = pd.date_range(start=monthly_average.mid_date.values[0], periods=monthly_average.shape[0], freq="1M")

    # Update the coordinates of the resulting DataArray
    monthly_average = monthly_average.assign_coords(mid_date=new_mid_date)

    return monthly_average



##### METHOD 1 #####
# 1: open datasets
# 2: concatenate
# 3: save as zarr
# 4: operation

# Define path to save the zarrs to
pathsave = 'Datacubes/'


# Load the first datacube so we can append on it
ds = xr.open_zarr(cubes[0],
            chunks=({'mid_date': 100, 'y': 100, 'x': 100}),
            drop_variables=variables_drop
            ).sel(x=slice(xmin, xmax),
                    y=slice(ymax, ymin))

# Download the rest of the cubes and append them to the first one
for n in range(1, len(cubes)):

        # Get the cube's URL
        url = cubes[n]
        
        # Load datacube according to prerequisites (space and variables)
        ds_temp = xr.open_zarr(url,
                    chunks=({'mid_date': 100, 'y': 100, 'x': 100}),
                    drop_variables=variables_drop
                    ).sel(x=slice(xmin, xmax),
                          y=slice(ymax, ymin))
    

        
        # Concatenate the datacubes along the time dimension
        ds = xr.concat((ds, ds_temp), dim = 'mid_date')

# Write the files as zarrs
write_job = ds.to_zarr(f"{pathsave}Datacube.zarr", mode='w', compute=False)
with ProgressBar():
        print(f"Writing to {pathsave}")
        write_job = write_job.compute()

# Open the zarr with dask
ds = xr.open_zarr(f"{pathsave}Datacube.zarr", chunks=({'mid_date': 100, 'y': 100, 'x': 100}))

# Apply the custom function to each 2D pixel of each variable
result = {}
for var_name, var_data in ds.data_vars.items():
    result_var = var_data.map_blocks(custom_function, template=var_data)
    result[var_name] = (["mid_date", "y", "x"], result_var)
    

# Create a new xarray dataset with the results
result_ds = xr.Dataset(result, coords={"mid_date": result["acquisition_date_img1"][1].coords["mid_date"].values, "y": ds["y"].values, "x": ds["x"].values})

# Overwrite the original data
write_job = result_ds.to_zarr(f"{pathsave}Datacube.zarr", mode='w', compute=False)
with ProgressBar():
        print(f"Writing to {pathsave}")
        write_job = write_job.compute()



#### METHOD 2 ####
# 1: open datasets
# 2: operation
# 3: concatenate
# 4: zarr



# Define path to save the zarrs to
pathsave = 'Datacubes/'


# Load the first datacube so we can append on it
ds = xr.open_zarr(cubes[0],
            chunks=({'mid_date': 100, 'y': 100, 'x': 100}),
            drop_variables=variables_drop
            ).sel(x=slice(xmin, xmax),
                    y=slice(ymax, ymin))


# Apply the custom function to each 2D pixel of each variable
result = {}
for var_name, var_data in ds.data_vars.items():
    result_var = var_data.map_blocks(custom_function, template=var_data)
    result[var_name] = (["mid_date", "y", "x"], result_var)
    

# Overwrite ds with the new, smaller dataset
ds = xr.Dataset(result, coords={"mid_date": result["acquisition_date_img1"][1].coords["mid_date"].values, "y": ds["y"].values, "x": ds["x"].values})

# Download the rest of the cubes and append them to the first one
for n in range(1, len(cubes)):

        # Get the cube's URL
        url = cubes[n]
        
        # Load datacube according to prerequisites (space and variables)
        ds_temp = xr.open_zarr(url,
                    chunks=({'mid_date': 100, 'y': 100, 'x': 100}),
                    drop_variables=variables_drop
                    ).sel(x=slice(xmin, xmax),
                          y=slice(ymax, ymin))
        
        ### Reduce the size of the dataset
        # Apply the custom function to each 2D pixel of each variable
        # Apply the custom function to each 2D pixel of each variable
        result = {}
        for var_name, var_data in ds_temp.data_vars.items():
            result_var = var_data.map_blocks(custom_function, template=var_data)
            result[var_name] = (["mid_date", "y", "x"], result_var)
            

        # Create a new xarray dataset with the results
        ds_temp = xr.Dataset(result, coords={"mid_date": result["acquisition_date_img1"][1].coords["mid_date"].values, "y": ds_temp["y"].values, "x": ds_temp["x"].values})

        
        # Concatenate the datacubes along the time dimension
        ds = xr.concat((ds, ds_temp), dim = 'mid_date')



# Write the files as zarrs
write_job = ds.to_zarr(f"{pathsave}{n}.nc", mode='w', compute=False)
with ProgressBar():
        print(f"Writing to {pathsave}")
        write_job = write_job.compute()