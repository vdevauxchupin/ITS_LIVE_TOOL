import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pickle
from pathlib import Path
import pyproj
import matplotlib.colors as mcolors

ds2 = xr.open_dataset('A:\PHD\Scripts\Data_CUBE\Datacubes\MalaspinaGlacierCube_32608.nc',chunks={'mid_date':100, 'x':100, 'y':100}).v.values
ds = xr.open_dataset('A:\PHD\Scripts\Clean_Datacube\Cubes/temp\MalaspinaGlacierCube_32607.nc', chunks={'mid_date':100, 'x':100, 'y':100}).v

#ds = xr.concat([ds, ds2], dim='mid_date', coords='minimal')

boundary_points = pickle.load(open('A:\PHD\Scripts\Clean_Datacube/boundary.p','rb'))
Xs,Ys = np.meshgrid(ds.x,ds.y)
X = np.vstack((Xs.ravel(),Ys.ravel())).T
tform = pyproj.Transformer.from_crs(crs_from=3338, crs_to=32607, always_xy=True)
fx, fy = tform.transform(boundary_points[:,0], boundary_points[:,1])
import pyproj
import matplotlib.path as path
p = path.Path(np.vstack((fx,fy)).T)
glacier_mask = p.contains_points(X)
glacier_mask = glacier_mask.reshape(ds.shape[1],ds.shape[2])

d_mask = np.broadcast_to(glacier_mask, (ds.shape[0],ds.shape[1],ds.shape[2]))

ds.values[(d_mask == True) & (np.isnan(ds.values))] = False

# Expand the mask to match the shape of the 3D array


# Calculate amount of values per year
nanyear = (ds==True).groupby('mid_date.year').sum(dim='mid_date').sum(dim='x').sum(dim='y')

# Calculate amount of values per month
nanmonth = (ds==True).groupby('mid_date.month').sum(dim='mid_date').sum(dim='x').sum(dim='y')

# count amount of slices per year
ds['slices_year'] = ds['mid_date'].dt.year
slices_per_year = np.array([np.sum(ds['slices_year']==i) for i in range(ds.mid_date.dt.year.min().values, 1 + ds.mid_date.dt.year.max().values)])
ratio_nan_year = nanyear.values/(np.sum(glacier_mask==True)*slices_per_year.values)

# count amount of slices per month
ds['slices_month'] = ds['mid_date'].dt.month
slices_per_month = np.array([np.sum(ds['slices_month']==i).values for i in range(1,13)])
ratio_values_month = nanmonth/(np.sum(glacier_mask==True)*slices_per_month)

# calculate spatial emptiness
ratio_space = (ds==True).sum(dim='mid_date')/len(ds.mid_date)

# plot spatial completedness
# Create a custom colormap based on Viridis
cmap = plt.cm.viridis
colors = cmap(np.linspace(0, 1, 256))

# Set the color for 0 to black
colors[0] = [0, 0, 0, 1]

# Create the modified colormap
modified_cmap = mcolors.ListedColormap(colors)

plt.figure()
plt.imshow(ratio_space.values, cmap = modified_cmap)
plt.colorbar()


# Calculate completedness per year


# Amount of slices per year
plt.figure(figsize=(15,15))
# Generate the histogram showing amount of slices per year
hist, bins, _ = plt.hist(ds['slices_years'], bins=range(1984, 2019))

# Calculate the bin centers
bin_centers = 0.5 * (bins[:-1] + bins[1:])

# Increase frequency of xticks
plt.xticks(bin_centers, [int(x) for x in bin_centers], rotation=45)  # Convert to integer labels

plt.xlabel('Year')
plt.ylabel('Frequency')
plt.title('Year Histogram')
plt.show()