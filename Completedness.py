import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pickle
from pathlib import Path
import pyproj
import matplotlib.colors as mcolors


ds2 = xr.open_dataset('A:\PHD\Scripts\Data_CUBE\Datacubes\MalaspinaGlacierCube_32608.nc',chunks={'mid_date':100, 'x':100, 'y':100}).v.values
ds = xr.open_dataset('Cubes/temp\MalaspinaGlacierCube_32607_ref.nc', chunks={'mid_date':100, 'x':100, 'y':100}).v

#ds = xr.concat([ds, ds2], dim='mid_date', coords='minimal')

boundary_points = pickle.load(open('boundary.p','rb'))
Xs,Ys = np.meshgrid(ds.x,ds.y)
X = np.vstack((Xs.ravel(),Ys.ravel())).T
tform = pyproj.Transformer.from_crs(crs_from=3338, crs_to=32607, always_xy=True)
fx, fy = tform.transform(boundary_points[:,0], boundary_points[:,1])
import pyproj
import matplotlib.path as path
p = path.Path(np.vstack((fx,fy)).T)
glacier_mask = p.contains_points(X)
glacier_mask = glacier_mask.reshape(ds.shape[1],ds.shape[2])

ds.values[:,glacier_mask==False] = np.nan


# Calculate amount of NaNs per year
#nanyear = ds.isnull().groupby('mid_date.year').sum(dim='mid_date').sum(dim='x').sum(dim='y')
nanyear = ds.notnull().groupby('mid_date.year').sum(dim='mid_date').sum(dim='x').sum(dim='y')


# Calculate amount of NaNs per month
#nanmonth = ds.isnull().groupby('mid_date.month').sum(dim='mid_date').sum(dim='x').sum(dim='y')
nanmonth = ds.notnull().groupby('mid_date.month').sum(dim='mid_date').sum(dim='x').sum(dim='y')

# count amount of slices per year
ds['slices_year'] = ds['mid_date'].dt.year
slices_per_year = np.array([np.sum(ds['slices_year']==i) for i in range(ds.mid_date.dt.year.min().values, 1 + ds.mid_date.dt.year.max().values)])
ratio_nan_year = nanyear/(np.count_nonzero(glacier_mask)*slices_per_year)

# count amount of slices per month
ds['slices_month'] = ds['mid_date'].dt.month
slices_per_month = np.array([np.sum(ds['slices_month']==i).values for i in range(1,13)])
ratio_nan_month = nanmonth/(np.count_nonzero(glacier_mask)*slices_per_month)

# calculate spatial emptiness
ratio_space = ds.count(dim='mid_date')/len(ds.mid_date)

# plot spatial completedness
# Create a custom colormap based on Viridis
cmap = plt.cm.viridis
colors = cmap(np.linspace(0, 1, 256))

# Set the color for 0 to black
colors[0] = [0, 0, 0, 1]

# Create the modified colormap
modified_cmap = mcolors.ListedColormap(colors)

plt.figure(figsize=(10, 10))
plt.imshow(ratio_space.values, cmap=modified_cmap)

cbar = plt.colorbar(aspect=40)
cbar.ax.tick_params(labelsize=14)  # Set the font size for the colorbar label
cbar.set_label('Valid data length/Timeseries length', fontsize=23)
plt.title('Dataset Spatial Completeness', fontsize=20)
plt.xlabel('Easting', fontsize=20)
plt.ylabel('Northing', fontsize=20)

plt.xticks(fontsize=17)
plt.yticks(fontsize=17)

plt.savefig('plots/spatial_complete.png', dpi=300)

plt.figure(figsize=(10,10))
# Create the bar plot
plt.bar(ratio_nan_year.year, ratio_nan_year.values)

# Set the plot title and labels
plt.title('Ratio of Valid glacier pixels per Year', fontsize = 25)
plt.xlabel('Years', fontsize = 20)
plt.ylabel('Nb Valid glacier pixels / Nb glacier pixels', fontsize = 20)


plt.xticks(fontsize=17)
plt.yticks(fontsize=17)


plt.savefig('plots/ratio_temporal_completedness_year.png', dpi=300)  #





plt.figure(figsize=(10,10))
# Create the bar plot
plt.bar(ratio_nan_month.month, ratio_nan_month.values)

# Set the plot title and labels
plt.title('Ratio of Valid glacier pixels per Month', fontsize = 25)
plt.xlabel('Months', fontsize = 20)
plt.ylabel('Nb Valid glacier pixels / Nb glacier pixels', fontsize = 20)
# Generate tick positions and labels for the x-axis
tick_positions = range(1, 13)  # Tick positions from 1 to 12
tick_labels = ['{:02d}'.format(i) for i in tick_positions]  # Format as '01' to '12'


plt.xticks(fontsize=17)
plt.yticks(fontsize=17)


plt.savefig('plots/ratio_temporal_completedness_month.png', dpi=300)  #




plt.figure(figsize=(10,10))
# Create the bar plot
plt.bar(nanmonth.month, nanmonth.values)

# Set the plot title and labels
plt.title('Amount of Valid glacier pixels per Month', fontsize = 25)
plt.xlabel('Months', fontsize = 20)
plt.ylabel('Nb Valid glacier pixels', fontsize = 20)


plt.xticks(fontsize=17)
plt.yticks(fontsize=17)



plt.savefig('plots/temporal_completedness_month.png', dpi=300)  #




plt.figure(figsize=(10,10))
# Create the bar plot
plt.bar(nanyear.year, nanyear.values)

# Set the plot title and labels
plt.title('Amount of Valid glacier pixels per Year', fontsize=25)
plt.ylabel('Nb Valid glacier pixels', fontsize=20)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
# Show the plot
plt.savefig('plots/temporal_completedness_year.png', dpi=300)  #

