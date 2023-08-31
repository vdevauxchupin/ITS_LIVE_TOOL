# Fetch the path for essential python scripts

import os
import sys
path =  os.getcwd()
# Avoid Windows users to have issues with how paths are written
path = path.replace('\\','/')

# Import python scripts from notebooks folder
sys.path.append(path + '/scripts')



import netCDF4 as ncf
import torch
import pickle
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
# Import all the necessary packages
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import xarray as xr
import numpy as np
import pandas as pd
import glob
import warnings
import shutil
from pathlib import Path
import pickle


qual_threshold = 0.01
var = 'v'

# Number of modes for the SVD
n_modes = 10
# Number of iterations for the steepest descent
n_iter = 500

TYPE = 'D'
amount = 5


ds = xr.open_dataset("Cubes/temp/1984_2022.nc").v
xs = ds.x.values
ys = ds.y.values
dst = ds.time


# Calculate the quality of each slice
inds_quality = np.array([np.count_nonzero(~np.isnan(ds.values[i]))/(ds.values.shape[1]*ds.values.shape[2]) for i in range(ds.values.shape[0])])

# Keep only good quality images
ds = ds.values[inds_quality>qual_threshold]
dst = dst.values[inds_quality>qual_threshold]

# FROM NOW ON ds IS NOT THE DATASET ANYMORE, BUT THE DATASET'S VARIABLE

# Sort the indices by date
inds = np.argsort(dst)

# Sort dates and variable
dst = dst[inds]
ds = ds[inds]

# Convert the NaNs to -1
ds[np.isnan(ds)==True] = -1
# Create a timedelta
dt = dst[1:] - dst[:-1]
# Convert to float
dt = dt / np.timedelta64(1, 's')
# Grab dimensions of the variable
n_t,n_row,n_col = ds.shape


# Load the boundary points
boundary_points = pickle.load(open('boundary.p','rb'))
Xs,Ys = np.meshgrid(xs,ys)
X = np.vstack((Xs.ravel(),Ys.ravel())).T

import pyproj
import matplotlib.path as path
tform = pyproj.Transformer.from_crs(crs_from=3338, crs_to=3413, always_xy=True)
fx, fy = tform.transform(boundary_points[:,0], boundary_points[:,1])


p = path.Path(np.vstack((fx,fy)).T)
glacier_mask = p.contains_points(X)

plt.scatter(*X[glacier_mask].T)
# Create a glacier mask
glacier_mask = np.hstack([glacier_mask]*n_t).reshape(n_t,n_row,n_col)
mask_f = glacier_mask.reshape(ds.shape[0],-1).T





off_ice_means = np.zeros(ds.shape[0])
for i in range(ds.shape[0]):
    
    # If there are no valid off-ice pixels, off_ice_means will be NaN. We don't want that
    # So we take the average of the previous averages
    off_ice_avg = np.median(ds[i][(ds[i]!=-1) & (glacier_mask[i]==0)])
    if np.isnan(off_ice_avg):
         off_ice_means[i] = np.mean(off_ice_means[:i])
    else:
        off_ice_means[i] = off_ice_avg
    ds[i][ds[i]!=-1] -= off_ice_means[i]
    
Uf = (ds).reshape(ds.shape[0],-1).T


# How many non-orthogonal modes to compute
l = 30

# Regularization strength on column-space (spreads out mode coefficients)
lamda_u = 10.0

# Regularization strength on row-space (makes it more likely to use more modes to describe a velocity field)
lamda_v = 1.0

# Penalizes mode gradients
lamda_x = 1000.0

# Penalizes variation in mode strength through time
lamda_t = 10000.0

# avoid divide by zero in time step (which is used in time-gradient regularization)
dt_0 = torch.from_numpy(dt + 1.0)

# Do some reshaping and rescaling
R = torch.from_numpy(Uf)
Rhat = R.ravel()
I_0 = Rhat!=-1.
I_1 = torch.from_numpy(mask_f.ravel())
I = I_0 & I_1
Rbar = Rhat[I]

Rmean = Rbar.mean()
Rstd = Rbar.std()

Rbar = (Rbar - Rmean)/Rstd

# Define low-rank matrix factors 
U_ = torch.randn(R.shape[0],l,requires_grad=True)
V_ = torch.randn(l,R.shape[1],requires_grad=True)

# Initialize
U_.data[:] *= 1e-2
V_.data[:] *= 1e-2

# Define an optimizer for gradient descent
optimizer = torch.optim.Adam([U_,V_],1e-1)

# Do 500 iterations of gradient descent
for i in range(n_iter):
    optimizer.zero_grad()
    
    # Compute the predicted matrix
    G = U_ @ V_

    # Compute the gradients of the right factor (the mode coefficients through time)
    dudt = (V_[:,1:] - V_[:,:-1])/dt_0
    
    # Reshape the left factor (columns hold the modes) into a grid and takes its gradients.
    U_grid = U_.T.reshape(l,n_row,n_col)    
    dudrow = U_grid[:,1:] - U_grid[:,:-1]
    dudcol = U_grid[:,:,1:] - U_grid[:,:,:-1]
    
    # Flatten and mask the predictions
    Gbar = G.ravel()[I]

    # Compute a variety of negative log likelihoods
    
    # Data misfit
    E_misfit = ((Rbar - Gbar)**2).mean() 
    #E_misfit = (torch.abs(Rbar - Gbar)).mean() 
    
    # norm regularization
    E_reg = lamda_u*(U_**2).sum()/len(Rbar) + lamda_v*(V_**2).sum()/len(Rbar)
    
    # Spatial gradient regularization
    E_space = lamda_x/len(Rbar)*((dudrow**2).sum() + (dudcol**2).sum())  
    
    # Time gradient regularization
    E_time = lamda_t/len(Rbar)*(dudt**2).sum()
    
    # Sum to form total cost
    E = E_misfit + E_reg + E_space + E_time
    
    # Backpropagate gradients
    E.backward()
    
    # Update factors
    optimizer.step()
    print(i,E_misfit.item(),E_reg.item(),E_space.item(),E_time.item())

u,s,v = torch.svd_lowrank(U_.detach()@V_.detach(),n_modes)
U_recon = ((u * s) @ v.T).T.reshape(n_t,n_row,n_col)

U_recon[glacier_mask==False] = np.nan




path =  os.getcwd()
# Avoid Windows users to have issues with how paths are written
path = path.replace('\\','/')

U_recon[glacier_mask==False] = np.nan

datacube_recon = xr.Dataset(
    data_vars=dict(
        v=(["time", "y", "x"], U_recon*365)),
    coords=dict(
        time=(["time"],dst),
        y=(["y"], ys),
        x=(["x"], xs),
    ),
)



# Save .nc file
datacube_recon.to_netcdf(f'Cubes/final/Interpolated_1984_2022.nc')
