#%% Import
import os
import sys
path =  os.getcwd()
# Avoid Windows users to have issues with how paths are written
path = path.replace('\\','/')

# Import python scripts from notebooks folder
sys.path.append(path + '/scripts')


import matplotlib.pyplot as plt
import pyproj
from datetime import timedelta, date
import pickle
import netCDF4 as ncf
import torch
os.environ['KMP_DUPLICATE_LIB_OK']='True'
# Import all the necessary packages
import matplotlib.path as mpath
import xarray as xr
import numpy as np
import pandas as pd
import glob
import warnings
import shutil
from pathlib import Path
import time



t0 = time.time()
#%% Import Dataset and crop it to TOI (Time Of Interest)
# Load Dataset

ds = xr.open_dataset('A:/PHD/Scripts/Clean_Datacube/Cubes/temp/MalaspinaGlacierCube_32607.nc')

# Attribute the variables to numpy arrays
vel = ds.v.values
mid_date = ds.mid_date.values
dt_start = ds.acquisition_img2.values
dt_end = ds.acquisition_img1.values
xs = ds.x.values
ys = ds.y.values

# Close xarray dataset and free memory
ds.close()
del ds

# Sort the arrays 
indx = np.argsort(mid_date)
vel = vel[indx]
vel = vel[14000:]

# Round all time arrays to days, because it doesn't make sense to have nanoseconds
mid_date = np.array(mid_date[indx], dtype = 'datetime64[D]')[14000:]
dt_start = np.array(dt_start[indx], dtype = 'datetime64[D]')[14000:]
dt_end = np.array(dt_end[indx], dtype = 'datetime64[D]')[14000:]

t1 =  time.time()-t0


#%% 
# SVD Interpolation
qual_threshold = 0.01
var = 'v'

# Number of modes for the SVD
n_modes = 10
# Number of iterations for the steepest descent
n_iter = 500


ds = vel.copy()
# Convert the NaNs to -1
ds[np.isnan(ds)==True] = -1
# Get time
dst = mid_date
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
tform = pyproj.Transformer.from_crs(crs_from=3338, crs_to=32607, always_xy=True)
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

t0 = time.time()

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

#U_recon[glacier_mask==False] = np.nan



t2 = time.time()-t0 


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
datacube_recon.to_netcdf(f'Cubes/final/SVD.nc')


#%% 
# Inversion

input_array = vel.copy()

# Decide on interpolation interval
day_interval = 5
# Create the date array with the new interval dates
dates_nonum = np.arange(mid_date[0], mid_date[-1], timedelta(days=day_interval)).astype(np.datetime64)

# Convert to numerical
dates = (dates_nonum - np.datetime64('1970-01-01T00:00:00Z'))/np.timedelta64(1, 's')
dt_start = (dt_start - np.datetime64('1970-01-01T00:00:00Z'))/np.timedelta64(1, 's')
dt_end = (dt_end - np.datetime64('1970-01-01T00:00:00Z'))/np.timedelta64(1, 's')


# Set all the points outside of the boundary of the glacier to NaN
boundary_points = pickle.load(open('boundary.p','rb'))
Xs,Ys = np.meshgrid(xs,ys)
X = np.vstack((Xs.ravel(),Ys.ravel())).T
tform = pyproj.Transformer.from_crs(crs_from=3338, crs_to=32607, always_xy=True)
fx, fy = tform.transform(boundary_points[:,0], boundary_points[:,1])
import pyproj
import matplotlib.path as path
p = path.Path(np.vstack((fx,fy)).T)
glacier_mask = p.contains_points(X)
glacier_mask = glacier_mask.reshape(input_array.shape[1],input_array.shape[2])
input_array[:,glacier_mask==False] = np.nan

# Grab the indices of the points inside the glacier
idx_valid = np.array(np.where(glacier_mask==True))


# create the output array
output_array = np.full((len(dates), input_array.shape[1], input_array.shape[2]), np.nan)


lamb = 100
import time
from scipy.optimize import nnls


# create the output array
output_array = np.full((len(dates), input_array.shape[1], input_array.shape[2]), np.nan)

t0 = time.time()

# Where the parallelization happens, we iterate through the different pixels of input_array
for i in range(idx_valid.shape[1]):


        # --------------- DESIGN MATRIX --------------- 



        # Select the corresponding pixel in the array
        pt = input_array[:,idx_valid[0][i],idx_valid[1][i]]

        # get indices of non_nan values in the timeseries
        non_nan_idx = np.array(np.where(~np.isnan(pt))[0])

        # Initialize matrix
        A = np.zeros((non_nan_idx.shape[0],dates.shape[0]))

        # We have to iterate through the satellite pairs that actually gave a measurement
        from tqdm import tqdm

        for j in tqdm(range(len(non_nan_idx))):
    # current contents of your for loop

            # Find the middate that is the closest to dt_start (supequal)
            start = np.argmin(np.abs(dates-dt_start[non_nan_idx[j]]))
            
            # Find the middate that is closest to dt_end (infequal)
            end = np.argmin(dt_end[non_nan_idx[j]] - dates[dates <= dt_end[non_nan_idx[j]]])

            # Divide 1 by the amount of middates between d_start and d_end 
            if end == A.shape[1]-2: # If the mid_date is at the end of the array (acquisition im2 equals last mid_date)
                A[j, start:] = 1/(1+A.shape[1]-start)
            else: # If the measurement is in A's bounds temporally (we can have a satellite pair with the 2nd pair being outside of our mid_dates)
                A[j, start:end+1] = 1/(1+end-start) # Attribute to each pixel in the timescale of the satellite pair, the 1/amount of pixel in the pair



        # --------------- INVERSION --------------- #



        M=np.shape(A)[1]
        # This is standard inversion (I took the formulas out of Parker's book)

        # Simpson rule weights (these are needed for integration (calculating the norm). Actually this might not even matter
        dg = np.ones(M)*1/3
        dg[1:M-1] += 1/3
        dg[1:M-2:2] += 2/3
        W = np.sqrt(np.diag(dg))

        # this is just a differencing matrix
        delta = np.diag(np.ones(M)) - np.diag(np.ones(M-1),-1)
        delta[0,0] = 0

        # vRR^Tv is the roughness norm
        R = np.matmul(W,delta)

        # This expression comes from differencing the minimizing functional with respect to v_i
        LHS = np.matmul(np.transpose(A),A) + 1/lamb*np.matmul(np.transpose(R),R)

        # put the interpolated (modeled) velocities in the output array
        output_array[:,idx_valid[0][i],idx_valid[1][i]], res = nnls(LHS, np.matmul(np.transpose(A),pt[non_nan_idx]))




t3 = time.time()-t0 

# Some values can get very weird
output_array[output_array>10000] = np.nan


datacube_recon = xr.Dataset(
    data_vars=dict(
        v=(["time", "y", "x"], output_array)),
    coords=dict(
        time=(["time"],dates_nonum),
        y=(["y"], ys),
        x=(["x"], xs),
    ),
)



# Save .nc file
datacube_recon.to_netcdf(f'Cubes/final/Inv.nc')

#%% 
# SVD on Inv

ds = xr.open_dataset('Cubes/final/Inv.nc')

# Attribute the variables to numpy arrays
vel = ds.v.values
mid_date = ds.time.values
xs = ds.x.values
ys = ds.y.values

# Close xarray dataset and free memory
ds.close()
del ds

# Sort the arrays 
indx = np.argsort(mid_date)
vel = vel[indx]

# Round all time arrays to days, because it doesn't make sense to have nanoseconds
mid_date = np.array(mid_date[indx], dtype = 'datetime64[D]')

t1 =  time.time()-t0


# SVD Interpolation
qual_threshold = 0.01
var = 'v'

# Number of modes for the SVD
n_modes = 10
# Number of iterations for the steepest descent
n_iter = 1000


ds = vel.copy()
# Filter out weird values
ds[ds>10000] = np.nan
# Convert the NaNs to -1
ds[np.isnan(ds)==True] = -1
# Get time
dst = mid_date
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
tform = pyproj.Transformer.from_crs(crs_from=3338, crs_to=32607, always_xy=True)
fx, fy = tform.transform(boundary_points[:,0], boundary_points[:,1])


p = path.Path(np.vstack((fx,fy)).T)
glacier_mask = p.contains_points(X)

plt.scatter(*X[glacier_mask].T)
# Create a glacier mask
glacier_mask = np.hstack([glacier_mask]*n_t).reshape(n_t,n_row,n_col)
mask_f = glacier_mask.reshape(ds.shape[0],-1).T

    
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

t0 = time.time()

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

#U_recon[glacier_mask==False] = np.nan



t2 = time.time()-t0 


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
datacube_recon.to_netcdf(f'Cubes/final/SVD_Inv.nc')





# %%
# Compare the datasets
amount = 5
TYPE = 'D'
inv = xr.open_dataset('Cubes/final/Inv.nc').resample(time=f'{amount}{TYPE}').median(dim="time", skipna = True)
svd = xr.open_dataset('Cubes/final/SVD.nc').resample(time=f'{amount}{TYPE}').median(dim="time", skipna = True)
svd_inv = xr.open_dataset('Cubes/final/SVD_Inv.nc').resample(time=f'{amount}{TYPE}').median(dim="time", skipna = True)

#### --- Original Cube
ds = xr.open_dataset('A:/PHD/Scripts/Clean_Datacube/Cubes/temp/MalaspinaGlacierCube_32607.nc')

boundary_points = pickle.load(open('boundary.p','rb'))
Xs,Ys = np.meshgrid(ds.x.values,ds.y.values)
X = np.vstack((Xs.ravel(),Ys.ravel())).T
tform = pyproj.Transformer.from_crs(crs_from=3338, crs_to=32607, always_xy=True)
fx, fy = tform.transform(boundary_points[:,0], boundary_points[:,1])
import pyproj
import matplotlib.path as path
p = path.Path(np.vstack((fx,fy)).T)
glacier_mask_32607 = p.contains_points(X)
glacier_mask_32607 = glacier_mask_32607.reshape(ds.v.shape[1],ds.v.shape[2])



# Attribute the variables to numpy arrays
vel = ds.v.values
mid_date = ds.mid_date.values
dt_start = ds.acquisition_img2.values
dt_end = ds.acquisition_img1.values
xs = ds.x.values
ys = ds.y.values

# Close xarray dataset and free memory
ds.close()
del ds

# Sort the arrays 
indx = np.argsort(mid_date)
vel = vel[indx]
vel = vel[14000:]

# Round all time arrays to days, because it doesn't make sense to have nanoseconds
mid_date = np.array(mid_date[indx], dtype = 'datetime64[D]')[14000:]


##### --- SVD TOTAL
t =  ds = xr.open_dataset('Cubes/final/Interpolated_1984_2022.nc').time.values
t_mask = np.logical_and(t>np.datetime64('2016-01-01'), t<np.datetime64('2020-01-02'))
ds = ds = xr.open_dataset('Cubes/final/Interpolated_1984_2022.nc').sel(time=t_mask).resample(time=f'{amount}{TYPE}').median(dim="time", skipna = True)

boundary_points = pickle.load(open('boundary.p','rb'))
Xs,Ys = np.meshgrid(ds.x.values,ds.y.values)
X = np.vstack((Xs.ravel(),Ys.ravel())).T
tform = pyproj.Transformer.from_crs(crs_from=3338, crs_to=3413, always_xy=True)
fx, fy = tform.transform(boundary_points[:,0], boundary_points[:,1])
import pyproj
import matplotlib.path as path
p = path.Path(np.vstack((fx,fy)).T)
glacier_mask_3413 = p.contains_points(X)
glacier_mask_3413 = glacier_mask_3413.reshape(ds.v.shape[1],ds.v.shape[2])

# Attribute the variables to numpy arrays
vel2 = ds.v.values
mid_date2 = ds.time.values

# Close xarray dataset and free memory
ds.close()
del ds

# Sort the arrays 
indx = np.argsort(mid_date2)
vel2 = vel2[indx]
vel2 = vel2

# Round all time arrays to days, because it doesn't make sense to have nanoseconds
mid_date2 = np.array(mid_date2[indx], dtype = 'datetime64[D]')


plt.figure()
plt.plot(inv.time.values, inv.v.values[:,140, 185], label='inv')
plt.plot(svd.time.values, svd.v.values[:,140, 185], label='svd_limited')
plt.plot(svd_inv.time.values, svd_inv.v.values[:,140, 185], label='svd_inv')
plt.plot(mid_date, vel[:,140,185], label='original')
plt.plot(mid_date2, vel2[:,200,200], label='svd_total')

plt.legend()

# Compute stds
std_inv = np.std(inv.v.values[:,glacier_mask_32607], axis = 1)
std_svd_lim = np.std(svd.v.values[:,glacier_mask_32607==True], axis = 1)
std_svd_inv = np.std(svd_inv.v.values[:,glacier_mask_32607==True], axis = 1)
std_svd_tot = np.std(vel2[:,glacier_mask_3413==True], axis = 1)
std_original = np.nanstd(vel[:,glacier_mask_32607==True], axis = 1)

plt.figure()
plt.plot(inv.time.values, std_inv, label='inv')
plt.plot(svd.time.values, std_svd_lim, label='svd_limited')
plt.plot(svd_inv.time.values, std_svd_inv, label='svd_inv')
plt.plot(mid_date2, std_svd_tot, label='original')
plt.plot(mid_date, std_original, label='svd_total')
plt.legend()

# Mean
m_inv = np.mean(inv.v.values[:,glacier_mask_32607==True], axis = 1)
m_svd_lim = np.mean(svd.v.values[:,glacier_mask_32607==True], axis = 1)
m_svd_inv = np.mean(svd_inv.v.values[:,glacier_mask_32607==True], axis = 1)
m_svd_tot = np.mean(vel2[:,glacier_mask_3413==True], axis = 1)
m_original = np.nanmean(vel[:,glacier_mask_32607==True], axis = 1)

plt.figure()
plt.plot(inv.time.values, m_inv, label='inv')
plt.plot(svd.time.values, m_svd_lim, label='svd_limited')
plt.plot(svd_inv.time.values, m_svd_inv, label='svd_inv')
plt.plot(mid_date2, m_svd_tot, label='original')
plt.plot(mid_date, m_original, label='svd_total')
plt.legend()



# %%
