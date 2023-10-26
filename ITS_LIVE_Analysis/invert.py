# AUTOGENERATED! DO NOT EDIT! File to edit: ../nbs/03_inversion.ipynb.

# %% auto 0
__all__ = ['urls', 'get_extents', 'design_matrices', 'Inv_reg']

# %% ../nbs/03_inversion.ipynb 4
import numpy as np
import pyproj
import matplotlib.path as path
import s3fs
import zarr
import matplotlib.pyplot as plt
import scipy
from datetime import timedelta
from tqdm import tqdm
import xarray as xr
import re
import pandas as pd
import geopandas as gpd
import matplotlib.path as mplp
import ipyleaflet as ipyl
from ipyleaflet import WMSLayer
import ipywidgets as ipyw
import json
import pandas as pd
from ipyleaflet import Map, WMSLayer, basemaps
from ipywidgets import HTML
from owslib.wms import WebMapService

# %% ../nbs/03_inversion.ipynb 6
urls = []

# %% ../nbs/03_inversion.ipynb 13
# Gather the points on the boundary of the glacier (should later be implemented as a GUI loading the points of a glacier's periphery)
#boundary_points = pickle.load(open('boundary.p','rb'))

# Modify the urls so they can be opened by zarr (replace 'http' by 's3' and delete '.s3.amazonaws.com')
#urls = [url.replace('http','s3') for url in urls]
try: 
    urls = [re.sub(r'http', 's3', url) for url in urls]
    urls = [re.sub(r'\.s3\.amazonaws\.com', '', url) for url in urls]
    # Create storing arrays for the coordinates on-glacier
    X_valid = []
    Y_valid = []
    X_tot = []
    Y_tot = []
    
    # Create an empty directoryimport pickle to hold many variables all tied to the datacubes
    data_dict = {}
    
    # We iterate through the different datacubes so they can each have one instance of the variables below
    for url in urls:
        zarr_store = None # To store the datacube's information and access its variables
        dates = None # To store the dates at which the inversion will give values
        A_m = None # 1st part of the design matrix
        reg_mat_Inv = None # Regularization in time, 2nd part of the design matrix
        mission = None # If you want to invert specifically for one mission in particular ('S1','L8','L9', etc...)
        index_sort = None # Indices representing the sorted dates (from older to most recent)
        inds_mission = None # Indices representing the sorted dates per mission chosen
        ind_tot = None # Indices representing the indices of the pixels on the GOI
        valid_idx = None # Easting and Northing values of the indices above
        proj_cube = None # Projection of the datacube
        
        # Create a dictionary entry for the URL with the desired subsets
        data_dict[url] = {
            'zarr_store': zarr_store,
            'dates': dates,
            'A_m': A_m,
            'reg_mat_Inv': reg_mat_Inv,
            'mission': mission,
            'index_sort': index_sort,
            'inds_mission': inds_mission,
            'dates': dates,
            'ind_tot': ind_tot,
            'valid_idx': valid_idx,
            'proj_cube': proj_cube
        }
except: 
    pass

# %% ../nbs/03_inversion.ipynb 15
def get_extents(url, X_tot, Y_tot, X_valid, Y_valid, data_dict, mission, lamb, derivative, day_interval):

    # Open the zarr files
    fs = s3fs.S3FileSystem(anon=True)
    store = zarr.open(s3fs.S3Map(url, s3=fs))
   
    # Update the dictionnary
    data_dict[url]['zarr_store'] = store

    # Get the cube's projection
    proj_cube = store.attrs['projection']

    # Load X and Y of the dataset
    X = store['x'][:]
    Y = store['y'][:]

    # Store the arrays in the total list
    X_tot.append(X)
    Y_tot.append(Y)

    # Load dimensions
    shape_arr = store['v'].shape
    
    Xs, Ys = np.meshgrid(X, Y)
    points = np.array((Xs.flatten(), Ys.flatten())).T

    idx_valid = []
    
    for b in range(len(gdf_list)):
        mpath = mplp.Path(list(gdf_list[b]['geometry'].to_crs(proj_cube).boundary.explode(index_parts = True).iloc[0].coords))
        glacier_mask = mpath.contains_points(points).reshape(Xs.shape)
        # Grab the indices of the points inside the glacier
        idx_valid.append(np.array(np.where(glacier_mask==True)))
        
    idx_valid = np.hstack(idx_valid)
    # Store the valid indices
    data_dict[url]['valid_idx'] = idx_valid
    
    # Store the cube projection
    data_dict[url]['proj_cube'] = proj_cube
    
    # Store the coordinates of the valid Xs and Ys
    X_valid.append([Xs[idx_valid[0][i], idx_valid[1][i]] for i in range(len(idx_valid[0]))])
    Y_valid.append([Ys[idx_valid[0][i], idx_valid[1][i]] for i in range(len(idx_valid[0]))])
    
    return X_tot, Y_tot, X_valid, Y_valid

# %% ../nbs/03_inversion.ipynb 17
def design_matrices(url, mission, lamb, derivative, day_interval):

    # If you passed 'mission' as an argument, it grabs the appropriate values
    if mission:
        # Get the indices of the mission
        filt1 = np.where(data_dict[urls[url]]['zarr_store']['satellite_img1'][:] == mission)
        filt2 = np.where(data_dict[urls[url]]['zarr_store']['satellite_img2'][:] == mission)
        inds_mission = np.intersect1d(filt1[0],filt2[0])

        # Grab only the indices corresponding to the missions
        mid_dates = np.datetime64('1970-01-01') + np.array(data_dict[urls[url]]['zarr_store']['mid_date'][:], dtype='timedelta64[D]')[inds_mission]
        im1 = np.datetime64('1970-01-01') + np.array(data_dict[urls[url]]['zarr_store']['acquisition_date_img1'][:], dtype='timedelta64[D]')[inds_mission]
        im2 = np.datetime64('1970-01-01') + np.array(data_dict[urls[url]]['zarr_store']['acquisition_date_img2'][:], dtype='timedelta64[D]')[inds_mission]
    else:
        # If 'None' was passed as a mission argument, we grab all the available data.
        inds_mission = None
        mid_dates = np.datetime64('1970-01-01') + np.array(data_dict[urls[url]]['zarr_store']['mid_date'][:], dtype='timedelta64[D]')
        im1 = np.datetime64('1970-01-01') + np.array(data_dict[urls[url]]['zarr_store']['acquisition_date_img1'][:], dtype='timedelta64[D]')
        im2 = np.datetime64('1970-01-01') + np.array(data_dict[urls[url]]['zarr_store']['acquisition_date_img2'][:], dtype='timedelta64[D]')
    
    # Get some arrays
    index_sort = np.argsort(np.datetime64('1970-01-01') + np.array(data_dict[urls[url]]['zarr_store']['mid_date'][:], dtype='timedelta64[D]'))
    mid_dates = mid_dates[index_sort]
    im1 = im1[index_sort]
    im2 = im2[index_sort]


    # Check which im is the smallest (first image, it changes depending on ITS_LIVE's version)
    if im2[0] < im1[0]:
        temp = im1
        im1 = im2
        im2 = temp


    # Create the date array with the new interval dates
    dates_nonum = np.arange(mid_dates[0], mid_dates[-1], timedelta(days=day_interval)).astype(np.datetime64)

    # Convert to numerical
    dates = (dates_nonum - np.datetime64('1970-01-01T00:00:00Z'))/np.timedelta64(1, 's')
    dt_start = (im1 - np.datetime64('1970-01-01T00:00:00Z'))/np.timedelta64(1, 's')
    dt_end = (im2 - np.datetime64('1970-01-01T00:00:00Z'))/np.timedelta64(1, 's')

    # --------------- DESIGN MATRICES --------------- 

    # Initialize matrix
    A_m = np.zeros((mid_dates.shape[0],dates.shape[0]))

    # We have to iterate through the satellite pairs that actually gave a measurement
    for j in range(1, len(mid_dates)):
    # current contents of your for loop

        # Find the middate that is the closest to dt_start (supequal)
        start = np.argmin(np.abs(dates-dt_start[j]))

        # Find the middate that is closest to dt_end (infequal)
        end = np.argmin(dt_end[j] - dates[dates <= dt_end[j]])

        # Divide 1 by the amount of middates between d_start and d_end 
        if end == A_m.shape[1]-1: # If the mid_date is at the end of the array (acquisition im2 equals last mid_date)
            A_m[j, start:] = 1/(1+A_m.shape[1]-start)
        else: # If the measurement is in A's bounds temporally (we can have a satellite pair with the 2nd pair being outside of our mid_dates)
            A_m[j, start:end+1] = 1/(1+end-start) # Attribute to each pixel in the timescale of the satellite pair, the 1/amount of pixel in the pairmid_dates.shape


    # Initialize regularization matrix
    if derivative == 1:
        reg_mat_Inv = np.zeros((A_m.shape[1] -1, A_m.shape[1]))

        for j in range(A_m.shape[1] -1):
            reg_mat_Inv[j, j] = -lamb/day_interval
            reg_mat_Inv[j, j+1] = lamb/day_interval

    elif derivative == 2:
        # Initialize 2nd derivative regularization matrix
        reg_mat_Inv = np.zeros((A_m.shape[1] -1, A_m.shape[1]))

        for j in range(A_m.shape[1] -2):
            reg_mat_Inv[j, j] = lamb/(day_interval**2)
            reg_mat_Inv[j, j+1] = -2*lamb/(day_interval**2)
            reg_mat_Inv[j, j+2] = lamb/(day_interval**2)
            
    data_dict[urls[url]]['A_m'] = A_m
    data_dict[urls[url]]['reg_mat_Inv'] = reg_mat_Inv
    data_dict[urls[url]]['mission'] = mission
    data_dict[urls[url]]['index_sort'] = index_sort
    data_dict[urls[url]]['inds_mission'] = inds_mission
    data_dict[urls[url]]['dates'] = dates_nonum
            
    return data_dict

# %% ../nbs/03_inversion.ipynb 18
def Inv_reg(vObs, data, fillvalue):
 
    # Grab observed velocities
    vObs = vObs[data['index_sort']]
    
    # Filter out the missions we don't want
    if mission:
        vObs = vObs[data['inds_mission']]  
    
    # Mask the NaNs so we don't compute the inversion for empty rows
    mask = np.logical_not(np.equal(vObs, fillvalue))
    
    # Create a masked velocity vector
    vObs_masked = vObs[mask]
    
    # Append regularization terms to dObs
    vObs_masked= np.hstack((vObs_masked, np.zeros((data['reg_mat_Inv'].shape[0]))))
     
    # Assemble the design matrix
    A_des = np.vstack((data['A_m'][mask], data['reg_mat_Inv']))
    
    # Invert the velocities
    vInv = np.linalg.solve(A_des.T@A_des,A_des.T@vObs_masked)
    
    return vInv