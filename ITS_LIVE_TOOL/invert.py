__all__ = ['mission', 'lamb', 'derivative', 'day_interval', 'sdate', 'edate', 'GPU', 'make_input_dict', 'create_data_dict',
           'get_extents', 'design_matrices', 'Inv_reg']

# Import the packaged
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
from ipyleaflet import Map, WMSLayer, basemaps, GeoData
from ipywidgets import HTML
from owslib.wms import WebMapService
import os
try:
    import torch
except:
    print('torch import failed. Is it installed ?')

# Create a dictionnary holding the urls of the datacubes with their extents
def make_input_dict(coords, gpdf, urls):
    
    mod_urls = [re.sub(r'http', 's3', url) for url in urls]
    mod_urls = [re.sub(r'\.s3\.amazonaws\.com', '', url) for url in mod_urls]
    
    d = {'coords': coords,
     'gpdf': gpdf,
     'urls': mod_urls
        }
    return d

# Create a dictionnary holding all the variables for each datacube
def create_data_dict(urls, mission, lamb, derivative, day_interval, data_map):
    
    # Modify the urls so they can be opened by zarr (replace 'http' by 's3' and delete '.s3.amazonaws.com')
    urls = [re.sub(r'http', 's3', url) for url in urls]
    urls = [re.sub(r'\.s3\.amazonaws\.com', '', url) for url in urls]

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
        valid_idx = None # Easting and Northing values of the indices above
        proj_cube = None # Projection of the datacube
        mask_dates = None # Mask that filters out dates outside of desired date range
        reg_mat_space_Inv = None # Space regularization matrix
        reg_mat_time_Inv = None # Time regularization matrix

        # Create a dictionary entry for the URL with the desired subsets
        data_dict[url] = {
            'zarr_store': zarr_store,
            'A_m': A_m,
            'mission': mission,
            'index_sort': index_sort,
            'inds_mission': inds_mission,
            'dates': dates,
            'valid_idx': valid_idx,
            'proj_cube': proj_cube,
            'mask_dates': mask_dates
        }

    # Create storing arrays for the coordinates on-glacier
    X_valid = []
    Y_valid = []
    X_tot = []
    Y_tot = []
            
    return X_tot, Y_tot, X_valid, Y_valid, data_dict, urls



# Gather the extends of each datacube and map the glacier's extent in the datacube's referentials
def get_extents(url, X_tot, Y_tot, X_valid, Y_valid, mission, lamb, derivative, day_interval, data_map, data_dict):

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
    
    for b in range(len(data_map.added_glaciers)):
        mpath = mplp.Path(data_map.added_glaciers[b]['geometry'].to_crs(proj_cube).boundary.explode(index_parts = True).iloc[0].coords[:])
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


# Map how the cubes are organized in space and time
def cubes_intersection(X_tot, Y_tot, X_valid, Y_valid, data_dict, urls, spatial_regularization, mission, lamb, derivative, day_interval, sdate, edate, nb_pts_tot, glacier_centerline):

    # CUBES INTERSECTION IN SPACE
    
    # Create Eastings and Northings arrays based on the Eastings and Northings of the datacubes
    X_arr = np.unique(np.hstack(X_tot))
    Y_arr = np.unique(np.hstack(Y_tot))

    # Crop to the GOI (so we avoid over-filling our matrix with NaNs)
    x_min = np.where(np.min(np.hstack(X_valid)) == X_arr)[0][0]
    x_max = np.where(np.max(np.hstack(X_valid)) == X_arr)[0][0]
    y_min = np.where(np.min(np.hstack(Y_valid)) == Y_arr)[0][0]
    y_max = np.where(np.max(np.hstack(Y_valid)) == Y_arr)[0][0]

    # Define the boundaries of the GOI.
    X_MIN = min(x_min-1, x_max+1)
    X_MAX = max(x_min-1, x_max+1)
    Y_MIN = min(y_min-1, y_max+1)
    Y_MAX = max(y_min-1, y_max+1)

    # If the GOI is not entirely covered by the datacubes, we need to crop the GOI so it fits the datacube. If not, then we expand the GOI to fit the datacube.
    if X_MIN < 0:
        X_MIN = 0
    if X_MAX > len(X_arr):
        X_MAX = len(X_arr)
    if Y_MIN < 0:
        Y_MIN = 0
    if Y_MAX > len(Y_arr):
        Y_MAX = len(Y_arr)


    # And now search the indices corresponding to the coordinates 
    x_matches = np.hstack([[np.where(i == X_arr[X_MIN:X_MAX])[0][0] for i in row] for row in X_valid]).astype(int)
    y_matches = np.hstack([[np.where(i == Y_arr[Y_MIN:Y_MAX])[0][0] for i in row] for row in Y_valid]).astype(int)
    coords_matches = list(zip(y_matches, x_matches))

    # Create an array representing the glacier
    template = np.zeros((len(Y_arr[Y_MIN:Y_MAX]), len( X_arr[X_MIN:X_MAX])))
    template[y_matches, x_matches] = 1

    # Determine to which datacube each match belongs to
    arr_belong = np.zeros((len(y_matches)))
    p = 0
    for i in range(len(urls)):
        arr_belong[p:data_dict[urls[i]]['valid_idx'].shape[1]+p] = i
        p += data_dict[urls[i]]['valid_idx'].shape[1]
    cube_belong = np.full((len(Y_arr[Y_MIN:Y_MAX]), len( X_arr[X_MIN:X_MAX])), -1)
    cube_belong[y_matches, x_matches] = arr_belong


    # Concatenate valid indices for the cubes
    neighbor_idx = []
    for i in range(len(urls)):
        neighbor_idx.append(data_dict[urls[i]]['valid_idx'])
    neighbor_idx = np.hstack(neighbor_idx)

    # Create a 2D array the size of template that links neighbor_idx entries and their position in the 2D template
    position = np.zeros(template.shape, dtype = int)
    for i in range(len(y_matches)):
        position[y_matches[i], x_matches[i]] = i    

    # Create a list that gathers the coordinates of central points and their neighboring points, their respective cubes and position in flattened template
    P = [[]]

    # Define the offsets for the surrounding cells
    offsets = [(-1, 0), (1, 0), (0, -1), (0, 1)]

    for i in range(len(y_matches)):
       # Iterate over each surrounding cell
        sum_surrounding = 0
        pt = []
        pt.append([y_matches[i], x_matches[i], cube_belong[y_matches[i], x_matches[i]], position[y_matches[i], x_matches[i]]])

        if spatial_regularization:
            for offset in offsets:

                surrounding_i = y_matches[i] + offset[0]
                surrounding_j = x_matches[i] + offset[1]

                # Check if the coordinates are within the array bounds
                if (0 <= surrounding_i < template.shape[0] and
                    0 <= surrounding_j < template.shape[1]):
                    # Check if the value of the surrounding cell is 1
                    if template[surrounding_i, surrounding_j] == 1:
                        sum_surrounding += 1
                        pt.append([surrounding_i, surrounding_j, cube_belong[surrounding_i, surrounding_j], position[surrounding_i, surrounding_j]])
        P.append(pt)

    # Get rid of first entry (which is empty)
    P = P[1:]
    
    # Reproject coordinates of centerline in coordinates of datacube
    geometry = glacier_centerline.main_centerline.to_crs(epsg=int(data_dict[urls[0]]['proj_cube'])).iloc[0]['geometry']

    # Extract coordinates
    coords = list(geometry.coords)

    # Explode the coordinates
    x_coords, y_coords = zip(*coords)

    # Create empty lists to store the indices
    x_indices_centerline = []
    y_indices_centerline = []

    # Find the index of the closest value in X for each x coordinate
    for x in x_coords:
        idx = (np.abs(X_arr[X_MIN:X_MAX] - x)).argmin()
        x_indices_centerline.append(idx)

    # Find the index of the closest value in Y for each y coordinate
    for y in y_coords:
        idx = (np.abs(Y_arr[Y_MIN:Y_MAX] - y)).argmin()
        y_indices_centerline.append(idx)
    
    # Plot the array gathering the intersection of the cubes, and the glacier's centerline
    plt.pcolormesh(X_arr[X_MIN:X_MAX], Y_arr[Y_MIN:Y_MAX], template)
    plt.scatter(X_arr[X_MIN:X_MAX][x_indices_centerline], Y_arr[Y_MIN:Y_MAX][y_indices_centerline])

    # CUBES INTERSECTION IN TIME
    min_date = []
    max_date = []

    for url in range(len(urls)):
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

        # If sdate is later than the first available date, we find its corresponding index
        try:
            sdate_ind = np.where(mid_dates >= np.datetime64(sdate))[0][0]
        except:
            sdate_ind = 0

        # If edate is sooner than the last available date, we find its corresponding index
        try:
            edate_ind = np.where(mid_dates > np.datetime64(edate))[0][0]
        except:
            edate_ind = None

        # Create a False/True mask where True if the date is in the desired range
        mask_dates = np.full(mid_dates.shape, False)
        mask_dates[sdate_ind:edate_ind] = True

        # Keep only the values within the desired range
        mid_dates = mid_dates[mask_dates]
        im1 = im1[mask_dates]
        im2 = im2[mask_dates]

        # Check which im is the smallest (first image, it changes depending on ITS_LIVE's version)
        if im2[0] < im1[0]:
            temp = im1
            im1 = im2
            im2 = temp

        min_date.append(np.min(im1))
        max_date.append(np.max(im2))

    # Determine the min and max of the dates available in the datacubes
    min_date = np.min(min_date)
    max_date = np.max(max_date)

    # Create the design matrices for each datacube
    # Distance between two cell centers
    space_interval = np.sqrt((X_arr[1]-X_arr[0])**2 + (Y_arr[1]-Y_arr[0])**2)
    
    return space_interval, min_date, max_date, template, P, neighbor_idx, x_indices_centerline, y_indices_centerline, X_arr, X_MIN, X_MAX, Y_arr, Y_MIN, Y_MAX, x_matches, y_matches


# Create the design matrices for each datacube
def design_matrices(urls, min_date, max_date, mission, lamb, derivative, day_interval, space_interval, sdate, edate, data_dict):

    for url in range(len(urls)):
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
    
        # If sdate is later than the first available date, we find its corresponding index
        try:
            sdate_ind = np.where(mid_dates >= np.datetime64(sdate))[0][0]
        except:
            sdate_ind = 0
        
        # If edate is sooner than the last available date, we find its corresponding index
        try:
            edate_ind = np.where(mid_dates > np.datetime64(edate))[0][0]
        except:
            edate_ind = None
        
        # Create a False/True mask where True if the date is in the desired range
        mask_dates = np.full(mid_dates.shape, False)
        mask_dates[sdate_ind:edate_ind] = True
    
        # Keep only the values within the desired range
        mid_dates = mid_dates[mask_dates]
        im1 = im1[mask_dates]
        im2 = im2[mask_dates]
    
        # Check which im is the smallest (first image, it changes depending on ITS_LIVE's version)
        if im2[0] < im1[0]:
            temp = im1
            im1 = im2
            im2 = temp
    
        # Create the date array with the new interval dates, from the first date available to two time intervals after ther last date available (python bounds thing)
        dates_nonum = np.arange(min_date, max_date+(day_interval*2), timedelta(days=day_interval)).astype(np.datetime64)
    
        # Convert to numerical
        dates = (dates_nonum - np.datetime64('1970-01-01T00:00:00Z'))/np.timedelta64(1, 's')
        dt_start = (im1 - np.datetime64('1970-01-01T00:00:00Z'))/np.timedelta64(1, 's')
        dt_end = (im2 - np.datetime64('1970-01-01T00:00:00Z'))/np.timedelta64(1, 's')
        delta_t = dates[1]-dates[0] # Calculate delta time between two dates
    
        # --------------- DESIGN MATRICES --------------- 
    
        # Initialize matrix
        A_m = np.zeros((mid_dates.shape[0],dates.shape[0]))
    
        # We have to iterate through the satellite pairs that actually gave a measurement
        for j in range(1, len(mid_dates)):
        # current contents of your for loop
            beg_dates = (dates[dates<=dt_start[j]] - dt_start[j]) # Get the dates inferior or equal to the 1st acquisition's date
            beg_diff =  delta_t + beg_dates[-1] # Calculate the difference between the closest date supequal to the 1st acquisition's date
            
            # If the 1st acquisition's date falls exactly on a date, then their difference in time is 0. Otherwise, we calculate it.
            if beg_diff == 0:
                start = len(beg_dates) - 1
            else:
                start = len(beg_dates)
                A_m[j, start-1:start] = beg_diff
                
            
            
            end_ind = np.where(dates == dates[dates>=dt_end[j]][0])[0][0] # Get the closest date supequal to the 2nd acquisition's date
            end_diff = delta_t - (dates[end_ind] - dt_end[j]) # Calculate the difference between this date and the 2nd acquisition's date
            
            # If the 2nd acquisition's date falls exactly on a date, then their difference in time is 0. Otherwise, we calculate it.
            if end_diff == 0:
                end = end_ind
            else:
                end = end_ind - 1
                A_m[j, end:end+1] = end_diff
            
            A_m[j, start:end] = delta_t # For every whole interval, attribute the value of the delta time
            A_m[j] /= ((dt_end[j]-dt_start[j])) # Divide everything by the time difference between acquisition 1 and 2. Results are convex (add up to 1)

        # Store A_m in its dict
        data_dict[urls[url]]['A_m'] = A_m
        data_dict[urls[url]]['mission'] = mission
        data_dict[urls[url]]['index_sort'] = index_sort
        data_dict[urls[url]]['inds_mission'] = inds_mission
        data_dict[urls[url]]['mask_dates']= mask_dates

    # Initialize regularization matrices (time and space)
    reg_mat_space = []
    if derivative == 1:
        reg_mat_time = np.zeros((A_m.shape[1] - derivative, A_m.shape[1]))
        reg_mat_space.append(np.diag(np.diag(np.ones((A_m.shape[1], A_m.shape[1]))*(-lamb/space_interval), 0))[:-derivative]) # Negative
        reg_mat_space.append(np.diag(np.diag(np.ones((A_m.shape[1], A_m.shape[1]))*(lamb/space_interval), 0))[:-derivative]) # Positive

        for j in range(A_m.shape[1] -1):
            reg_mat_time[j, j] = -lamb/day_interval
            reg_mat_time[j, j+1] = lamb/day_interval

    elif derivative == 2:
        # Initialize centered differences regularization matrix
        reg_mat_time = np.zeros((A_m.shape[1], A_m.shape[1]))
        reg_mat_space.append(np.diag(np.diag(np.ones((A_m.shape[1], A_m.shape[1]))*(-2*lamb/(space_interval**2)), 0))) # Negative
        reg_mat_space.append(np.diag(np.diag(np.ones((A_m.shape[1], A_m.shape[1]))*(lamb/(space_interval**2)), 0))) # Positive

        for j in range(A_m.shape[1] -2):
            reg_mat_time[j, j] = lamb/(day_interval**2)
            reg_mat_time[j, j+1] = -2*lamb/(day_interval**2)
            reg_mat_time[j, j+2] = lamb/(day_interval**2)

    return reg_mat_space, reg_mat_time, dates_nonum, dates, data_dict

# Create regularization matrices for the cross patter 
def extend_Matrices(reg_mat_time, nb_pts_tot, reg_mat_space, spatial_regularization, dates, template):

    # Get the shapes of the regularization matrices
    SRT = reg_mat_time.shape
    time_reg_mat = np.zeros((nb_pts_tot*SRT[0], nb_pts_tot*SRT[1]))

    # Stack time_reg_temp on the diagonal of time_reg
    for i in range(0,nb_pts_tot):
        time_reg_mat[i*SRT[0]:(i+1)*SRT[0], i*SRT[1]:(i+1)*SRT[1]] = reg_mat_time

    if spatial_regularization:
        # Initialize spatial regulization matrix
        space_reg_mat = np.zeros((2*SRT[0], time_reg_mat.shape[1]))

        space_reg_mat[:, :SRT[1]] = np.tile(reg_mat_space[0], (2,1))
        space_reg_mat[:SRT[1], SRT[1]:3*SRT[1]] = np.tile(reg_mat_space[0], 2)
        space_reg_mat[SRT[1]:2*SRT[0], 3*SRT[1]:] = np.tile(reg_mat_space[0], 2)
        SRS = space_reg_mat.shape
    else: 
        space_reg_mat = None
        SRS = None

    # Update the size of temporal reg mat
    SRT = time_reg_mat.shape

    # Initialize velocity matrices
    len_pt_inverted = len(dates)
    vxInv = np.full((len_pt_inverted, template.shape[0], template.shape[1]), np.nan)
    vyInv = np.full((vxInv.shape), np.nan)
    
    return vxInv, vyInv, len_pt_inverted, SRT, SRS, space_reg_mat, time_reg_mat

# Function to run the inversion for each point and organize the data in the host arrays
def looper(i, name_cube, vxInv, vyInv, GPU, spatial_regularization, nb_pts_tot, time_reg_mat, space_reg_mat, len_pt_inverted, data_dict, urls, P, neighbor_idx, x_matches, y_matches, X_arr, X_MIN, X_MAX, Y_arr, Y_MIN, Y_MAX, device, dates_nonum, x_indices_centerline, y_indices_centerline, SRT, SRS):

    
    # Get all the points used in the inversion. 1st point is the one we are inverting for
    vxObs = [np.array(data_dict[urls[P[i][v][2]]]['zarr_store']['vx'][:,  neighbor_idx[:,P[i][v][3]][0], neighbor_idx[:,P[i][v][3]][1]][data_dict[urls[P[i][v][2]]]['index_sort']][data_dict[urls[P[i][v][2]]]['mask_dates']], dtype=np.float64) for v in range(len(P[i]))]
    vyObs = [np.array(data_dict[urls[P[i][v][2]]]['zarr_store']['vy'][:,  neighbor_idx[:,P[i][v][3]][0], neighbor_idx[:,P[i][v][3]][1]][data_dict[urls[P[i][v][2]]]['index_sort']][data_dict[urls[P[i][v][2]]]['mask_dates']], dtype=np.float64) for v in range(len(P[i]))]
    
    # Grab the nodata value
    fillvalue = vxObs[0].min()
    
    # Mask the non-valid values (where vObs has no observations)
    mask = [np.logical_not(np.equal(vxObs[v], fillvalue)) for v in range(len(vxObs))]
    
    # Mask the non-valid inversion timestamps (where there is no observation)
    mask_inversion = data_dict[urls[P[i][0][2]]]['A_m'][mask[0]].sum(axis=0) == 0

    # Mask the observed vectors
    vxObs_masked = [vxObs[v][mask[v]] for v in range(len(vxObs))]
    vyObs_masked = [vyObs[v][mask[v]] for v in range(len(vyObs))]

    # Get the length of each point
    len_pts_tot =  [len(vxObs_masked[v]) for v in range(len(vxObs))]

    # Stack the observed vectors
    vxObs_masked = np.hstack(vxObs_masked)
    vyObs_masked = np.hstack(vyObs_masked)

    # Count the amount of points
    nb_pts = len(P[i])

    # Assemble the design matrix, masked 
    # Initialize the design matrix
    if spatial_regularization:
        A_des = np.zeros((len(vxObs_masked) + SRT[0] + SRS[0], SRT[1]))
        # Fill-in the values of A_m depending on each point
        cr = 0
        cc = 0
        for v in range(nb_pts):
            A_des[cr:cr+len_pts_tot[v], cc:cc+len_pt_inverted] = data_dict[urls[P[i][v][2]]]['A_m'][mask[v]]
            cr += len_pts_tot[v]
            cc += len_pt_inverted

        # Append time and spatial regularization matrices, fitted to the amount of points we have
        A_des[cr:cr + SRT[0], :] = time_reg_mat
        A_des[cr + SRT[0]:, :] = space_reg_mat

    else:
        A_des = np.zeros((len(vxObs_masked) + SRT[0], SRT[1]))
        cr = len_pts_tot[0]
        cc = len_pt_inverted
        A_des[:cr] = data_dict[urls[P[i][0][2]]]['A_m'][mask[0]]
        A_des[cr:cr + SRT[0], :] = time_reg_mat


    # Invert the velocities
    vxInv[: , y_matches[i], x_matches[i]] = Inverter(GPU, spatial_regularization, vxObs_masked, len_pt_inverted, A_des, device, mask_inversion, SRS, SRT)
    vyInv[: , y_matches[i], x_matches[i]] = Inverter(GPU, spatial_regularization, vyObs_masked, len_pt_inverted, A_des, device, mask_inversion, SRS, SRT)
    
    
    # Save the amount of iterations in a text file
    if i%100 == 0:
        with open("Counter.txt", "w") as text_file:
            text_file.write(f"Counter: {i}")

    # Save the matrices along-the-way in case the algorithm fails
    if i%10000 == 0 and i != 0:

        # Create output folder
        path_save = "Datacubes"
        os.makedirs(path_save, exist_ok=True)
    
        # Gather the data in a dataset
        new_ds = xr.Dataset(
        {
            "vx": (["time", "y", "x"], vxInv),
            "vy": (["time", "y", "x"], vyInv),
        },
        coords={
            "time": dates_nonum,
            "x": X_arr[X_MIN:X_MAX],
            "y": Y_arr[Y_MIN:Y_MAX]
        },
        attrs=data_dict[urls[0]]['zarr_store'].attrs,
        ).chunk({'time': 1, 'x': 100, 'y': 100})
    
        # Add the coordinates of the centerline
        new_ds['x_coords_centerline'] = ('point', x_indices_centerline)
        new_ds['y_coords_centerline'] = ('point', y_indices_centerline)
    
        from dask.diagnostics import ProgressBar
        save_data = new_ds.to_netcdf(f"{path_save}/{name_cube}.nc", mode='w', compute = False) 
        with ProgressBar():
            print(f"Writing to {path_save}")
            save_data.compute()
            
    return vxInv, vyInv
            
# Function to invert a point's velocities           
def Inverter(GPU, spatial_regularization, vObs_masked, len_pt_inverted, A_des, device, mask_inversion, SRS, SRT):

    # Extend the observations with 0s to match the length of spacetime regularizations
    if spatial_regularization:
        vObs_masked = np.hstack((vObs_masked, np.zeros((SRT[0] + SRS[0]))))
    else:
        vObs_masked = np.hstack((vObs_masked, np.zeros((SRT[0]))))

    if GPU:
        # Migrate velocity vector to torch
        vObs_masked = torch.from_numpy(vObs_masked).to(device).double()
        # Migrate design matrix to GPU
        A_des = torch.from_numpy(A_des).to(device)
        vInv = torch.linalg.solve(A_des.T@A_des,A_des.T@(vObs_masked))[:len_pt_inverted].cpu()
    else:
        #A_des = scipy.sparse.coo_array(A_des)
        #vInv = scipy.sparse.linalg.lsqr(A_des, vObs_masked)[0][:len_pt_inverted]
        #vInv = scipy.sparse.linalg.lsqr(A_des.T@A_des,A_des.T@(vObs_masked))
        vInv = np.linalg.solve(A_des.T@A_des,A_des.T@vObs_masked)[:len_pt_inverted]
        
    vInv[mask_inversion] = np.nan

    return vInv

# Save the dataset as a netcdf
def save_dataset(name_cube, vxInv, vyInv, dates_nonum, X_arr, Y_arr, X_MIN, X_MAX, Y_MIN, Y_MAX, x_indices_centerline, y_indices_centerline, data_dict, urls):
    
    
    # Calculate velocity magnitude
    mag = np.sqrt(vxInv**2+vyInv**2)

    # Plot the average velocity magnitude in time
    plt.scatter(dates_nonum, np.nanmean(np.nanmean(mag, axis = 2), axis = 1))

    # Create output folder
    path_save = "Datacubes"
    os.makedirs(path_save, exist_ok=True)

    # Gather the data in a dataset
    new_ds = xr.Dataset(
    {
        "vx": (["time", "y", "x"], vxInv),
        "vy": (["time", "y", "x"], vyInv),
    },
    coords={
        "time": dates_nonum,
        "x": X_arr[X_MIN:X_MAX],
        "y": Y_arr[Y_MIN:Y_MAX]
    },
    attrs=data_dict[urls[0]]['zarr_store'].attrs,
    ).chunk({'time': 1, 'x': 100, 'y': 100})

    # Add the coordinates of the centerline
    new_ds['x_coords_centerline'] = ('point', x_indices_centerline)
    new_ds['y_coords_centerline'] = ('point', y_indices_centerline)

    from dask.diagnostics import ProgressBar
    save_data = new_ds.to_netcdf(f"{path_save}/{name_cube}.nc", mode='w', compute = False) 
    with ProgressBar():
        print(f"Writing to {path_save}")
        save_data.compute()