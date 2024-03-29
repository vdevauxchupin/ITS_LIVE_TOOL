# %%
# AUTOGENERATED! DO NOT EDIT! File to edit: ../nbs/06_inversion.ipynb.

# %% auto 0
__all__ = ['mission', 'lamb', 'derivative', 'day_interval', 'sdate', 'edate', 'GPU', 'make_input_dict', 'create_data_dict',
           'get_extents', 'design_matrices', 'Inv_reg']

# %% ../nbs/06_inversion.ipynb 4
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
from scipy.interpolate import make_interp_spline
import os
from shapely.geometry import shape, Polygon, Point
import matplotlib.patches as patches
import time
from ITS_LIVE_TOOL import datacube_tools
dc = datacube_tools.DATACUBETOOLS()
try:
    import torch
except:
    print('torch import failed. Is it installed ?')
    
    

def convert_pt_epsg(pt_coords, source_crs, target_crs):
    # Create a transformer object to convert coordinates from source to target CRS
    transformer = pyproj.Transformer.from_crs(source_crs, target_crs, always_xy=True)
    # Project the coordinates of the corner from source to target CRS
    projected_coords = transformer.transform(pt_coords[0], pt_coords[1])
    
    return projected_coords


def subdivide_array(array, size_subarrays, X_template, Y_template):
    subdivided_arrays = []
    min_max_coords = []
    rows, cols = array.shape

    if rows % size_subarrays != 0:
        imax = rows//size_subarrays + 1
    else: 
        imax = rows/size_subarrays

    if cols % size_subarrays != 0:
        jmax = rows//size_subarrays + 1
    else: 
        jmax = rows/size_subarrays    

    for i in range(0, int(imax)):
        if i*size_subarrays < min((i+1)*size_subarrays, rows-1):
            for j in range(0, int(jmax)):
                if j*size_subarrays < min((j+1)*size_subarrays, cols-1):
                    sub_array = array[i*size_subarrays:min((i+1)*size_subarrays, rows-1), j*size_subarrays:min((j+1)*size_subarrays, cols-1)]
                    subdivided_arrays.append(sub_array)#.copy())
                    min_x, max_x = j*size_subarrays, min((j+1)*size_subarrays, cols-1)
                    min_y, max_y = i*size_subarrays, min((i+1)*size_subarrays, rows-1)
                    min_max_coords.append((min_y, min_x, max_y, max_x))

    min_max_coords_geo = []
    for coord in min_max_coords:
        coord = coord
        min_y, min_x, max_y, max_x = coord
        min_max_coords_geo.append((Y_template[min_y], X_template[min_x], Y_template[max_y], X_template[max_x]))

    return subdivided_arrays, min_max_coords, min_max_coords_geo



def plot_subdivisions(array, Xarr, Yarr, min_max_coords_geo, X_lims_geo, Y_lims_geo, valid_subcubes):
    plt.figure(figsize=(8, 8))
    plt.pcolormesh(Xarr, Yarr, array, cmap='viridis')
    for c in valid_subcubes:
        min_y, min_x, max_y, max_x = min_max_coords_geo[c]
        rect = patches.Rectangle((min_x, min_y), max_x - min_x, max_y - min_y, linewidth=1, edgecolor='r', facecolor='none')
        plt.gca().add_patch(rect)
    for X in X_lims_geo:
        plt.axvline(X, color = 'green')
    for Y in Y_lims_geo:
        plt.axhline(Y, color = 'green')
    plt.title('Subdivisions of Array')
    plt.axis('off')
    plt.show()

    
def select_satellites_dates(ds, arg_dict):
    """
    Select mid_date values from the dataset where satellite_img1 or satellite_img2 is in values,
    and between the corresponding start and end dates.
    
    Parameters:
        ds (xarray.Dataset): The dataset containing satellite_img1 and satellite_img2 variables.
        satellites (list): List of satellite values to select.
        sdates (list): List of start dates.
        edates (list): List of end dates.

    Returns:
        xarray.DataArray: Selected mid_date values.
    """
    # Convert date strings to datetime objects
    sdates = pd.to_datetime(arg_dict['sdate'])
    edates = pd.to_datetime(arg_dict['edate'])
    
    # Create an empty boolean array with the same length as mid_date
    selected_indices = xr.zeros_like(ds['mid_date'], dtype=bool)
    
    # Loop through each value and update the selection condition
    for satellite, sdate, edate in zip(arg_dict['mission'], sdates, edates):
        satellite_condition = ((ds['satellite_img1'] == satellite) | (ds['satellite_img2'] == satellite))
        date_condition = ((pd.to_datetime(ds['mid_date']) >= sdate) & (pd.to_datetime(ds['mid_date']) <= edate))
        selected_indices |= (satellite_condition & date_condition)
        
    ds = ds.sel(mid_date = selected_indices)
    # Use advanced indexing to select mid_date values corresponding to the selected indices
    selected_dates = ds.sortby('mid_date')
    
    return selected_dates



# Calculate the nansum of grouped indices
def custom_nansum(arr, indices):
    
    arrays = [arr[idx] for idx in indices]
    
    # Calculate sum ignoring NaNs
    summed = np.nansum(arrays, axis=0)
    
    # Create mask where all arrays have NaN
    all_nan_mask = np.all(np.isnan(arrays), axis=0)
    
    # Set NaN where all arrays have NaN
    summed[all_nan_mask] = np.nan
    
    return summed



def prepare_data(arg_dict):
    
    start_time = time.time()
    
    
    
    ###### ----- ######

    
    
    print('Calculating overlay glaciers/polygon......')
    
    
    # There are multiple cases of what happens next:
    # 1: User drew polygon intersecting with RGI glacier shape
    # 2: User did not draw any polygon so we take the RGI glacier shape
    # 3: User drew polygon over an area with data but no RGI glacier shape
    
    
    # Case 1
    try: 
        polygon = shape(arg_dict['data_map'].draw_control.last_draw.get('geometry')) # Grab the polygon from the map
        selection = gpd.GeoDataFrame(index=[0], crs='epsg:4326', geometry=[polygon]) # Convert to geopandas dataframe
        
        # For every glacier selected, calculate the overlay:
        overlays = []  # List to store the overlay GeoDataFrames

        for glacier in arg_dict['data_map'].added_glaciers:
            overlay = gpd.overlay(selection , glacier.to_crs('epsg:4326'), how='intersection')
            overlays.append(overlay)
        
        # Merge the intersection of glaciers and shape together
        merged_df = pd.concat(overlays, ignore_index=True)
        merged_gdf = gpd.GeoDataFrame(merged_df, crs='epsg:4326', geometry='geometry')
        
        # 'overlays' is now our new layer with which we work
        flag_overlay = True
    
    except:
        # Case 2
        try: 
            # If above fails, we just take the entire glaciers
            overlays = []  # List to store the overlay GeoDataFrames
            for glacier in arg_dict['data_map'].added_glaciers:
                overlay = glacier.to_crs('epsg:4326')
                overlays.append(overlay)
            # Merge the intersection of glaciers and shape together
            merged_df = pd.concat(overlays, ignore_index=True)
            merged_gdf = gpd.GeoDataFrame(merged_df, crs='epsg:4326', geometry='geometry')        
            flag_overlay = False
        
        # Case 3
        except:
            overlays = []  # List to store the overlay GeoDataFrames
            merged_gdf = selection # The correct shape is directly the one drew by the user      
            flag_overlay = False           
    
    # Iterate through the corners of the shape
    minx, miny, maxx, maxy = merged_gdf.total_bounds
    
    # List comprehension to create a list of corner coordinates
    corner_coordinates = [(minx, maxy), (maxx, maxy), (maxx, miny), (minx, miny)]
    
    end_eval_overlay = time.time()
    execution_time = round(end_eval_overlay - start_time, 2)
    print("......Execution time:", execution_time, "seconds")
    
    
    
    ####### ----- #######
    
    
    
    print("Fetch datacubes information......")

    # Fetch the URL for each corner coordinate and store in a list
    urls = [
        dc.find_datacube_catalog_entry_for_point((lon, lat), '4326')[0]['properties']['zarr_url']
        for lon, lat in corner_coordinates
    ]

    # Assuming urls is your list of URLs of each corner
    urls = [re.sub(r'http://|\.s3\.amazonaws\.com', '', url) for url in urls]

    # List to store unique 'x' and 'y' coordinates
    X_tot, Y_tot, X_valid, Y_valid, X_valid_index, Y_valid_index, min_date, max_date = [], [], [], [], [], [], [], []


    # Loop over each corner URL and its corresponding coordinates
    for url, pt in zip(urls, corner_coordinates):
        # Open the Zarr store
        fs = s3fs.S3FileSystem(anon=True)
        store = zarr.open(s3fs.S3Map(url, s3=fs))

        # Extract 'x' and 'y' arrays and projection attribute
        x, y = store['x'][:], store['y'][:]
        projection = store.attrs['projection']

        X_tot.append(x)
        Y_tot.append(y)

        # Convert the coordinates of the corner to the projection of the point
        projected_coords = convert_pt_epsg(pt, '4326', projection)

        # Find the 'x' and 'y' values to which the coordinates of the point correspond
        x_valid_index = np.abs(x - projected_coords[0]).argmin()
        y_valid_index = np.abs(y - projected_coords[1]).argmin()
        X_valid.append(x[x_valid_index])
        Y_valid.append(y[y_valid_index])
        

    # Determine the min and max of the dates available in the datacubes
    min_date = np.min([np.datetime64(sdate) for sdate in arg_dict['sdate']])
    max_date = np.max([np.datetime64(edate) for edate in arg_dict['edate']])

    # Create the timeframe on which we interpolate the dates for all the datacubes
    dates_total = np.arange(min_date, max_date+(arg_dict['days_interval']*2), timedelta(days=arg_dict['days_interval'])).astype(np.datetime64)

    
    end_eval_valid = time.time()
    execution_time = round(end_eval_valid - end_eval_overlay, 2)
    print("......Execution time:", execution_time, "seconds")

    # Determine the min and max of the dates available in the datacubes
    min_date = np.min(min_date)
    max_date = np.max(max_date)

    # Create the timeframe on which we interpolate the dates for all the datacubes
    dates_total = np.arange(min_date, max_date+(arg_dict['days_interval']*2), timedelta(days=arg_dict['days_interval'])).astype(np.datetime64)
    
    end_eval_valid = time.time()
    execution_time = round(end_eval_valid - end_eval_overlay, 2)
    print("......Execution time:", execution_time, "seconds")
    
    
    
    ###### ----- ######
    
    
    
    print("Starting creating the template (this step can take a few minutes) ......")

    # Create Eastings and Northings arrays based on the Eastings and Northings of all datacubes
    X_arr, Y_arr = np.unique(np.hstack(X_tot)), np.unique(np.hstack(Y_tot))

    # Calculate the overall minimum and maximum valid values across all datacubes
    x_min_all, x_max_all = np.min(X_valid), np.max(X_valid)
    y_min_all, y_max_all = np.min(Y_valid), np.max(Y_valid)

    # Define the overall boundaries of the region of interest
    X_MIN_all = max(np.searchsorted(X_arr, x_min_all) - 1, 0)
    X_MAX_all = min(np.searchsorted(X_arr, x_max_all) + 1, len(X_arr))
    Y_MIN_all = max(np.searchsorted(Y_arr, y_min_all) - 1, 0)
    Y_MAX_all = min(np.searchsorted(Y_arr, y_max_all) + 1, len(Y_arr))

    # Create the eastings and northings of the template
    Y_template = Y_arr[Y_MIN_all:Y_MAX_all]
    X_template = X_arr[X_MIN_all:X_MAX_all]

    # Create a template that gathers the shape of interest
    template = np.zeros((len(Y_template), len(X_template)))

    # Evaluate which points are included in the shape of interest
    Xs, Ys = np.meshgrid(X_template, Y_template)

    # Create an empty mask with the same shape as 'template'
    mask = np.zeros_like(template, dtype=bool)
    
    # Reproject the CRS of merged_gdf to projection
    merged_gdf = merged_gdf.to_crs(epsg=int(projection))

    # Iterate through the polygons in 'merged_gdf' and check if each point in 'template' falls within any polygon
    counter = 0
    for polygon in merged_gdf.geometry:
        for i in range(template.shape[0]):
            for j in range(template.shape[1]):
                counter += 1
    with tqdm(total=counter) as pbar:
        for polygon in merged_gdf.geometry:
            for i in range(template.shape[0]):
                for j in range(template.shape[1]):    
                    point = Point(Xs[i][j], Ys[i][j])
                    if polygon.contains(point):
                        mask[i][j] = True
                    pbar.update(1)
                    
    # Grab the indices of the points inside the glacier
    template = np.zeros(mask.shape) 
    template[mask==True] = 1
    template = template 
    
    X_lims = []
    Y_lims = []

    for c in range(len(urls)):
        x_indices = np.where(np.isin(X_template, X_tot[c]))
        y_indices = np.where(np.isin(Y_template, Y_tot[c]))
        X_lims.append((x_indices[0][0], x_indices[0][-1]))
        Y_lims.append((y_indices[0][0], y_indices[0][-1]))
    X_lims = np.unique(X_lims)
    Y_lims = np.unique(Y_lims)

    # Grab the indices of values that are only 1 apart
    X_lims = np.delete(X_lims, np.where(np.diff(X_lims) == 1))
    Y_lims = np.delete(Y_lims, np.where(np.diff(Y_lims) == 1))

    # Filter out 0s and values equal to template.shape - 1
    X_lims = X_lims[(X_lims != 0) & (X_lims != X_template.shape[0] - 1)]
    Y_lims = Y_lims[(Y_lims != 0) & (Y_lims != Y_template.shape[0] - 1)]

    X_lims_geo = [X_template[X] for X in X_lims]
    Y_lims_geo = [Y_template[Y] for Y in Y_lims]
    
    end_eval_ice = time.time()
    execution_time = round(end_eval_ice - end_eval_valid, 2)
    print("......Execution time:", execution_time, "seconds")
    
    
    
    ###### ----- ######
    
    
    
    print("Starting subdividing the area of interest......")
    
    # Subdivide template
    subdivided_arrays, min_max_coords, min_max_coords_geo = subdivide_array(template, arg_dict['size_subarrays'], X_template, Y_template)

    
    # Calculate coordinates of the subdivisions in template, subarrays, and datacube referentials
    coords_templates = [] # List to store corner coordinates of subarrays
    coords_subarrays_geo = []
    coords_subarrays = []
    cube_belong = [] # List to store datacube ownership of subarray
    for coords in min_max_coords:
        ymin, xmin, ymax, xmax = coords
        coords_y = [ymin] + [ymax]
        coords_x = [xmin] + [xmax]
        for y_lim in Y_lims:
            if ymin <= y_lim <= ymax:
                coords_y = coords_y + [y_lim]

        for x_lim in X_lims:
            if xmin <= x_lim <= xmax:    
                coords_x = coords_x + [x_lim]

        # Sort the corner coordinates
        coords_x.sort()
        coords_y.sort()

        corners = []
        corners_geo = []
        corners_subarrays = []
        cube_centroid = []
        for ymin, ymax in zip(coords_y[:-1], coords_y[1:]):
            for xmin, xmax in zip(coords_x[:-1], coords_x[1:]):
                corners.append([ymin, xmin, ymax, xmax]) # append the list of subarray corners
                corners_subarrays.append([ymin-coords_y[0],xmin-coords_x[0],ymax-coords_y[0], xmax-coords_x[0]])
                corners_geo.append([Y_template[ymin], X_template[xmin], Y_template[ymax-1], X_template[xmax-1]]) # Calculate the corners in geo
                centroid_x = X_template[(xmax+xmin)//2] # Calculate the coordinate of subarray centroid in geo
                centroid_y = Y_template[(ymax+ymin)//2]
                url_centroid = dc.find_datacube_catalog_entry_for_point((centroid_x, centroid_y), projection)[0]['properties']['zarr_url'] # Get datacube containing centroid
                url_centroid = re.sub(r'http://|\.s3\.amazonaws\.com', '', url_centroid)
                index_centroid = urls.index(url_centroid) # Grab the index in list of urls corresponding to datacube
                cube_centroid.append(index_centroid) # Append the list of datacubes subarrays

        cube_belong.append(cube_centroid)       
        coords_templates.append(corners)
        coords_subarrays_geo.append(corners_geo)
        coords_subarrays.append(corners_subarrays)

    # We first check which subcubes have no ice, and we squeeze them out of the arrays
    valid_subcubes = []
    valid_pixels = 0

    for sa in range(len(subdivided_arrays)):
        subarray = subdivided_arrays[sa] # Grab a subarray
        corners = coords_templates[sa] # Grab corners
        corners_geo = coords_subarrays_geo[sa] # Grab corners geo
        corners_subarrays = coords_subarrays[sa]
        urls_subarray = cube_belong[sa] # Grab urls in subarray
        flags = []
        for c in range(len(urls_subarray)):
            cc = corners_subarrays[c] # Get coordinates of the subarray
            cg = corners_geo[c] # Get coordinates of the subarray in geo
            ins3xr, small_ins3xr, bbox_centrer_point_cubexy = dc.get_subcube_for_bounding_box([cg[1], cg[0], cg[3], cg[2]], projection, variables=['landice']) # Grab ice mask its_live
            subarray[cc[0]:cc[2],cc[1]:cc[3]] = small_ins3xr.landice.values*subarray[cc[0]:cc[2],cc[1]:cc[3]]
            flags.append(1 in subarray[cc[0]:cc[2],cc[1]:cc[3]])
            valid_pixels += np.count_nonzero(subarray[cc[0]:cc[2],cc[1]:cc[3]] == 1)

        if any(flags):
            valid_subcubes.append(True)
        else:
            valid_subcubes.append(False)

    valid_subcubes = [index for index, value in enumerate(valid_subcubes) if value]
    
    plt.figure()
    plt.imshow(template)

    plot_subdivisions(template, X_template, Y_template, min_max_coords_geo, X_lims_geo, Y_lims_geo, valid_subcubes)
    
    end_subdiv = time.time()
    execution_time = round(end_subdiv - end_eval_ice, 2)
    print("......Execution time:", execution_time, "seconds")
    
    
    var_dict = {
    'template': template,
    'X_template': X_template,
    'Y_template': Y_template, 
    'subarrays': subdivided_arrays,
    'coords_template': coords_templates,
    'coords_subarrays_geo': coords_subarrays_geo,
    'coords_subarrays': coords_subarrays,
    'cube_belong': cube_belong, 
    'urls': urls,
    'valid_subcubes': valid_subcubes,
    'dx': np.abs(X_template[1]-X_template[0]),
    'timestep': timedelta(days=arg_dict['days_interval']),
    'regular_dates': dates_total,
    'projection': projection,
    'valid_pixels': valid_pixels
    }
    
    return var_dict
    

    
    
def prepare_subcube(subcube_number, var_dict, arg_dict):
    # Grab necessary data from pre-defined arrays
    subarray = var_dict['subarrays'][subcube_number]  # Grab a subarray
    corners = var_dict['coords_template'][subcube_number]  # Grab corners
    corners_geo = var_dict['coords_subarrays_geo'][subcube_number]  # Grab corners geo
    corners_subarrays = var_dict['coords_subarrays'][subcube_number]
    urls_subarray = var_dict['cube_belong'][subcube_number]  # Grab urls in subarray
    subcubes = []

    ### --- Grab data and sort it --- ###
    for c in range(len(urls_subarray)):
        cc = corners_subarrays[c]  # Get coordinates of the subarray
        cg = corners_geo[c]  # Get coordinates of the subarray in geo
        ins3xr, small_ins3xr, bbox_centrer_point_cubexy = dc.get_subcube_for_bounding_box([cg[1], cg[0], cg[3], cg[2]], var_dict['projection'], variables=['v', 'vx', 'vy', 'mid_date', 'acquisition_date_img1', 'acquisition_date_img2', 'satellite_img1', 'satellite_img2',])  # Grab ice mask its_live
        selected_data = select_satellites_dates(small_ins3xr, arg_dict) # Keep only dates in the date range and with the correct satellite ID
        nan_indices = np.where(np.all(np.isnan(selected_data['v']), axis=(1, 2)))[0]  # Get indices where slices are NaNs entirely
        indices = np.arange(0, len(selected_data.mid_date))  # Get indices of mid_date
        indices_to_keep = np.setdiff1d(indices, nan_indices)  # Get indices of mid_date that are not full of NaNs
        selected_data = selected_data.isel(mid_date=indices_to_keep)  # Remove entries at mid_date indices where all the slices are nans
        subcubes.append(selected_data)

    # Concatenate mid_dates from all subcubes
    mid_dates_tot = np.concatenate([dataset.mid_date.values for dataset in subcubes])
    im1_tot = np.concatenate([dataset.acquisition_date_img1.values for dataset in subcubes])
    im2_tot = np.concatenate([dataset.acquisition_date_img2.values for dataset in subcubes])

    # Create arrays filled with NaNs for vxInv_tot, vyInv_tot, vInv_tot
    vxInv_tot = np.full((len(mid_dates_tot), subarray.shape[0], subarray.shape[1]), np.nan)
    vyInv_tot = np.full((vxInv_tot.shape), np.nan)
    vInv_tot = np.full((vxInv_tot.shape), np.nan)

    # Fill-up the arrays
    time_inds = 0
    for c in range(len(subcubes)):
        ymin, xmin, ymax, xmax = corners_subarrays[c]
        lt = subcubes[c].mid_date.shape[0]
        vxInv_tot[time_inds:lt + time_inds, ymin:ymax, xmin:xmax] = subcubes[c].vx.values
        vyInv_tot[time_inds:lt + time_inds, ymin:ymax, xmin:xmax] = subcubes[c].vy.values
        vInv_tot[time_inds:lt + time_inds, ymin:ymax, xmin:xmax] = subcubes[c].v.values
        time_inds += lt
    
    # Get unique mid_dates and their indices
    unique_values, unique_indices, inverse_indices = np.unique(mid_dates_tot, return_index = True, return_inverse=True)
    # Group together equal mid_dates indices
    grouped_indices = [(np.where(mid_dates_tot == mid_date)[0]) for mid_date in unique_values]
    
    # Create host arrays
    vxInv_host = np.full((len(unique_values),  subarray.shape[0], subarray.shape[1]), np.nan)
    vyInv_host = vxInv_host.copy()
    vInv_host = vxInv_host.copy()
    
    # Assemble the values from all the datacubes
    for t in range(len(unique_values)):
        vxInv_host[t] = custom_nansum(vxInv_tot, grouped_indices[t])
        vyInv_host[t] = custom_nansum(vyInv_tot, grouped_indices[t])
        vInv_host[t] = custom_nansum(vInv_tot, grouped_indices[t])
    
    # Sort the arrays based on mid_dates
    sorted_indices = np.argsort(unique_values)
    mid_dates_tot = unique_values[sorted_indices]
    vxInv_host = vxInv_host[sorted_indices]
    vyInv_host = vyInv_host[sorted_indices]
    vInv_host = vInv_host[sorted_indices]
    im1_tot = im1_tot[unique_indices][sorted_indices]
    im2_tot = im2_tot[unique_indices][sorted_indices]

    # Create dictionary with arrays for inversion
    subarray_dict = {
        'vxInv_tot': vxInv_host,
        'vyInv_tot': vyInv_host,
        'vInv_tot': vInv_host,
        'im1': im1_tot,
        'im2': im2_tot,
        'mid_date': mid_dates_tot
    }
    
    # Create dictionnary for indices
    indices_dict = {
        'subarray': var_dict['subarrays'][subcube_number] ,  # Grab a subarray
        'corners': var_dict['coords_template'][subcube_number] ,  # Grab corners
        'corners_geo': var_dict['coords_subarrays_geo'][subcube_number],  # Grab corners geo
        'corners_subarrays': var_dict['coords_subarrays'][subcube_number],
        'urls_subarray': var_dict['cube_belong'][subcube_number],  # Grab urls in subarray
    }

    return subarray_dict, indices_dict    
    
    
    
def likely_cutout(xp, yp, var_dict, arg_dict):
    
    for idx, coord in enumerate(var_dict['coords_template']):
        ymin, xmin, ymax, xmax = coord[0]
        if (xmin <= xp < xmax) and (ymin <= yp < ymax):
            print("The point belongs to the subarray:", idx)
            break
        
    subarray_dict, indices_dict = prepare_subcube(idx, var_dict, arg_dict) # Get ITS_LIVE data
    i_ice, j_ice = np.where(indices_dict['subarray'] == 1) # Get indices of on-ice pixels
    ymin_inv = indices_dict['corners'][0][0] # Get y index of the subarray considered in the template reference
    xmin_inv = indices_dict['corners'][0][1] # Get x index of the subarray considered in the template reference

    xp += -xmin_inv
    yp += -ymin_inv

    neighbors = grab_cross_around_pixel(subarray_dict, indices_dict['subarray'], yp, xp)

    # Create dataframe
    df = pd.DataFrame({'v': neighbors['v'][0],
                       'time': subarray_dict['mid_date']})

    # Calculate baseline
    df['baseline'] = np.array([int((subarray_dict['im2'][i] - subarray_dict['im1'][i]) / 1e9 / 3600 / 24)
                 for i in range(len(subarray_dict['im2']))
                ])

    # Throw out NaNs
    df = df.dropna(subset=['v'])
    
    # Calculate cumulative averages
    unique_baseline_values = np.unique(df['baseline'])
    cumulative_averages = [np.mean(df.loc[df['baseline'] <= i, 'v']) for i in unique_baseline_values]
    derivative = np.diff(cumulative_averages) / np.diff(unique_baseline_values) # Calculate diff of curve derivative

    plt.figure(figsize=(8, 6))
    plt.plot(unique_baseline_values, cumulative_averages, marker='o', linestyle='-', color='b')
    plt.xlabel('Unique Baseline Values')
    plt.ylabel('Cumulative Average of v')
    plt.title('Cumulative Average of v by Baseline')
    # Finding the index of the most negative value in the derivative array
    most_negative_index = np.argmin(derivative)
    # Plotting the red vertical line
    plt.axvline(x=unique_baseline_values[most_negative_index], color='r', linestyle='--')
    plt.grid(True)
    plt.show()
    print(f'Most likely cutout for baseline of {unique_baseline_values[most_negative_index]} days')

    return df, neighbors, subarray_dict
    
    
    
    
    
    
    
    
    
def grab_cross_around_pixel(data_dict, subarray, i, j):
    time, rows, cols = data_dict['vxInv_tot'].shape
    neighbors_vx = [data_dict['vxInv_tot'][:, i, j]]
    neighbors_vy = [data_dict['vyInv_tot'][:, i, j]]
    neighbors_v = [data_dict['vInv_tot'][:, i, j]]
    directions = [(0, -1), (0, 1), (-1, 0), (1, 0)]  # (row_change, col_change)
    
    for di, dj in directions:
        ni, nj = i + di, j + dj
        if 0 <= ni < rows and 0 <= nj < cols and subarray[ni, nj] != 0:
            neighbors_vx.append(data_dict['vxInv_tot'][:, ni, nj])
            neighbors_vy.append(data_dict['vyInv_tot'][:, ni, nj])
            neighbors_v.append(data_dict['vInv_tot'][:, ni, nj])
        else:
            neighbors_vx.append(np.zeros((time)))
            neighbors_vy.append(np.zeros((time)))
            neighbors_v.append(np.zeros((time)))
    
    neighbors = {
                'vx': neighbors_vx,
                'vy': neighbors_vy,
                'v': neighbors_v
            }
    
    return neighbors


def Interpolator(arg_dict, var_dict, neighbors, subarray_dict, lambda_space, lambda_time, pt = False):
    
    vx_tot = []
    vy_tot = []
    im1_tot = []
    im2_tot = []
    mid_dates_tot = []
    v_tot = []
    
    for point in range(5):
        
        # Create a mask that will mask out: nans, dates outside of range
        nanmask = ~np.isnan(neighbors['v'][point])
        
        # Create a mask taking out erratic values based on cutouts from user
        baseline_days = np.array(
            [int(
                (subarray_dict['im2'][nanmask][i] - subarray_dict['im1'][nanmask][i]) / 1e9 / 3600 / 24)
             for i in range(len(subarray_dict['im2'][nanmask]))
            ])
        baseline_mask = np.where(np.logical_and(baseline_days > var_dict['cutout_min'], baseline_days < var_dict['cutout_max']))[0]

        # Store the point values in the lists
        im1_tot.append(subarray_dict['im1'][nanmask][baseline_mask])
        im2_tot.append(subarray_dict['im2'][nanmask][baseline_mask])
        vx_tot.append(neighbors['vx'][point][nanmask][baseline_mask])
        vy_tot.append(neighbors['vy'][point][nanmask][baseline_mask])
        v_tot.append(neighbors['v'][point][nanmask][baseline_mask])
        mid_dates_tot.append(subarray_dict['mid_date'][nanmask][baseline_mask])
    
    # Now that we have all the points, stack them
    im1_tot = np.concatenate(im1_tot)
    im2_tot = np.concatenate(im2_tot)
    
    # Reference epoch
    epoch =  np.datetime64('1970-01-01T00:00:00Z')
    
    # Amount of Seconds in 1 day
    hour_year = (24*365)

    # Convert to numerical
    
    im1_num = ((np.array(im1_tot, dtype='datetime64[h]') - epoch)/np.timedelta64(1, 'h')/hour_year).astype(np.float64) # We divide by 1s because epoch is in seconds
    im2_num = ((np.array(im2_tot, dtype='datetime64[h]') - epoch)/np.timedelta64(1, 'h')/hour_year).astype(np.float64) 
    
    # Generate the SBAS pairs
    im_unique = np.unique(np.append(im1_num, im2_num))
    im_unique = np.sort(im_unique)
    
    # Calculate the time difference between each SBAS pair
    diff_im = np.diff(im_unique)


    # Dates
    dates = (im_unique[:-1] + im_unique[1:])/2

    # Dates at date format
    dates_nonum = (dates*hour_year*np.timedelta64(1, 'h'))+epoch
    im_unique_nonum = (im_unique*hour_year*np.timedelta64(1, 'h'))+epoch
    
    # Check if arrays are empty. If yes, then computation will fail
    #if all(not arr.size for arr in vx_tot) or all(not arr.size for arr in vy_tot) == True:
    if len(vx_tot[0]) < 4:
        flag_empty_input = True
        vxInv = np.full((var_dict['regular_dates'].shape[0]), np.nan)
        vyInv = np.full((var_dict['regular_dates'].shape[0]), np.nan)
        return vxInv, vyInv
        
    else:        
        # --------------- DESIGN MATRICES --------------- 


        # The design matrix of each point is composed of A_SBAS, reg_mat_time, reg_mat_space
        # Only A_SBAS can change from point to point. The rest are the same

        # Calculate the SBAS matrices
        A_SBAS_tot = [] # intialize storage
        d_v_tot = []
        d_vx_tot = []
        d_vy_tot = []

        # Initialize a counter
        c = 0

        # loop over each point in the cross 
        for p in range(5):

            # Initialize matrix
            A_SBAS = np.zeros((vx_tot[p].shape[0], dates_nonum.shape[0]), dtype=np.float64) # Depends on the length of the velocity vector
            d_v = np.zeros((v_tot[p].shape[0], 1), dtype=np.float64)
            d_vx = np.zeros((vx_tot[p].shape[0], 1), dtype=np.float64)
            d_vy = np.zeros((vy_tot[p].shape[0], 1), dtype=np.float64)

            # Fill-in the temporal overlap with inversion temporal grid
            for i in range(A_SBAS.shape[0]):

                ac1 = np.where(im1_num[c+i] == im_unique)[0][0] # Select index where im1 intersects unique dates
                ac2 = np.where(im2_num[c+i] == im_unique)[0][0] # Select index where im2 intersects unique dates
                A_SBAS[i,  ac1:ac2] = diff_im[ac1:ac2]        # Fill-in the SBAS matrix: plug-in the differences for which the image overlaps the inversion temporal grid
                d_v[i] = v_tot[p][i] * (im2_num[c+i] - im1_num[c+i])
                d_vx[i] = vx_tot[p][i] * (im2_num[c+i] - im1_num[c+i])
                d_vy[i] = vy_tot[p][i] * (im2_num[c+i] - im1_num[c+i])


            # Update the counter
            c += vx_tot[p].shape[0]

            A_SBAS_tot.append(A_SBAS)
            d_v_tot.append(d_v)
            d_vx_tot.append(d_vx)
            d_vy_tot.append(d_vy)

        # Initialize regularization matrices (time and space)
        reg_mat_space = []

        # Initialize centered differences regularization matrix
        reg_mat_time = np.zeros((A_SBAS.shape[1] - 2, A_SBAS.shape[1]))
        reg_mat_space.append(np.diag(np.diag(np.ones((A_SBAS.shape[1], A_SBAS.shape[1]))*(-4*lambda_space/(var_dict['dx'])), 0))) # Negative
        reg_mat_space.append(np.diag(np.diag(np.ones((A_SBAS.shape[1], A_SBAS.shape[1]))*(lambda_space/(var_dict['dx'])), 0))) # Positive

        # We then apply central differences
        for j in range(A_SBAS.shape[1]-2):
            reg_mat_time[j, j] = lambda_time * 2 / (diff_im[j] * (diff_im[j] + diff_im[j + 1]))
            reg_mat_time[j, j + 1] = -lambda_time * 2 / (diff_im[j] * diff_im[j + 1])
            reg_mat_time[j, j + 2] = lambda_time * 2 / (diff_im[j + 1] * (diff_im[j] + diff_im[j + 1]))



        # --------------- Assemble the design matrix for the inversion ---------------
        nb_points = len(A_SBAS_tot)
        # First, grab the shape of all the matrices we need
        rows_A_SBAS_tot = sum(arr.shape[0] for arr in A_SBAS_tot) # Amount of rows for all the SBAS matrices
        rows_time = nb_points * reg_mat_time.shape[0] # Amount of rows for all time regularization matrices
        cols_SBAS = A_SBAS_tot[0].shape[1] # Amount of columns for an SBAS matrix
        cols_des = nb_points * cols_SBAS # Total amount of columns total

        if arg_dict['spatial_regularization']:
            rows_space = reg_mat_space[0].shape[0] # If spatial regularization grab the amount of rows
        else:
            rows_space = 0

        A_des = np.zeros((rows_A_SBAS_tot + rows_time + rows_space, cols_des)) # Initialize the design matrix
        latest = 0 # Initialize a counter

        # Iterate through all the points
        for p in range(nb_points):
            A_des[latest: latest + A_SBAS_tot[p].shape[0], p*cols_SBAS: (p+1)*cols_SBAS] = A_SBAS_tot[p] # Plug-in the SBAS matrix of the point p
            latest += A_SBAS_tot[p].shape[0] # Update the counter
            A_des[rows_A_SBAS_tot + reg_mat_time.shape[0]*p : rows_A_SBAS_tot + reg_mat_time.shape[0]*(p+1), cols_SBAS*p : cols_SBAS*(p+1)] = reg_mat_time # Plug-in temporal matrix

        if arg_dict['spatial_regularization']:
            # Initialize spatial regulization matrix
            space_reg_mat = np.zeros((reg_mat_space[0].shape[0], cols_des)) # Initialize spatial reg matrix
            space_reg_mat[:, :reg_mat_time.shape[1]] = reg_mat_space[0] # Plug-in coefficients for point 1
            space_reg_mat[:, reg_mat_time.shape[1]:] = np.tile(reg_mat_space[1], nb_points-1) # Plug-in coefficients for the other points
            A_des[rows_A_SBAS_tot + rows_time:, :] = space_reg_mat # Plug-in spatial reg matrix    

        # --------------- Invert the velocities --------------- 
        start = time.time()
        vxInv, vx = Inverter(arg_dict['sparse'], d_vx_tot, arg_dict['GPU'], arg_dict['device'], A_des, arg_dict['spatial_regularization'], var_dict['regular_dates'], dates_nonum, rows_time, rows_space, A_SBAS_tot)
        vyInv, vy = Inverter(arg_dict['sparse'], d_vy_tot, arg_dict['GPU'], arg_dict['device'], A_des, arg_dict['spatial_regularization'], var_dict['regular_dates'], dates_nonum, rows_time, rows_space, A_SBAS_tot)
        #vInv = Inverter(d_v_tot)
        #print(time.time()-start)

        # --------------- Remap the SBAS onto a regular grid with spline interpolation --------------- 

        # This snippet looks at the correspondance of image pairs compared to the regular grid.
        # We look at which date index in the regular grid is the closest to a date of SBAS, less than one timestep away.
        # If we can't find one, we simply leave it as NaN.
        # This ensures the spline interpolation does NOT put a value when we have no measurement available.

        # Initialize an array to NaN dates without SBAS
        nan_coeffs = np.full(var_dict['regular_dates'].shape, np.nan)

        # For each entry in dates_nonum
        for i in range(len(im_unique_nonum[1:])):

            # Find the absolute difference between the current entry and all entries in array1
            diffs = np.abs(var_dict['regular_dates'] - im_unique_nonum[i])

            # Find the index of the closest entry in array1 that is less than a timestep away
            inds = np.where(diffs <= var_dict['timestep'])[0]

            # If such an index exists
            if inds.size > 0:
                # Find the index of the minimum difference
                min_ind = inds[np.argmin(diffs[inds])]

                # Put a '1' in nan_coeffs at the corresponding index
                nan_coeffs[min_ind] = 1

        # Mask the interpolated values
        vxInv = vxInv * nan_coeffs
        vyInv = vyInv * nan_coeffs    

        if pt == False:
            return vxInv, vyInv
        else:
            return vxInv, vyInv, v_tot[0], mid_dates_tot[0]
    
    
    
def Inverter(sparse, velocity, GPU, device, A_des, spatial_regularization, regular_dates, dates_nonum, rows_time, rows_space, A_SBAS_tot):
    
    # Velocities to invert
    vObs_masked = np.concatenate(velocity)

    # Extend the observations with 0s to match the length of spacetime regularizations
    if spatial_regularization:
        vObs_masked = np.vstack((vObs_masked, np.zeros((rows_time + rows_space, 1))))
    else:
        vObs_masked = np.vstack((vObs_masked,  np.zeros((rows_time, 1))))
        #vObs_masked = vObs_masked
    
    if GPU:
        # Migrate velocity vector to torch
        vObs_masked = torch.from_numpy(vObs_masked).to(device).double()
        # Migrate design matrix to GPU
        A_des = torch.from_numpy(A_des).to(device)
        vInv = torch.linalg.solve(A_des.T@A_des,A_des.T@(vObs_masked))[:A_SBAS_tot[0].shape[1],0].cpu()
    else:
        if sparse:
            A_des_sparse = scipy.sparse.lil_array(A_des)
            A_des_sparse = scipy.sparse.coo_array(A_des_sparse)
            vInv = scipy.sparse.linalg.lsqr(A_des_sparse, vObs_masked)
            vInv = vInv[0][:A_SBAS_tot[0].shape[1]]
        else:
            vInv = np.linalg.solve(A_des.T@A_des,A_des.T@vObs_masked)[:A_SBAS_tot[0].shape[1],0]
            
    if vInv.shape[0] < 4:
        closest_indices = np.argmin(np.abs(regular_dates[:, None] - dates_nonum), axis=0)
        interpolated_values = np.full(len(regular_dates), np.nan)
        interpolated_values[closest_indices] = vInv
    else:       
        # Interpolate on a regular grid
        # Create a cubic spline function
        spline_func = make_interp_spline(dates_nonum.astype('int64'), vInv)

        # Now you can use cubic_spline_func to interpolate at new_dates
        interpolated_values = spline_func(regular_dates.astype('datetime64[s]').astype('int64')) # we convert regular_dates to the same format as dates_nonum
        
        # Mask any speed higher than 20km/yr (very unlikely)
        interpolated_values[(interpolated_values > 2e4) | (interpolated_values < -2e4)] = np.nan
        
    return interpolated_values, vInv




def save_dataset(plotting, arg_dict, var_dict, vxInv, vyInv):
    # Plot the average velocity magnitude in time
    if plotting == True:
        # Calculate velocity magnitude
        mag = np.sqrt(vxInv**2+vyInv**2)
        plt.scatter(var_dict['regular_dates'], np.nanmean(np.nanmean(mag, axis = 2), axis = 1))

    # Create output folder
    path_save = f"Datacubes/{arg_dict['name_cube']}"
    os.makedirs(path_save, exist_ok=True)

    # Gather the data in a dataset
    new_ds = xr.Dataset(
    {
        "vx": (["time", "y", "x"], vxInv),
        "vy": (["time", "y", "x"], vyInv),
        "glacier_mask": (["y", "x"], var_dict['template']),
    },
    coords={
        "time": var_dict['regular_dates'],
        "x": var_dict['X_template'],
        "y": var_dict['Y_template']
    },
    attrs={"urls": var_dict['urls'], "projection": var_dict['projection']},
    ).chunk({'time': 100, 'x': 100, 'y': 100})

    from dask.diagnostics import ProgressBar
    save_data = new_ds.to_netcdf(f"{path_save}/{arg_dict['name_cube']}.nc", mode='w', compute = False) 
    with ProgressBar():
        print(f"Writing to {path_save}")
        save_data.compute()
