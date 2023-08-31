import json
from shapely.geometry import MultiLineString, LineString, Point, MultiPoint
import numpy as np
import pandas as pd
import xarray as xr
import pyproj



# Function to split the geojson lines transects
def line_splitter(name_transect, delta_d):
    # read in the GeoJSON file
    with open(f'{name_transect}.geojson', 'r') as f:
        data = json.load(f)

    # create a shapely MultiLineString object from the GeoJSON data
    multilinestring = MultiLineString(data['features'][0]['geometry']['coordinates'])

    # initialize the list to store the coordinates of all points generated
    all_points_coords = []

    # iterate through each line in the MultiLineString
    for line_coords in multilinestring:
        # create a shapely LineString object from the line coordinates
        line = LineString(line_coords)

        # set the distance delta
        delta_d = 100
        
        # generate the equidistant points
        distances = np.arange(0, line.length, delta_d)
        points = [line.interpolate(distance) for distance in distances] + [line.boundary[1]]

        # append the coordinates of the points to the list of all points coordinates
        all_points_coords += [list(point.coords)[0] for point in points]

        # convert to a dataframe
        all_points_coords = pd.DataFrame(all_points_coords, columns=['x', 'y'])

        #CRS
        epsg = f"EPSG:{data['crs']['properties']['name'].split(':')[-1]}"
        print(f'transect projection: {epsg}')

    return all_points_coords, epsg




# Function to extract the points' values from the desired dataset
def extract_points_from_datacube(ds, variable, df_points, epsg_pts, amount, type):

    # Define a list of (lon, lat) coordinates to extract
    points_list = [(df_points.x.values[i], df_points.y.values[i]) for i in range(len(df_points.x.values))]

    # If epsg of points and dataset not the same, convert points
    if epsg_pts != f"EPSG:{ds.attrs['projection']}":

        # Create a pyproj transformer object to transform the coordinates
        transformer = pyproj.Transformer.from_crs(epsg_pts, f"EPSG:{ds.attrs['projection']}", always_xy=True)
        # Reproject the coordinates using the transformer object
        points_list = [transformer.transform(x, y) for x, y in points_list]


    # Extract the values at each coordinate and concatenate them into a single array
    points_values = []
    for xx, yy in points_list:
        points_values.append(ds[variable].sel(x=xx, y=yy, method='nearest').values)
    points_values = np.stack(points_values, axis=1)

    COL = np.array([np.absolute(ds.x.values-xx).argmin() for xx in [points_list[i][0] for i in range(len(df_points))]])
    ROW = np.array([np.absolute(ds.y.values-yy).argmin() for yy in [points_list[i][1] for i in range(len(df_points))]])
    # The resulting array will have shape (time, len(points_list))
    # Create a pandas DataFrame from the points_values array

    dist = [0]
    dist = np.array((np.append(dist, np.round(np.cumsum([np.sqrt((points_list[i][0]-points_list[i-1][0])**2+(points_list[i][1]-points_list[i-1][1])**2) for i in range(1,len(points_list))])))), dtype=int)

    # Create the dataframe, linearize its time axis every 5 days (lowest dt in the datacube), median over overlapping values, ignore nans and interpolate over time
    df = pd.DataFrame(points_values, index=pd.to_datetime(ds.time.values), columns=dist).groupby(pd.Grouper(freq=f'{amount}{type}')).median().interpolate('time')

    return df, ROW, COL
