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
from ipyleaflet import Map, WMSLayer, basemaps, GeoData, AwesomeIcon, Marker, Polygon
from ipywidgets import HTML, widgets
from owslib.wms import WebMapService
import ipywidgets as widgets
from ipywidgets import Label, VBox
from owslib.wfs import WebFeatureService
from requests import Request
import urllib.request, json 

from ITS_LIVE_TOOL import obj_setup

def create_glacier_obj(label, rgi_id, utm_zone):

    glacier = obj_setup.Glacier(label, rgi_id, utm_zone, 'manual')

    return glacier

def create_glacier_point_obj(point_coords, label, rgi_id):

    point = obj_setup.Glacier_Point(label, label, rgi_id, point_coords)
    return point


def create_glacier_centerline_obj(label, rgi_id):

    centerline = obj_setup.Glacier_Centerline(label, rgi_id)
    return centerline