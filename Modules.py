# Imports required

import subprocess
import sys
import os
import os, shutil
import glob
from glob import glob
from os import listdir
from os.path import exists

# Install the whl files
############################################################
sys.path.append("Required")
WHL_Dir = 'Required/'

try:
	import fiona
except:
	subprocess.check_call([sys.executable, '-m', 'pip', 'install', WHL_Dir+'Fiona-1.8.4-cp27-cp27m-win32.whl'])
import fiona
from fiona.crs import from_epsg

try:
	import gdal
except:
	subprocess.check_call([sys.executable, '-m', 'pip', 'install', WHL_Dir+'GDAL-2.2.4-cp27-cp27m-win32.whl'])
import gdal

try:
	import pyproj
except:
	subprocess.check_call([sys.executable, '-m', 'pip', 'install', WHL_Dir+'pyproj-1.9.6-cp27-cp27m-win32.whl'])
import pyproj

try:
	import rasterio
except:
	subprocess.check_call([sys.executable, '-m', 'pip', 'install', WHL_Dir+'rasterio-1.0.13-cp27-cp27m-win32.whl'])
import rasterio


# subprocess.check_call([sys.executable, '-m', 'pip', 'install', WHL_Dir+'Shapely-1.6.4.post2-cp27-cp27m-win32.whl'])
import shapely
from shapely import geometry
from shapely.geometry import shape
from shapely.geometry import Polygon, mapping
from shapely.geometry import Point

# subprocess.check_call([sys.executable, '-m', 'pip', 'install', WHL_Dir+'Rtree-0.8.3-cp27-cp27m-win32.whl'])
import rtree

#############################################################
try:
	import matplotlib
except:
	subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'matplotlib'])

import matplotlib.pyplot as plt

try:
	import numpy as np
except:
	subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'numpy'])

import numpy as np

try:
	import pandas as pd
except:
	subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'pandas'])

import pandas as pd


# GDAL (ORG included w/ GDAL) FIONA RASTERIO PYPROJ manual installation

import gdal
import ogr
import osr
import fiona
import rasterio
import rasterio.plot
import pyproj

# IMPORTS FOR MAPPING
import rasterio
from rasterio.transform import from_origin
from rasterio.warp import reproject, Resampling
import matplotlib.pyplot as plt
from rasterio.transform import from_bounds, from_origin
from rasterio.warp import reproject, Resampling
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.merge import merge
from rasterio.plot import show
from rasterio.warp import transform

try:
	import scipy
except:
	subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'scipy'])

import scipy

# try:
# 	import shapely
# except:
# 	subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'shapely'])

# import shapely

# from shapely import geometry
# from shapely.geometry import shape
# from shapely.geometry import Polygon, mapping
# from shapely.geometry import Point

# DOWNLOAD MAPS

##################################################
try:
	import rtree
except:
	subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'rtree'])
import rtree

try:
	import geopandas as gpd
except:
	subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'geopandas'])

import geopandas as gpd
from geopandas import GeoDataFrame

# try:
# 	import rtree
# except:
# 	subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'rtree'])

# import rtree


try:
	import tqdm
except:
	subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'tqdm'])

import tqdm

try:
	import requests
except:
	subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'requests'])

import requests

try:
	import collections
except:
	subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'collections'])

import collections

try:
	import tarfile
except:
	subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'tarfile'])

import tarfile

import datetime

import math

import urllib2,urllib

import copy
from copy import deepcopy

# Used for adding the online basemap
# Failed needs glob import
try:
	import cartopy
except:
	subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'WHL_Files/Cartopy-0.17.0-cp27-cp27m-win32.whl'])
import cartopy

try:
	import contextily as ctx
except:
	subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'contextily'])
import contextily as ctx

import shapefile

import time

import zipfile

import math

import csv

import collections

#### For google ee
import re

try:
  # if setuptools is available, use it to take advantage of its dependency
  # handling
  from setuptools import setup                          # pylint: disable=g-import-not-at-top
except ImportError:
  # if setuptools is not available, use distutils (standard library). Users
  # will receive errors for missing packages
  from distutils.core import setup                      # pylint: disable=g-import-not-at-top
import ee
import ee.mapclient
from ee import batch
import time

from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size

from matplotlib import cm
from matplotlib import colors
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

