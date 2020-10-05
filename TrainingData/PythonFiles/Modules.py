# Imports required

import subprocess
import sys
import os
import os, shutil
import glob
from glob import glob
from os import listdir
from os.path import exists

############################################
# start with basics for gdal and pandas
try:
	import numpy as np
except:
	subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'numpy'])

import numpy as np

try:
	import matplotlib
except:
	subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'matplotlib'])

import matplotlib.pyplot as plt

try:
	import pandas as pd
except:
	subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'pandas'])

import pandas as pd

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

# GDAL (ORG included w/ GDAL) FIONA RASTERIO PYPROJ manual installation

import ogr
import osr
import rasterio.plot
import pyproj

# IMPORTS FOR MAPPING
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

from scipy import stats

from scipy.stats import alpha
from scipy.stats import anglit
from scipy.stats import arcsine
from scipy.stats import argus
from scipy.stats import beta
from scipy.stats import betaprime
from scipy.stats import bradford
from scipy.stats import burr
from scipy.stats import burr12
from scipy.stats import cauchy
from scipy.stats import chi
from scipy.stats import chi2
from scipy.stats import cosine
from scipy.stats import crystalball
from scipy.stats import dgamma 
from scipy.stats import dweibull
from scipy.stats import erlang
from scipy.stats import expon
from scipy.stats import exponnorm
from scipy.stats import exponweib
from scipy.stats import exponpow
from scipy.stats import f
from scipy.stats import fatiguelife
from scipy.stats import fisk
from scipy.stats import foldcauchy
from scipy.stats import foldnorm
from scipy.stats import frechet_r
from scipy.stats import frechet_l
from scipy.stats import genlogistic
from scipy.stats import gennorm
from scipy.stats import genpareto
from scipy.stats import genexpon
from scipy.stats import genextreme
from scipy.stats import gausshyper
from scipy.stats import gamma
from scipy.stats import gengamma
from scipy.stats import genhalflogistic
# Will not work
# from scipy.stats import geninvgauss
from scipy.stats import gilbrat
from scipy.stats import gompertz
from scipy.stats import gumbel_r
from scipy.stats import gumbel_l
from scipy.stats import halfcauchy
from scipy.stats import halflogistic
from scipy.stats import halfnorm
from scipy.stats import halfgennorm
from scipy.stats import hypsecant
from scipy.stats import invgamma
from scipy.stats import invgauss
from scipy.stats import invweibull
from scipy.stats import johnsonsb
from scipy.stats import johnsonsu
from scipy.stats import kappa4
from scipy.stats import kappa3
from scipy.stats import ksone
from scipy.stats import kstwobign
from scipy.stats import laplace
from scipy.stats import levy
from scipy.stats import levy_l
from scipy.stats import levy_stable
from scipy.stats import logistic
from scipy.stats import loggamma
from scipy.stats import loglaplace
from scipy.stats import lognorm
from scipy.stats import lomax
from scipy.stats import maxwell
from scipy.stats import mielke
from scipy.stats import moyal
from scipy.stats import nakagami
from scipy.stats import ncx2
from scipy.stats import ncf
from scipy.stats import norm
from scipy.stats import norminvgauss
from scipy.stats import pareto
from scipy.stats import pearson3
from scipy.stats import powerlaw
from scipy.stats import powerlognorm
from scipy.stats import powernorm
from scipy.stats import rdist
from scipy.stats import rayleigh
from scipy.stats import rice
from scipy.stats import recipinvgauss
from scipy.stats import semicircular
from scipy.stats import skewnorm
from scipy.stats import t
from scipy.stats import trapz
from scipy.stats import triang
from scipy.stats import truncexpon
from scipy.stats import truncnorm
from scipy.stats import tukeylambda
from scipy.stats import uniform
from scipy.stats import vonmises
from scipy.stats import vonmises_line
from scipy.stats import wald
from scipy.stats import weibull_min
from scipy.stats import weibull_max
from scipy.stats import wrapcauchy
from scipy.stats import rv_continuous


import warnings

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

from datetime import timedelta
from datetime import datetime

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
from matplotlib import pyplot


import itertools

