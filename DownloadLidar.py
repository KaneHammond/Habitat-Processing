import Modules
from Modules import *

#####################################################################################
# Path appending and creation

sys.path.append("Shapefile")
SF_Dir = 'Shapefile/'

dir = 'DEM_Output'
if not os.path.exists(dir):
    os.makedirs(dir)
Out_Dir = 'DEM_Output/'

dir = 'DEM_Output/Ph4QL3/'
if not os.path.exists(dir):
    os.makedirs(dir)
# Define output paths for this QL/Phase
f1out = 'DEM_Output/Ph4QL3/'

dir = 'DEM_Output/Ph5QL3/'
if not os.path.exists(dir):
    os.makedirs(dir)
# Define output paths for this QL/Phase
f2out = 'DEM_Output/Ph5QL3/'

dir = 'DEM_Output/LiDarData/'
if not os.path.exists(dir):
    os.makedirs(dir)
# Define output path for raw data
Lidar_Out = 'DEM_Output/LiDarData/'

dir = 'DEM_Output/ZipImg'
if not os.path.exists(dir):
    os.makedirs(dir)
DzipOut = 'DEM_Output/ZipImg/'

#####################################################################################
# Open shapefiles containing Lidar MetaData and prep geopandasdataframe

file1 = gpd.read_file(SF_Dir+"IndexLiDARJamesRiverPh4QL3.shp")

file2 = gpd.read_file(SF_Dir+"IndexLiDARJamesRiverPh5QL3.shp")

# Convert to different crs
file1 = file1.to_crs(epsg=3857)
file2 = file2.to_crs(epsg=3857)

# Write a dataframe with both data sets
df = pd.DataFrame(file1)
df2 = pd.DataFrame(file2)
# Combine the 2
All = pd.concat([df, df2])
# Reset the index
All = All.reset_index(drop=True)
# Write a new GPD item
gdf = gpd.GeoDataFrame(All, geometry = All.geometry)
# Set a crs since this is being created from concat df
# gdf.crs = {'init' :'epsg:3857'}
# fiona correctly sets crs, setting by dictionary results in error ** USE FIONA **
gdf.crs = fiona.crs.from_epsg(3857)
###################################################################################

## VIEW BASIC INFO FOR SHAPEFILE
## VIEW BASIC INFO FOR SHAPEFILE

# print file1.head()
# print file1.columns.values

# RETURNED (file1.columns.values):
# [u'NAME' u'Project' u'Phase' u'Block' u'TileNum' u'Filename' u'Path'
#  u'QL_Level' u'Contours' u'LAS_Class' u'H_Project' u'H_Datum' u'V_Datum'
#  u'Vert_RMSEz' u'Pt_Density' u'NPS' u'DEM_PstSp' u'MinCntrInt'
#  u'Start_Date' u'End_Date' u'ASCII_Path' u'tilename' u'Comments'
#  u'SRC_SRS' 'geometry']

# file1 and file2 will share same data info. The gdf is the concat of the 2
# completed by use of dataframe conversion to reduce dependancies on other imports

######################################################################################

# VIEW THE DATASET

# Plot enable section for basic mapping
Plotf1 = 'True'
# Plotf1 = 'T'
Plotf2 = 'True'
# Plotf2 = 'T'
Plot_g = 'True'
# Plot_g = 'T'
###############

def Plotfile1(Plotf1):
	if Plotf1=='True':
		# Plot file1
		## Will show the shapefile extent plain without basemap
		file1N = 'James River Ph4 QL3'
		# create figure and axes for Matplotlib
		fig, ax = plt.subplots(1, figsize=(10, 6))
		# create map
		file1.plot(column='Phase', ax=ax, legend=True)
		# add a title
		ax.set_title(file1N, size=25)
		# create an annotation for the data source
		ax.annotate('Source: ftp://swc:water@lidarftp.swc.nd.gov', 
			xy = (.35, .135),
			xycoords='figure fraction', fontsize=7, verticalalignment='top',
			horizontalalignment='right', transform=ax.transAxes, color='gray')
		# set xy labels
		ax.set_xlabel('Meters (epsg:3857)')
		ax.set_ylabel('Meters (epsg:3857)')
		fig.savefig(f1out+file1N+'.png', dpi=300)
		# plt.show()

def Plotfile2(Plotf2):
	if Plotf2=='True':
		# Plot file2
		## Will show the shapefile extent plain without basemap
		file2N = 'James River Ph5 QL3'
		# create figure and axes for Matplotlib
		fig, ax = plt.subplots(1, figsize=(10, 6))
		# create map
		file2.plot(column='Phase', ax=ax, legend=True)
		# add a title
		ax.set_title(file2N, size=25)
		# create an annotation for the data source
		ax.annotate('Source: ftp://swc:water@lidarftp.swc.nd.gov', 
			xy = (.35, .135),
			xycoords='figure fraction', fontsize=7, verticalalignment='top',
			horizontalalignment='right', transform=ax.transAxes, color='gray')
		# set xy labels
		ax.set_xlabel('Meters (epsg:3857)')
		ax.set_ylabel('Meters (epsg:3857)')
		fig.savefig(f2out+file2N+'.png', dpi=300)
		# plt.show()
		# sys.exit()

def PlotGDF(Plot_g):
	if Plot_g=='True':
		## Will show the merge of the two
		gdfN = 'James River Phase Overlay'
		# create figure and axes for Matplotlib
		fig, ax = plt.subplots(1, figsize=(10, 6))
		# create map
		gdf.plot(column='Phase', ax=ax, legend=True)
		# add a title
		ax.set_title(gdfN, size=25)
		# create an annotation for the data source
		ax.annotate('Source: ftp://swc:water@lidarftp.swc.nd.gov', 
			xy = (.35, .135),
			xycoords='figure fraction', fontsize=7, verticalalignment='top',
			horizontalalignment='right', transform=ax.transAxes, color='gray')
		# set xy labels
		ax.set_xlabel('Meters (epsg:3857)')
		ax.set_ylabel('Meters (epsg:3857)')
		fig.savefig(Out_Dir+gdfN+'.png', dpi=300)

# Call upon given defs
Plotfile1(Plotf1)
Plotfile2(Plotf2)
PlotGDF(Plot_g)

# # Plot enable section for basemap mapping
Plotf1B = 'True'
# Plotf1B = 'T'
Plotf2B = 'True'
# Plotf2B = 'T'
Plot_gB = 'True'
# Plot_gB = 'T'
# ##################

# This def will bring in a simple basemap
def add_basemap(ax, zoom, url='http://tile.stamen.com/terrain/tileZ/tileX/tileY.png'):
    # ctx is the contextily module
    xmin, xmax, ymin, ymax = ax.axis()
    basemap, extent = ctx.bounds2img(xmin, ymin, xmax, ymax, zoom=zoom, url=url)
    ax.imshow(basemap, extent=extent, interpolation='bilinear')
    # restore original x/y limits
    ax.axis((xmin, xmax, ymin, ymax))

def Plotfile1B(Plotf1B):
	if Plotf1B=='True':
		Plotf1BN = 'James River Ph4 QL3_BM'
		# create figure and axes for Matplotlib
		fig, ax = plt.subplots(figsize=(10, 10))
		# create map
		file1.plot(column='Phase', ax=ax, legend=True)
		# Provide basemap
		add_basemap(ax, zoom=10)
		ax.set_title(Plotf1BN, size=25)
		ax.set_xlabel('Meters (epsg:3857)')
		ax.set_ylabel('Meters (epsg:3857)')
		fig.savefig(f1out+Plotf1BN+'.png', dpi=300)

def Plotfile2B(Plotf2B):
	if Plotf2B=='True':
		Plotf2BN = 'James River Ph5 QL3_BM'
		# create figure and axes for Matplotlib
		fig, ax = plt.subplots(figsize=(10, 10))
		# create map
		file2.plot(column='Phase', ax=ax, legend=True)
		# Provide basemap
		add_basemap(ax, zoom=10)
		ax.set_title(Plotf2BN, size=25)
		ax.set_xlabel('Meters (epsg:3857)')
		ax.set_ylabel('Meters (epsg:3857)')
		fig.savefig(f2out+Plotf2BN+'.png', dpi=300)

def PlotGDFB(Plot_gB):
	if Plot_gB=='True':
		gdfBN = 'James River Phase Overlay_BM'
		# create figure and axes for Matplotlib
		fig, ax = plt.subplots(figsize=(10, 10))
		# create map
		gdf.plot(column='Phase', ax=ax, legend=True)
		# Provide basemap
		add_basemap(ax, zoom=10)
		ax.set_title(gdfBN, size=25)
		ax.set_xlabel('Meters (epsg:3857)')
		ax.set_ylabel('Meters (epsg:3857)')
		fig.savefig(Out_Dir+gdfBN+'.png', dpi=300)

# Call upon given defs
Plotfile1B(Plotf1B)
Plotfile2B(Plotf2B)
PlotGDFB(Plot_gB)

###################################################################################
## Clip the base lidar segments to prep only ones over lake

file3 = gpd.read_file(SF_Dir+"LakeSak-polygon.shp")
# Match crs data
file3 = file3.to_crs(epsg=3857)

# # Plot enable section for Sak map (With Basemap)
PlotSakQ = 'True'
# PlotSakQ = 'T'

def PlotSakB(PlotSakQ):
	if PlotSakQ=='True':
		PlotSakBN = 'Lake Sak Overlay_BM'
		# create figure and axes for Matplotlib
		fig, ax = plt.subplots(figsize=(10, 10))
		# create map
		file3.plot(ax=ax, legend=True, alpha=0.5, color = 'red')
		# Provide basemap
		add_basemap(ax, zoom=10)
		ax.set_title(PlotSakBN, size=25)
		ax.set_xlabel('Meters (epsg:3857)')
		ax.set_ylabel('Meters (epsg:3857)')
		fig.savefig(Out_Dir+PlotSakBN+'.png', dpi=300)

# Call the plot
PlotSakB(PlotSakQ)

# Define total segments in entire lidar overlay
TotalSegs = len(gdf)

# Clip the data
Clip = gpd.overlay(file3, gdf, how='intersection')

# Write the Clip to a shapefile 
Clip.to_file(SF_Dir+'DEM_Grid_Cut.shp', driver = 'ESRI Shapefile')

# sys.exit()

# Define total segments in selection
SelectSegs = len(Clip)

# # Plot enable section for clip mapping
# PlotLDS = 'True'
PlotLDS = 'd'
# ##################

def PlotLDselection(PlotLDS):
	if PlotLDS=='True':
		PlotLDSN = 'Selected Lidar Grid Sections'
		# create figure and axes for Matplotlib
		fig, ax = plt.subplots(figsize=(10, 10))
		# create map
		Clip.plot(column='Phase', ax=ax, legend=True)
		# Provide basemap
		add_basemap(ax, zoom=10)
		ax.set_title(PlotLDSN, size=25)
		ax.set_xlabel('Meters (epsg:3857)')
		ax.set_ylabel('Meters (epsg:3857)')
		ax.annotate('%i of %i items in selection' % (SelectSegs, TotalSegs), 
			xy = (.29, .26),
			xycoords='figure fraction', fontsize=7, verticalalignment='top',
			horizontalalignment='right', transform=ax.transAxes, color='black')
		fig.savefig(Out_Dir+PlotLDSN+'.png', dpi=300)

# Call the plot
PlotLDselection(PlotLDS)





###################################################################################
# Prep data in clip for requesting download

# print Clip.columns.values
# Returned:
# [u'NAME' u'Project' u'Phase' u'Block' u'TileNum' u'Filename' u'Path'
#  u'QL_Level' u'Contours' u'LAS_Class' u'H_Project' u'H_Datum' u'V_Datum'
#  u'Vert_RMSEz' u'Pt_Density' u'NPS' u'DEM_PstSp' u'MinCntrInt'
#  u'Start_Date' u'End_Date' u'ASCII_Path' u'tilename' u'Comments'
#  u'SRC_SRS' u'Name' u'descriptio' u'tessellate' 'geometry']

###########################
# NAME: The numerical value for the base file name. This stays the same. Ex - 02765270
# Project: This provides the name, phase, and quad. Ex - LiDAR_JamesRiver_Ph5_QL3
# Phase: This provides the phase number. Ex - 5
# Block: Gives the letter associated with block. Ex - C
# TileNum: Gives number of tile. Ex - 02765270 (Same as NAME)
# Filename: Provides a file name for download. Contains raw data?? Ex - r02765270.zip
# Path: The path to the phase and block of the data. Ex - Phase5/Block_C-1234
# QL_Level: The quadrant of given file. Ex = 3
# Contours: IDK it just says no. Ex - No
# LAS_Class: IDK it just says yes and then a number. Ex - Yes-6
# H_Project: Projection information. Ex - UTM, Zone 14 North, m
# H_Datum: Horizontal Coordinate system info. Ex - NAD83 (NSRS2007)
# V_Datum: Vertical datum coordinate system info. Ex - NAVD88-Geoid03, Units-meters
# Vert_RMSEz: Resolution? Ex - 15cm
# Pt_Density: Ex - None (none was given for this dataset)
# NPS: Point spacing. Ex - 1.4m
# DEM_PstSp:  Ex: 1m
# MinCntrInt:  Ex - 2ft
# Start_Date: Start date. Ex - 2015/11/15 00:00:00 or 11/5/2013
# End_Date
# ASCII_Path: File path. Ex - LiDAR_JamesRiver_Ph5_QL3/Phase5/Block_C-1234/ASCII_Grids_ft/gft02805256.asc
## ASC file: Text file containing geometric shapes and points information; 
## exported by various Autodesk programs; saves data in a text format, 
## which is more compatible with other programs.
# tilename: Provides an ars file for the given tile. Ex - r02865268.ras (Raster image)
# Comments: None 
# SRC_SRS: Coordinate reference code. Ex - epsg:3721 
# Name: Name of location. Ex - LakeSak
# descriptio: (Incomplete spelling). None
# tessellate: (cover (a plane surface) by repeated use of a 
# single shape, without gaps or overlapping) Ex - 1
# geometry: The coordinates


# Start lists of data for DEM download and prep
# print df[['Path', 'LAS_Class', 'ASCII_Path']]
# print df['TileNum']
# Data = df['ASCII_Path'].values.tolist()
# print Data[0]
# print df['NAME']
# Names = df['Filename'].values.tolist()
# print Names

##############################################################################
# Write file with max xy coordinates for next section

# Define bounds of each poly
Bounds = pd.DataFrame(Clip['geometry'].bounds)
# Define names
NamePoly = pd.DataFrame(Clip['NAME'])
# Join data
BoundsData = Bounds.join(NamePoly, how='outer')
# Write as csv
BoundsData.to_csv(Out_Dir+'LidarPolyBounds.csv')

###############################################################################

# Data Download section

df = pd.DataFrame(Clip)

# List of file numbers
Names = df['NAME'].values.tolist()

# List of project names
Project = df['Project'].values.tolist()

# List of phase and block locations
Path = df['Path'].values.tolist()

# File parameters
Fstart = 'dem'
Fend = '.zip'
# # Base website link:
BaseSite = 'ftp://swc:water@lidarftp.swc.nd.gov/'

# Loop download
ProgBarLimit = len(Names)
print '\nDownloading Zipped .img files...\n'
for i in tqdm.tqdm(range(ProgBarLimit)):
	# Define link for given file in loop
	Link = BaseSite+Project[i]+'/'+Path[i]+'/'+'DEM_Imagine/'+Fstart+Names[i]+Fend
	# Define name of out file for file (Zip Form)
	OutFname = DzipOut+Fstart+Names[i]+Fend
	# Download file into the ZipImg folder
	urllib.urlretrieve(Link, OutFname)
	# Clear ftp cache from last request to send new one (python 2.7 urllib issue)
	urllib.urlcleanup()

# Unzip Files in folder
# Unzip files
temp = list(os.listdir(DzipOut))
Loc = 'DEM_Output/ZipImg/'

ProgBarLimit = len(temp)
print '\nExtracting Zipped .img files...\n'
i = 0
for i in tqdm.tqdm(range(ProgBarLimit)):
	fname = Loc+temp[i]
	with zipfile.ZipFile(fname, 'r') as zip_ref:
	    zip_ref.extractall(Lidar_Out)
	i = i+1
