import Modules
from Modules import *

# DEM: "EPSG","26914"]

sys.path.append("DEM_Output/LiDarData")
LAS_Dir = 'DEM_Output/LiDarData/'

sys.path.append("Output")
OutData = "Output/"

dir = 'DEM_Output/DEM_Poly/'
if not os.path.exists(dir):
    os.makedirs(dir)
DEM_Poly = 'DEM_Output/DEM_Poly/'

########################################################################################
# Read all DEMs in the unzip folder. Then the DEM data will be processed to polygon
# shapefiles containing elevation contours (1 meter intervals).

# List of all dem files in the data folder
Files = []
for file in os.listdir(LAS_Dir):
	if file.endswith('.img'):
		Files.append(file)

# Define files with the elevation contour we desire 

gdal.UseExceptions()

Passed = []
ElevationCont = 567
ProgBarLimit = len(Files)

print '\nFiltering DEMs...\n'
for i in tqdm.tqdm(range(ProgBarLimit)):
	src_ds = gdal.Open(LAS_Dir+Files[i])
	prj=src_ds.GetProjection()
	# print prj
	# sys.exit()
	if src_ds is None:
		print 'Unable to open %s' % Files[i]
		sys.exit(1)
	try:
		srcband = src_ds.GetRasterBand(1)
	except RuntimeError, e:
		print 'No band %i found' % 1
		print e
		sys.exit(1)
	Max = float(srcband.GetMaximum())
	Min = float(srcband.GetMinimum())
	# if Min<=ElevationCont:
	# 	if Max>=ElevationCont:
	# 		Passed.append(Files[i])
	Passed.append(Files[i])
	src_ds = []
	
# Loof for polygon creation
ProgBarLimit = len(Passed)

print '\nWriting DEM Polygons...\n'
for i in tqdm.tqdm(range(ProgBarLimit)):

	# Define a folder name for the polygon data to be written
	# Example: 02485318_Poly where the numbers are the orriginal image name
	FolderName = Passed[i]
	# FN will be used as the file name itself
	FN = FolderName.replace('.img', '_Poly')
	FolderName = FN.replace('dem', '')

	# Write folder for the new shape data
	dir = DEM_Poly+FolderName
	if not os.path.exists(dir):
	    os.makedirs(dir)
	Out_Dir = DEM_Poly+FolderName

	gdal.UseExceptions()
	#  get raster datasource
	src_ds = gdal.Open(LAS_Dir+Passed[i])
	if src_ds is None:
		print 'Unable to open %s' % src_filename
		sys.exit(1)

	try:
		srcband = src_ds.GetRasterBand(1)
	except RuntimeError, e:
		# for example, try GetRasterBand(10)
		print 'Band ( %i ) not found' % srcband
		print e
		sys.exit(1)

	# Copy the spatial reference system of the file used to create the shp file
	srs = osr.SpatialReference()
	srs.ImportFromWkt( src_ds.GetProjectionRef() )

	#  create output datasource
	dst_layername = FN
	drv = ogr.GetDriverByName("ESRI Shapefile")
	dst_ds = drv.CreateDataSource(Out_Dir + '/' + dst_layername + ".shp" )
	dst_layer = dst_ds.CreateLayer(dst_layername, srs = srs)

	field_defn = ogr.FieldDefn("elev", ogr.OFTReal)
	dst_layer.CreateField(field_defn)

	gdal.Polygonize( srcband, None, dst_layer, 0, callback=None )
######################################################################################