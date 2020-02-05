import Modules
from Modules import *

# Identify basic file paths
PLD = 'ParsedLakeData/'

SHP_Loc = 'Shapefile/'

# List folders within directory, this is the segment folder location
SegmentFolders = listdir(PLD)

# List sub dir for NDVI data
SubDir1 = []
for aItem in SegmentFolders:
	InsertItem = listdir(PLD+aItem)
	SubDir1.append(InsertItem[1])
	
# List the ndvi image folders within each segment 
# folder (Within the sub dirs)
NDVI_Folders = []
Index = 0
for aItem in SubDir1:
	TempDir = PLD+SegmentFolders[Index]+'/'+aItem+'/'
	FolderItems = listdir(TempDir)
	for aItem in FolderItems:
		# Define folder path to image
		PathInsert = TempDir+aItem
		# Define singel tif in folder
		FolderItem = listdir(PathInsert)
		# Combine the folder path with file name
		InsertItem = PathInsert+'/'+FolderItem[0]
		NDVI_Folders.append(InsertItem)
	Index=Index+1

# MetaData = []
# for aTif in NDVI_Folders:
# 	Temp = []
# 	with rasterio.open(aTif) as imageFile:
# 		Width = imageFile.width
# 		Height = imageFile.height
# 		Bounds = imageFile.bounds
# 		Band = imageFile.read(1)
# 		DataType = Band.dtype
# 		MinVal = np.nanmean(Band)
# 		MaxVal = np.nanmax(Band)
# 		Mean = np.nanmean(Band)
# 		plt.hist(Band[~np.isnan(Band)], bins='auto')
# 		plt.show()
# 		sys.exit()

# Pull in nesting data to overlay with segment shapefile. This will allow
# for the determination of which DEM and NDVI data to use.

NestData = []

# Open the extent files specific to this analysis. This is done by
# using folders made in parsed lake data. 

# Paths to extent files
ExtentLoc = []
for Folder in SegmentFolders:
	# List folders within the parsed data folder location
	Temp = listdir(PLD+Folder)
	# Go through the subfolders for each one, listing the folders for NDVI,
	# Cloud data, shapefiles, etc.
	for SubFolder in Temp:
		FName = str(SubFolder)
		# Select folder with Shp ending
		if FName[-3::]=='Shp':
			Loc = PLD+Folder+'/'+str(SubFolder)
			Temp2 = listdir(Loc)
			# Filter through folders and select the extent shapefile
			for ShpFile in Temp2:
				# Parse the shapefile names to determine which one is correct
				FName = str(ShpFile)
				Split = FName.split('.')
				BaseName = Split[0]
				FileEnding = Split[-1]
				# Check the parsed data to filter
				if str(BaseName[-6::])=='Extent':
					if str(FileEnding)=='shp':
						# Add the file path to the extent location
						FileLocation = PLD+Folder+'/'+SubFolder+'/'+FName
						ExtentLoc.append(FileLocation)


# Define crs and write a geo loc based item 
# MainCrs = fiona.open(SHP_Loc+'DEM_Grid_Cut.shp', 'r').crs

# Nests = fiona.open(SHP_Loc+'Nests.shp', 'r')
# Nests.crs = 

# for aItem in NDVI_Folders:
# 	print aItem
# 	sys.exit()

# if Nests.crs != MainCrs:
# 	print 'Incorrect crs'

# Convert EPSG

NestDat = gpd.read_file(SHP_Loc+'Nests.shp')
NestProj = NestDat.copy()

NestProj['geometry']=NestProj['geometry'].to_crs(epsg=32614)

NestProj.crs = from_epsg(32614)

NestProj.to_file(SHP_Loc+'NestsConverted.shp')

# sys.exit()

# for Segment in ExtentLoc:
# 	FileUse = fiona.open(Segment)
# 	Temp = []
# 	for Nest in Nests:
# 		if 

points = ([pt for pt in fiona.open(SHP_Loc+'NestsConverted.shp')])

# Nest location and segment it is in list
NSR = []

for Segment in ExtentLoc:
	with fiona.open(Segment, 'r') as FileUse:
		for feature in FileUse:
			Segment = str(feature['properties']['NAME'])
			Geom = feature['geometry']
		# print Geom
		Temp = []
		for i, pt in enumerate(points):
			point = shape(pt['geometry'])
			# SHP = shape(FileUse['geometry'])
			# print SHP
			# print point
			# print shape(Geom)
			if point.within(shape(Geom)):
				NSR.append([i, Segment])
	FileUse.close()

print NSR



