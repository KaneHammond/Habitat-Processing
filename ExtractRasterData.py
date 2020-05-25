# Module purpose is to extract NDVI calculations from nesting point data
# as well as distance from shoreline to water (approximately)

import Modules
from Modules import *

# Identify basic file paths
TCI = 'TIF_Code_Images/'
PLD = 'ParsedLakeData/'
SHP_Loc = 'Shapefile/'
SegMeta = 'Analysis_Files/Analysis_Meta.txt'

TFolders = listdir(TCI)
SegmentFolders = listdir(PLD)

# Define a sorting def that will sort a given list by first value
def SortList(L):
	return(sorted(L, key=lambda x: x[0]))

################################ Let the user know it has been initiated ##########

# MUTE HERE 2

print '\nPreparing Image File Paths...\n'

################################# Open nest training data file #################
# If not written yet, it will be written here. This will be used to determine if
# new data is availble.

# CheckData = []

# with open('Analysis_Files/Distance.csv', 'r+') as File:
# 	for aRow in File:
# 		CheckData.append(aRow)

# print CheckData
# sys.exit()



# ADD IT HERE WORK WORK WORK
# You can make the list open/write and then compare said list to the nests on file
# if it is a complete match, it should be noted that this section is not required.
# If nests are identified, those nests will have to be taken out separately for 
# the training module and appended to existing training data. This may be best to be
# placed with the segmented nest information, as the addition of new segments to an 
# overall analysis will not be noticed at this point. Skip for now.
########################## CHECK COORDINATE SYSTEMS ########################

# If the CRS does not match between the shapefiles and the raster files,
# it is much simpler to modify the shapefiles crs. We will do this if needed.

# List all raster files in folders
BasePaths = list(os.listdir('TIF_Code_Images/'))

with rasterio.Env():
	for aPath in BasePaths:
		Files = list(os.listdir('TIF_Code_Images/'+aPath))
		for File in Files:
			FP = 'TIF_Code_Images/'+aPath+'/'+File
			with rasterio.open(FP, 'r') as src:
				# We only need to request the crs once
				src_crs = src.crs
			src.close()
			break
		break

# We now know the standard crs for the tif images through src_crs

################################## Image Band PREP ################################

# Image file locations

Image_Folders = []
for Folder in TFolders:
	# List folders within the TIF folder location
	Temp = listdir(TCI+Folder)
	# Append the folder and image files within that folder to the list
	Image_Folders.append([Folder, Temp])

# Image_Folders FORMAT: [Folder, [Files in Folder (image bands)]]

# Open folders to select extent files. THese will be used to relate imagery 
# to segments. Then that data will be referenced to the nests within them.

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
		# print SubFolder
		# sys.exit()
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
				# MODIFIED TO USE SEG FILE NO LONGER EXTENT
				if str(BaseName[0:3])=='Seg':
					if str(FileEnding)=='shp':
						FileLocation = PLD+Folder+'/'+SubFolder+'/'+FName
						# Open the shapefile
						crop_extent = gpd.read_file(FileLocation)
						# Save it to raster crs
						crop_extent = crop_extent.to_crs(src_crs)
						# Save the crop extent as the segment extent
						ExLoc = PLD+Folder+'/'+SubFolder+'/'+'Extent'+'.shp'
						Save = crop_extent.to_file(ExLoc)
						# Add the file path to the extent location
						ExtentLoc.append(ExLoc)


################################# NEST INFORMATION PREP #######################
# Convert EPSG

# NestDat = gpd.read_file(SHP_Loc+'Nests.shp')
NestDat = gpd.read_file(SHP_Loc+'Nests.shp')
# Save it to raster crs
NestMod = NestDat.to_crs(src_crs)
# Save the crop extent as the segment extent
Save = NestMod.to_file(SHP_Loc+'NestsConverted.shp')

# Iterate through the points

points = ([pt for pt in fiona.open(SHP_Loc+'NestsConverted.shp')])

# # List of Nest location and associated segment
# NSR = []
# print 'NSR...'
# i = 0
# ProgBarLimit = len(ExtentLoc)

# for z in tqdm.tqdm(range(ProgBarLimit)):
# # for Segment in ExtentLoc:
# 	# Open the segment to read it
# 	with fiona.open(ExtentLoc[z], 'r') as FileUse:
# 		for feature in FileUse:
# 			# Define the extent geometry and Id
# 			SegId = str(feature['properties']['Id'])
# 			Geom = feature['geometry']
# 		for i, pt in enumerate(points):
# 			# Iterate through points to identify which ones are
# 			# located within the segment.
# 			point = shape(pt['geometry'])
# 			if point.within(shape(Geom)):
# 				NSR.append([str(pt['properties']['Nest Id']), SegId])

# 	FileUse.close()
# 	# print Segment

# To check if needed
# df = pd.DataFrame(NSR)
# df.to_csv('NSR.csv')
# df=None

# ExtentLoc no longer needed, close it out
ExtentLoc = None
t = []
with open('NSR.csv', 'r') as File:
	for M in File:
		M = M.replace('\n', '')
		Sp1 = M.split(',')
		t.append([Sp1[1], Sp1[-1]])
NSR = None
t = t[1::]
NSR =[]
for x in t:
	NSR.append(x)

# sys.exit()

# Write a data set with date information for the nests
NestDates = []
for i, pt in enumerate(points):
	# Calculated hatch date
	EstDate = str(pt['properties']['Calc Hatch'])
	try:
		# Attempt to backtrack the lay date by its incubation.
		# This estimate will allow for the selection of the most
		# appropriate image file by date comparison.
		ConvDate = datetime.datetime.strptime(EstDate, '%m/%d/%Y')
		# Incubation period subtracted to estimate actual nest date
		PredNestDate = ConvDate-timedelta(30)
		# Write to list with nest id, add P for Predicted Nesting Date
		NestDates.append([str(pt['properties']['Nest Id']), PredNestDate, 'P'])
	except:
		# If fails, use the date the nest was found..
		DateFnd = pt['properties']['Init Date']
		DateFnd = datetime.datetime.strptime(DateFnd, '%m/%d/%Y')
		# Add F for Found Date
		NestDates.append([str(pt['properties']['Nest Id']), DateFnd, 'F'])

# Add the nest location (segment number) to the dated set
# FORMAT = [Nest Id, Seg, Date, F/P]
CompNest = []
for aNest in NestDates:
	# For each dated nest, compare it to the segment/nest dataset
	for SegNest in NSR:
		if SegNest[0]==aNest[0]:
			CompNest.append([SegNest[0], SegNest[1], aNest[1], aNest[-1]])

# Remove points from memory
points = None
# Nest Dates not needed 
NestDates = None

# Filter out the points we cannot use (Date Found with no float)

FilteredNests = []

for aItem in CompNest:
	if aItem[-1]=='P':
		FilteredNests.append(aItem)

# CompNest not needed anymore
CompNest = None
# NestDat not needed
NestDat = None
NestProj = None

################################# Segment information Prep ############################

# Determine which image files are needed for each nest, filter by path row first
# via segment location info
Raw = []
with open(SegMeta, 'r') as File:
	for aItem in File:
		Raw.append(aItem)
# Close the file
File.close()

# Split the data
Raw = Raw[0].split('-')
# Split the path row data (second part of list)
Raw = Raw[1].split(',')
# Write new list of corrected data. Format is currently PR:'segment'.'path'/'row'

SegmentMeta = []
for aItem in Raw:
	# Split the list item at ':'
	Sp1 = aItem.split(':')
	# Split this again at '.'
	Sp2 = Sp1[1].split('.')
	# Split the p/r
	# Sp3 = Sp2[1].split('/')
	# print Sp2
	# Define the pr code and segment for list
	# PRCode = str(Sp3[0]+Sp3[1])
	PRCode = str(Sp2[1])
	# Append items to list
	SegmentMeta.append([Sp2[0], PRCode])

# SegmentMeta Format: [Segment, pr code]

# Remove lists/items from memory
Sp1 = None
Sp2 = None
Sp3 = None
PRCode = None
Raw = None

##################################### Combine Segment and Image information #######################

# Define image paths for each segment

SegImagePaths = []

for Segment in SegmentMeta:
	# SegmentMeta Format: [Segment, pr code]
	TempImageFolders = []
	for aItem in Image_Folders:
		# Format of Image_Folders: [Image Folder, [Images]]
		# Split the image info
		Sp1 = aItem[0].split('_')
		# Define the pr of the image folder
		FolderPR = Sp1[1]
		# If they match for the seg, save info to list
		if Segment[1]==FolderPR:
			# Define full path using the image folder, will include all dates
			Path = 'TIF_Code_Images/'+aItem[0]+'/'
			TempImageFolders.append(Path)
	SegImagePaths.append([Segment[0], Segment[1], TempImageFolders])


# To check data
df = pd.DataFrame(SegImagePaths)
df.to_csv('SegImagePaths.csv')
df = None

# Close/erase items no longer needed
Sp1 = None
FolderPR = None
Path = None
TempImageFolders = None
SegmentMeta = None
Image_Folders = None

#################################### Write list for Nest/Image Paths #########################

print '\nCombining Image Data With Nest Data...\n'

NestImInfo = []

for Nest in FilteredNests:
	# Define individual variables for filtering and keeping
	NestDate = Nest[2]
	NestID = Nest[0]
	NestSegment = Nest[1]
	# TempImages will be the best fitted images based upon nesting month.
	TempImages = []
	for Paths in SegImagePaths:
		# Filter through image paths by segment and then date
		Segment = Paths[0]
		if Segment == NestSegment:
			# Define the image dates and compare each to the nest date
			ImagePaths = Paths[2]
			for Image in ImagePaths:
				# Extract image date
				Sp1 = Image.split('/')
				Sp2 = Sp1[1].split('_')
				DateCode = Sp2[-1]
				Date = datetime.datetime.strptime(DateCode, '%Y%m%d')
				# Use the date info for matching appropriate image to nest
				if Date.year == NestDate.year:
					if Date.month == NestDate.month:
						# Select the image(s) and determine the closest one to the
						# predicted nesting date.
						TempImages.append([Image, Date])

	# Calculate the best fit image by date
	Days = []
	# Days are the absolute value of days away from predicted nesting date
	for anImage in TempImages:
		DayOfIm = anImage[1].day
		Calc = abs(int(DayOfIm-NestDate.day))
		Days.append([Calc, anImage])

	# Sort the list by first value
	Days = SortList(Days)
	# Select most appropriate date
	try:
		SelectedImagePath = Days[0][1][0]
		# Make the list needed to process the nests. Append the nest id and the 
		# appropriate image folder path.
		NestImInfo.append([NestID, SelectedImagePath])
	except:
		# Assume a date error, i.e. the images for that month are not downloaded
		# print NestDate
		pass
Sp1 = None
Sp2 = None
Days = None
Calc = None
DayOfIm = None
TempImages = None
Date = None
NestDate = None
NestID = None
NestSegment = None
DateCode = None
SegImagePaths = None


################################# Write this meta data to lists ######################
# Useful for working on the program. Larger data sets take a long time to process this

# At this point, the only list needed to continue the rest of the analysis is
# NestImInfo. Write to list, mute rest of code for testing other sections!

df = pd.DataFrame(NestImInfo)
df.to_csv('Analysis_Files/NestImInfo.csv')
df = []
# sys.exit()

# MUTE HERE 2

# NestImInfo FORMAT: [NestID, Image folder path]

NestImInfo = []

with open('Analysis_Files/NestImInfo.csv', 'r') as File:
	for aItem in File:
		aItem = aItem.replace('\n', '')
		Sp1 = aItem.split(',')
		# Sp1[0] is the index from the csv. We do not neet that.
		NestImInfo.append([Sp1[1], Sp1[2]])
File.close()

#Drop first rec
NestImInfo = NestImInfo[1::]

################################### Filter out nests near interference ##################

def get_value_at_point(rasterfile, pos):
	gdata = gdal.Open(rasterfile)
	gt = gdata.GetGeoTransform()
	data = gdata.ReadAsArray().astype(np.float)
	gdata = None
	y = int((pos[0] - gt[0])/gt[1])
	x = int((pos[1] - gt[3])/gt[5])
	return data[y, x]

# Mute 1

# Segment max shoreline distances
DistMax = []

ProgBarLimit = len(NestImInfo)

print '\nHabitable Zone Calculations...\n'
for z in tqdm.tqdm(range(ProgBarLimit)):
# for aItem in NestPixQua:
	aItem = NestImInfo[z]
	# Select single nest point of interest by "shortList"
	pts = gpd.read_file(SHP_Loc+'NestsConverted.shp')
	shortlist = [aItem[0]]

	# Convert the Id variables to string to match the list used
	# for filtering nests within a given segment.
	pts['Nest Id'] = pts['Nest Id'].astype('str')

	# Mask only the points within this year. This removes all 
	# other nests and focuses only on those within the correct
	# date bounds.
	mask = ~pts['Nest Id'].isin(shortlist)
	pts = pts.loc[~mask]

	# Write columns containing x and y coordinate data.
	# This specific method of data extraction requires this
	# for looping the data.
	pts['x'] = pts['geometry'].x
	pts['y'] = pts['geometry'].y
	pts = pts[['x', 'y', 'Nest Id', 'geometry']]
	coords = [(x,y) for x, y in zip(pts.x, pts.y)]

	# Rename header Nest Id to work with loop, need to remove space
	pts = pts.rename(columns = {'Nest Id': 'NEST_ID'})


	################ Initially used a shapefile here. 

	# Define the pixel_qu tiff file
	Sp1 = aItem[-1].split('/')
	PixelQuFile = str(aItem[-1]+Sp1[1]+'_'+'pixel'+'_'+'qa.tif')

	################################### Check exact pixel value (pixel quality)

	# Open the filtered NDVI clip file
	src = rasterio.open(PixelQuFile)

	# View a the raster below to check
	# pyplot.imshow(src.read(1), cmap = 'jet')
	# plt.show()
	# sys.exit()

	# Loop through the selected nests and extract the NDVI 
	# values at those locations

	# pts['Raster Value'] = [x[0] for x in src.sample(coords)]
	pts['Raster Value'] = [x[0] for x in src.sample(coords)]
	################################### Check exact pixel value (pixel quality)
	# Convert to list
	ptslist = pts['Raster Value'].values.tolist()

	# Write list of accepted variables that may be found in raster data.
	# Zero typically means no data, while this may cause error in some cases
	L_W = [322, 324, 0]

	if ptslist[0]==322:
		# Only run if point is identified on land
		Continue = 0
		TempDist = []
		# Incr must be at least 1 to write file
		Incr = 1
		i = 0
		Count = 0
		while Continue==0:
			for index, row in pts.iterrows():
				ID = int(row.NEST_ID)
				PointBuff = row.geometry.buffer(Incr)
				# PixQua = aItem[-1]+

				# This is a shapely item, write it to a shapefile using ogr
				src_ds = gdal.Open(PixelQuFile)
				# Copy the spatial reference system of the file used to create the shp file
				srs = osr.SpatialReference()
				srs.ImportFromWkt(src_ds.GetProjectionRef())

				#  create output datasource
				dst_layername = 'CVar'
				drv = ogr.GetDriverByName("ESRI Shapefile")
				dst_ds = drv.CreateDataSource("Analysis_Files/WorkingFile.shp" )
				dst_layer = dst_ds.CreateLayer(dst_layername, srs = srs)
				defn = dst_layer.GetLayerDefn()
				feat = ogr.Feature(defn)
				geom = ogr.CreateGeometryFromWkb(PointBuff.wkb)
				feat.SetGeometry(geom)
				dst_layer.CreateFeature(feat)

				feat = geom = None
				src_ds = None
				srcband = None
				srs = None
				drv = None
				dst_ds = None
				dst_layer = None
				field_defn = None
				dst_layer = None
				# org complete
		
				# Use fiona to open the shapefile and clip the tif for calcualtions
				with fiona.open('Analysis_Files/WorkingFile.shp', "r") as shapefile:
					shapes = [feature["geometry"] for feature in shapefile]
				shapefile.close()
				import rasterio.mask
				# Add the pixel quality tiff here
				# with rasterio.open(aItem[-1]) as src:
				with rasterio.open(PixelQuFile) as src:
					out_image, out_transform1 = rasterio.mask.mask(src, shapes, crop=True)
					out_meta = src.meta
					# Define the parameters of the tif clip
					out_meta.update({"driver": "GTiff",
						 "height": out_image.shape[1],
						 "width": out_image.shape[2],
						 "transform": out_transform1})
				# Write a clipped tif file
				with rasterio.open('Analysis_Files/WorkingClipFile.tif', 'w', **out_meta) as dst:
					dst.write(out_image)
		
				dst.close()
				src.close()

				# Open the raster file and calculate raster values
				src = rasterio.open('Analysis_Files/WorkingClipFile.tif')
				array = src.read(1)

				# The array seems to include 0 into the general assessment of the area
				# possibly assigning 0 to the outer bounds that do not have data. As
				# arrays are not circular. Might have to allow zeros to pass.
				ArrayList = array.tolist()

				# Error check
				E = 0

				for AS in ArrayList:
					for V in AS:
						if V != 322:
							if V != 324:
								# In this data, 1 is placed to fill the missing
								# array values. Changed from 0 from the Google
								# imagery data.
								if V!=1:
									# If a value is outside of these three,
									# it is interference.
									TempDist.append('Fail')
									E = 1
									Continue = 1

				if E == 0:
					for AS in ArrayList:
						for V in AS:
							if V == 324:
								TempDist.append(Incr)
								Continue = 1
								# print array
								# sys.exit()
						
				i = i+1
				Incr = 15*i
				# print Incr
				array = None
				src = None

		
		DistMax.append([aItem[0], TempDist[0]])
		# print DistMax
		# Delete the working files
		os.remove('Analysis_Files/WorkingFile.shp')
		os.remove('Analysis_Files/WorkingClipFile.tif')


DistDf = pd.DataFrame(DistMax, columns = ['Nest Id', 'Distance'])

DistDf.to_csv('Analysis_Files/NestDistance.csv')

# Close out of DistDf. This section is mostly practical for working
# on sections of code and being able to run these portions separate.

DistDf = None

# sys.exit()

# Mute 1

######################### Train the data for NDVI ############################
# Use the Nest Distance data for this section. This data is already filtered by
# location and potential nest point location in a given satellite image i.e. in 
# a pixel with interference or defined as water instead of land.

################################ Filter through sucessful records #################

RawDist = []

with open('Analysis_Files/NestDistance.csv', 'r') as File:
	for aRow in File:
		RawDist.append(aRow)
File.close()

# Clean up the data. Do not include nests which failed the initial anlaysis
# due to interference.

NestsNeeded = []

for aRow in RawDist:
	aRow = aRow.replace('\n', '')
	aRow = aRow.split(',')
	if aRow[2]!='Fail':
		NestsNeeded.append(aRow[1])

NestsNeeded = NestsNeeded[1::]
RawDist = None

###################### Relate sucessful records to image paths ####################

NestIm2 = []


for aRow in NestImInfo:
	for aItem in NestsNeeded:
		if aRow[0]==aItem:
			NestIm2.append(aRow)

NestsNeeded = None
# NestIm2 FORMAT: ['Nest Id', 'Base Image Path (not band specific)']

# Add the actual image paths for NDVI to the list.

for aRow in NestIm2:
	FP = aRow[1]
	Sp1 = aRow[1].split('/')
	NDVIPath = FP + Sp1[1] + '_ndvi.tif'
	# B5Path = FP + Sp1[1] + '_B5.tif'
	aRow.append(NDVIPath)
	# aRow.append(B5Path)

# Clear remaining lists to free up memory
NestImInfo = None

##################### Calculate NDVI at nest points #############################
print '\nCalculating NDVI From Nesting Pixels...\n'

NestNDVI = []

ProgBarLimit = len(NestIm2)

for z in tqdm.tqdm(range(ProgBarLimit)):
	aItem = NestIm2[z]

	# Open the point shapefile 

	# Select single nest point of interest by "shortList"
	pts = gpd.read_file(SHP_Loc+'NestsConverted.shp')
	shortlist = [aItem[0]]

	# Convert the Id variables to string to match the list used
	# for filtering nests within a given segment.
	pts['Nest Id'] = pts['Nest Id'].astype('str')

	# Mask only the points within this year. This removes all 
	# other nests and focuses only on those within the correct
	# date bounds.
	mask = ~pts['Nest Id'].isin(shortlist)
	pts = pts.loc[~mask]

	# Write columns containing x and y coordinate data.
	# This specific method of data extraction requires this
	# for looping the data.
	pts['x'] = pts['geometry'].x
	pts['y'] = pts['geometry'].y
	pts = pts[['x', 'y', 'Nest Id', 'geometry']]
	coords = [(x,y) for x, y in zip(pts.x, pts.y)]

	# Rename header Nest Id to work with loop, need to remove space
	pts = pts.rename(columns = {'Nest Id': 'NEST_ID'})

	################ Define the NDVI pixel value
	# print aItem
	srcNDVI = rasterio.open(aItem[-1])
	# pyplot.imshow(srcNDVI.read(1), cmap='jet')
	# pyplot.show()
	# sys.exit()

	# Loop through the selected nests and extract the NDVI 
	# values at those locations

	pts['Raster Value'] = [x[0] for x in srcNDVI.sample(coords)]
	# print pts

	# Convert to list
	ptslist = pts['Raster Value'].values.tolist()

	# Scale Factor for this data needs to be applied
	NDVI = float(ptslist[0])*0.0001
	# print NDVI
	# sys.exit()

	# Close out what has been used
	srcNDVI = None
	ptslist = None

	NestNDVI.append([aItem[0], NDVI])

	pts = None
	mask = None
	shortlist = None
	coords = None
	ptslist = None
	B4val = None
	B5val = None

# Write this list to a dataframe and then to a csv
NDVIdf = pd.DataFrame(NestNDVI, columns = ['NestId', 'NDVI'])
NDVIdf.to_csv('Analysis_Files/NestNDVI.csv')

# Clear the list
NDVIdf = None
NestNDVI = None



