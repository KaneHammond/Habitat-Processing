# Module purpose is to extract NDVI calculations from nesting point data
# as well as distance from shoreline to water (approximately)

import Modules
from Modules import *

# Identify basic file paths
PLD = 'ParsedLakeData/'

SHP_Loc = 'Shapefile/'

# List folders within directory, this is the segment folder location
SegmentFolders = listdir(PLD)

def SortList(L):
	return(sorted(L, key=lambda x: x[0]))

# OPEN FOLDERS TO SELECT EXTENT FILES AND NDVI FILES

# Open the extent files specific to this analysis. This is done by
# using folders made in parsed lake data. 
# Paths to extent files
ExtentLoc = []
# Paths to NDVI folders
NDVI_Folders = []
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
		if FName[-4::]=='NDVI':
			NDVI_Loc = PLD+Folder+'/'+str(SubFolder)
			NDVI_Folders.append(NDVI_Loc)

# Define crs and write a geo loc based item 
# MainCrs = fiona.open(SHP_Loc+'DEM_Grid_Cut.shp', 'r').crs

# Nests = fiona.open(SHP_Loc+'Nests.shp', 'r')
# Nests.crs = 

# for aItem in NDVI_Folders:
# 	print aItem
# 	sys.exit()

# if Nests.crs != MainCrs:
# 	print 'Incorrect crs'

################################# NEST INFORMATION PREP #######################
# Convert EPSG

# NestDat = gpd.read_file(SHP_Loc+'Nests.shp')
NestDat = gpd.read_file(SHP_Loc+'Nests.shp')
NestProj = NestDat.copy()

NestProj['geometry']=NestProj['geometry'].to_crs(epsg=32614)

NestProj.crs = from_epsg(32614)

NestProj.to_file(SHP_Loc+'NestsConverted.shp')

# ITERATE THROUGH THE POINTS

points = ([pt for pt in fiona.open(SHP_Loc+'NestsConverted.shp')])

# List of Nest location and associated segment
NSR = []

for Segment in ExtentLoc:
	# print Segment
	with fiona.open(Segment, 'r') as FileUse:
		for feature in FileUse:
			# print str(feature['properties']['Id'])
			# Standard id from feature is set to zero for all of the polygons
			# in this segment file 4/16/20. Keep an eye on this when switching
			# between shapefiles.
			Segment = str(feature['properties']['Id'])
			Geom = feature['geometry']

		for i, pt in enumerate(points):
			# print pt['properties']['Id']
			# sys.exit()
			point = shape(pt['geometry'])
			if point.within(shape(Geom)):
				NSR.append([str(pt['properties']['Nest Id']), Segment])
	# sys.exit()
	FileUse.close()

# Sort the values by segments, leaving only single records per seg
Temp = []
for aItem in NSR:
	Temp.append(aItem[-1])
# Write unique values
TempMod = list(set(Temp))

# Write dataset with segment info with all nests in a given segment

NSR_Filtered = []

i = 0
while i<len(TempMod):
	NT = []
	for aItem in NSR:
		if aItem[-1]==TempMod[i]:
			NT.append(aItem[0])
	NT.insert(0, TempMod[i])
	NSR_Filtered.append(NT)
	i = i+1

# Write a data set with date information for the nests
NestDates = []
for i, pt in enumerate(points):
	# Calculated hatch date
	EstDate = str(pt['properties']['Calc Hatch'])
	try:
		ConvDate = datetime.datetime.strptime(EstDate, '%m/%d/%Y')
		# Incubation period subtracted to estimate actual nest date
		PredNestDate = ConvDate-timedelta(30)
		# Write to list with nest id, add P for Predicted Nest Date
		NestDates.append([str(pt['properties']['Nest Id']), PredNestDate, 'P'])
	except:
		# If fails, use the date the nest was found
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

##################################### NDVI DATA SELECTION ########################

# Set a list of years available from data.
YRaw = []
for aRow in NestDates:
	D = aRow[1]
	Y = D.year
	YRaw.append(Y)
# Filter for unique values
Temp = Count = collections.Counter(YRaw)
Years = Temp.keys()
# Sort the years found in the nest data. Then use these years for parsing the
# NDVI extraction section.
Years.sort()

# Open the coverage quality meta file. This will contain information
# explaining the intereference within a given segment.
MetaCovRaw = []

with open('Analysis_Files/SegmentImageInfo.csv', 'r') as File:
	for aRow in File:
		MetaCovRaw.append(aRow)
File.close()

# Drop the header
MetaCovRaw = MetaCovRaw[1::]

# Clean up the data
MetaCov = []
for aItem in MetaCovRaw:
	# Split the string item
	L_item = aItem.split(',')
	# Drop the index value
	L_item = L_item[1::]
	# Pull the date code
	DC = L_item[-1].split('_')
	DC = DC[-1]
	DC = int(DC)
	# Append the clean data to the MetaCov list element
	MetaCov.append([L_item[0], L_item[1], L_item[2], L_item[3], L_item[4], DC])

# MetaCov Format: [Segment, Clear land, Water, Interference, Landsat Code, Date Code]

############################ Extract NDVI Values for Nests ##############
# # Define a total DF to input all nest info data for NDVI
# NDVI_DF = pd.DataFrame(columns=['Segment', 'Nest Id', 'NDVI', 'Date', 'Year', 'SatImage'])


# #### Start Loop here *** Process NDVI MUTE 1


# ProgBarLimit = len(Years)

# print '\nExtracting Nest NDVI Values...\n'
# for z in tqdm.tqdm(range(ProgBarLimit)):
# # for aYear in Years:
# 	# print aYear
# 	aYear = Years[z]
# 	# Start by filtering each year, then specific months
# 	TempNests = []
# 	# Filter through dated nest info
# 	for aNest in CompNest:
# 		# Check year
# 		if aNest[2].year == aYear:
# 			# Select those that had a predicted date (meaning eggs were floated)
# 			if aNest[-1] == 'P':
# 				TempNests.append(aNest)

# 	# TempNests are the nest for this specific year

# 	# Select appropriate image for the nest and date
# 	# Define the appropriate dates within a year, i.e. their predicted nesting
# 	# date. Which is determined via backtracking the dates.

# 	# Define months and possible images to train the data.
# 	# Select from metadata file, then filter by nest date.
# 	# MetaCov Format: [Segment, Clear land, Water, Interference, Landsat Code, Date Code]
# 	ImageInfo = []
# 	for aItem in MetaCov:
# 		# Convert the date code to the datetime format, could have done it
# 		# in a shorter way...
# 		DCT = str(aItem[-1])
# 		DC = DCT[0:4]+'/'+DCT[4:6]+'/'+DCT[6:8]
# 		# DateItem is our converted datetime object for the satellite image
# 		DateItem = datetime.datetime.strptime(DC, '%Y/%m/%d')
# 		if DateItem.year==aYear:
# 			ImageInfo.append(aItem)

# 	# Determine which image was collected most closely with the associated 
# 	# nesting date (predicted)
# 	Nest_Image = []
# 	# TempNests are nests for this particular year
# 	for aNest in TempNests:
# 		# Estimated lay date
# 		EDL = aNest[2]
# 		# Write a temp list for the associated images in relation to date
# 		# for this nest
# 		TempImageDates = []
# 		NestSeg = aNest[1]
# 		# Filter by month first and segment
# 		for anImage in ImageInfo:
# 			# Check segment match
# 			# FORMAT ISSUE!!! The segments for images start at 0, the segment
# 			# files are supposed to start at 1 TEMP FIX WORK HERE
# 			# Add 1 to the segment value
# 			SegmentMod = str(int(anImage[0])+1)
# 			if SegmentMod==NestSeg:
# 				# Convert datecode back to date item
# 				DCT = str(anImage[-1])
# 				DC = str(DCT[0:4])+'/'+str(DCT[4:6])+'/'+str(DCT[6:8])
# 				# DateItem is our converted datetime object for the satellite image
# 				DateItem = datetime.datetime.strptime(DC, '%Y/%m/%d')
# 				# Check month match
# 				if EDL.year == DateItem.year:
# 					if EDL.month == DateItem.month:
# 						TempImageDates.append(anImage)

# 		# Filter the images and select the most closely related one
# 		TIDMOD = []
# 		for anImage in TempImageDates:
# 				# Convert datecode back to date item
# 				DCT = str(anImage[-1])
# 				DC = DCT[0:4]+'/'+DCT[4:6]+'/'+DCT[6:8]
# 				# DateItem is our converted datetime object for the satellite image
# 				DateItem = datetime.datetime.strptime(DC, '%Y/%m/%d')
# 				# Determine the absolute value of the diffrence between
# 				# predicted nesting datea and image date
# 				diff = abs(int(EDL.day-DateItem.day))
# 				# Append this information to the TempImageDates data.
# 				# The smallest difference will be the date selected
# 				T = copy.deepcopy(anImage)
# 				T.append(diff)
# 				TIDMOD.append(T)

# 		TIDMOD = SortList(TIDMOD)
# 		# The best image will be the last in the list
# 		BestFitImage = TIDMOD[-1]
# 		# The best fit image is the image file for that specific nest
# 		# FORMAT: [Nest Id, Image]
# 		Nest_Image.append([aNest[0], BestFitImage[4]])

# 	# Add the appropriate image information to the TempNests data
# 	# Define the images to be used
# 	ImUse = []
# 	for aItem in Nest_Image:
# 		ImUse.append(aItem[-1])
# 		# Filter through each item in nest image
# 		for Nest in TempNests:
# 			# Select a nest record matching nest image and append the image
# 			# information to the TempNests record. This will simplify the filter
# 			if aItem[0]==Nest[0]:
# 				Nest.append(aItem[-1])

# 	# Some nests may be using the same image file, we will need to clarify this
# 	# Filter the ImUse data, keep only the "keys", no double listings
# 	ImUse = collections.Counter(ImUse).keys()
# 	# Define the processing images and nests
# 	ProcRecs = []
# 	for aItem in ImUse:
# 		Temp = []
# 		for aNest in TempNests:
# 			if aNest[-1]==aItem:
# 				Temp.append(aNest[0])
# 		ProcRecs.append([aItem, Temp])
# 	# ProcRecs FORMAT:
# 	# [Image, [Nests]]

# 	# Write a dataset of only nest ids for this year, will be used in 
# 	# the nest mask below
# 	NestIds = []
# 	for aItem in TempNests:
# 		NestIds.append(aItem[0])

# 	for aRecord in TempNests:
# 		aNest = aRecord[0]
# 		# Open the Nest file
# 		pts = gpd.read_file(SHP_Loc+'NestsConverted.shp')
# 		# Open out nest converted file, this is in the correct
# 		# crs format.
# 		pts = gpd.read_file(SHP_Loc+'NestsConverted.shp')
# 		shortlist = [aNest]
# 		# Convert the Id variables to string to match the list used
# 		# for filtering nests within a given segment, this is given
# 		# in NRS filter.
# 		pts['Nest Id'] = pts['Nest Id'].astype('str')
# 		# Mask only the points within this year. This removes all 
# 		# other nests and focuses only on those within the correct
# 		# date bounds.
# 		mask = ~pts['Nest Id'].isin(shortlist)
# 		pts = pts.loc[~mask]
# 		# pts = pts['geometry']
# 		# Write columns containing x and y coordinate data.
# 		# This specific method of data extraction requires this
# 		# for looping the data.
# 		pts['x'] = pts['geometry'].x
# 		pts['y'] = pts['geometry'].y
# 		pts = pts[['x', 'y', 'Nest Id', 'geometry']]
# 		coords = [(x,y) for x, y in zip(pts.x, pts.y)]
# 		# Open the correct NDVI file. Define the landsat code in 
# 		# the accepted scene. Have to identify proper segment and image
# 		for item in TempNests:
# 			#FORMAT: [Nest Id, Seg, Date, P/F, Image]
# 			if item[0]==aNest:
# 				NestSegment = item[1]
# 				Image = item[-1]

# 		try:
# 			# Select correct NDVI image/segment
# 			for Path in NDVI_Folders:
# 				# Split the folder path info to extract the base 
# 				# landsat image name. This will include the info
# 				# for the Sensor, Path/Row, and date needed for
# 				# matching the correct NDVI clip
# 				SplitPath = Path.split('/')
# 				# Again define the segment, this time from the NDVI Folders.
# 				# This ensures the correct NDVI folder is searched
# 				# WORK HERE TEMP FIX ADD 1 to seq
# 				SegmentP=str(int(SplitPath[1][7::])+1)
# 				if SegmentP==NestSegment:
# 					# Open the correct landsat image
# 					NDVI_Files = listdir(Path+'/'+'TIF_Files')
# 					for aFile in NDVI_Files:
# 						SplitName = aFile.split('_')
# 						BaseFile = SplitName[0]+'_'+SplitName[1]+'_'+SplitName[2]
# 						# sys.exit()
# 						# the \n is in the orriginal dataset
# 						# I simply do not have it removed here and must
# 						# add it to ensure a match is made
# 						if str(BaseFile+'\n')==str(Image):
# 							CorrectNDVIFile = Path+'/'+'TIF_Files'+'/'+aFile

# 			# Open the filtered NDVI clip file
# 			src = rasterio.open(CorrectNDVIFile)

# 			# Loop through the selected nests and extract the NDVI 
# 			# values at those locations
		
# 			pts['Raster Value'] = [x[0] for x in src.sample(coords)]
# 			# Initially failed below. The code above was as follows:
# 			# pts['Raster Value'] = [x for x in src.sample(coords)]
# 			# It was modified to allow it to pass and complete the data
# 			# extraction. Leaving this here incase it needs to be used.
# 			# pts['Rater Value'] = probes.apply(lambda x: x['Raster Value'][0], axis=1)
# 			src.close()

# 			# Format a temporary dataframe
# 			Tempdf = pd.DataFrame()
# 			# (columns=['Segment', 'Nest_ID', 'NDVI'])
# 			# pts[['x', 'y', 'Id', 'geometry']]
# 			Tempdf['Nest Id'] = pts['Nest Id']
# 			Tempdf['NDVI'] = pts['Raster Value']
# 			Tempdf['Segment'] = NestSegment
# 			Tempdf['Year'] = aRecord[2].year
# 			Tempdf['Date'] = aRecord[2]
# 			Tempdf['SatImage'] = CorrectNDVIFile
# 			# append the temporary dataframe to the complete one, append
# 			# doesn't happen in place, must define the df as the append
# 			NDVI_DF = NDVI_DF.append(Tempdf, sort=False, ignore_index=True)
# 			# print Tempdf
# 			# Tempdf = None
# 			pts = None
# 		except:
# 			# Format a temporary dataframe
# 			Tempdf = pd.DataFrame()
# 			# (columns=['Segment', 'Nest_ID', 'NDVI'])
# 			# pts[['x', 'y', 'Id', 'geometry']]
# 			Tempdf['Nest Id'] = pts['Nest Id']
# 			Tempdf['NDVI'] = 'Failed'
# 			Tempdf['Segment'] = NestSegment
# 			Tempdf['Year'] = aRecord[2].year
# 			Tempdf['Date'] = aRecord[2]
# 			try:
# 				Tempdf['SatImage'] = CorrectNDVIFile
# 			except:
# 				Tempdf['SatImage'] = 'Nan'
# 			# append the temporary dataframe to the complete one, append
# 			# doesn't happen in place, must define the df as the append
# 			NDVI_DF = NDVI_DF.append(Tempdf, sort=False, ignore_index=True)

# 	# print NDVI_DF

# NDVI_DF.to_csv('Analysis_Files/NestNDVI.csv')
# NDVI_DF = None

# # END OF NDVI EXTRACT SECTION MUTE 1 START 199ish

######################## DATA PREP HABITABLE ZONE ############################

# Open the NDVI csv and select nests which passed the NDVI calculations

RAWDAT = []

with open('Analysis_Files/NestNDVI.csv') as File:
	for aRow in File:
		RAWDAT.append(aRow)
File.close()
# NDVIDATA FORMAT: [Segment, P/F/Nan (NDVI success), 
# Date, Year, Nest Id, folder path to tif]

# Drop header
RAWDAT = RAWDAT[1::]

NDVIDATA = []
for aRow in RAWDAT:
	SpltRw = aRow.split(',')
	# Remove the new line
	Path = SpltRw[6].replace('\n', '')
	NDVIDATA.append([SpltRw[1], SpltRw[2], SpltRw[3], SpltRw[4], SpltRw[5], Path])

######################### DEFINE HABITABLE ZONE FROM WATER ##################

# Distance calulation START MUTE 2

# Modify the NDVIDATA to write the foler path to the pixel quality data for
# each nest.

# NestPixQua FORMAT: [Segment, Date, Year, Nest Id, Folder path to pixel_qa]
NestPixQua = []
for aItem in NDVIDATA:
	if aItem[2] != "":
		if aItem[2] != "Failed":
			# Split the name of the previous path to extract image code values
			SpltFldr = aItem[-1].split('/')
			SpltCd = SpltFldr[-1].split('_')
			# Use image code segments to define the pixel quality file
			# WORK TEMP MOD ON SEGMENT 
			SegMod = str(int(aItem[0])-1)
			
			PixFile1 = SpltCd[0]+'_'+SpltCd[1]+'_'+SpltCd[2]+'_'+SegMod+'pixel_qa'+'.shp'
			PixFile2 = SpltCd[0]+'_'+SpltCd[1]+'_'+SpltCd[2]+'_'+SegMod+'Clip'+'.tif'
			Path1 = 'ParsedLakeData/Segment%s/Seg%sCLOUD/SHP_Files/%s' % (SegMod, SegMod, PixFile1)
			Path2 = 'ParsedLakeData/Segment%s/Seg%sCLOUD/TIF_Files/%s' % (SegMod, SegMod, PixFile2)
			# Append this data to the new list
			NestPixQua.append([aItem[0], aItem[1], aItem[3], aItem[4], Path1, Path2])

# Segment max shoreline distances
DistMax = []

# Call the segment and nest point data for extracting NDVI values
ProgBarLimit = len(NestPixQua)

print '\nMax Distance Calculations...\n'
for z in tqdm.tqdm(range(ProgBarLimit)):
# for aItem in NestPixQua:
	aItem = NestPixQua[z]

	# Select single nest point of interest by "shortList"

	# Open out nest converted file, this is in the correct
	# crs format.
	pts = gpd.read_file(SHP_Loc+'NestsConverted.shp')
	shortlist = [aItem[1]]
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


	# Open the pixel quality shapefile and extract water only
	src = gpd.read_file(aItem[-2])
	# Land = 322
	# Water = 324

	# Write list of all interference values
	RawVar = src['Value'].values.tolist()
	VariablesInit = collections.Counter(RawVar).keys()
	L_W = [322, 324, 0]
	Variables = (list(set(VariablesInit) - set(L_W)))
	src = None
	# Start raster value loop. This will use raster data to determine
	# what is present in a buffer around the nest. If int is found
	# before water, it is rejected.

	Continue = 0
	TempDist = []
	# Incr must be at least 1 to write file
	Incr = 1
	i = 0
	while Continue==0:
		for index, row in pts.iterrows():
			ID = int(row.NEST_ID)
			PointBuff = row.geometry.buffer(Incr)

			# This is a shapely item, write it to a shapefile using ogr
			src_ds = gdal.Open(aItem[-1])
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
			with rasterio.open(aItem[-1]) as src:
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
			# # Delete the band clips previously written
			# os.remove(NDVITIF+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+str(SegCount)+'_'+Band1+'_'+'Clip.tif')
			# os.remove(NDVITIF+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+str(SegCount)+'_'+Band2+'_'+'Clip.tif')
			src = rasterio.open('Analysis_Files/WorkingClipFile.tif')
			array = src.read(1)

			# Check if numbers are in array
			VarFound = 0
			for aVar in Variables:
				if int(aVar) in array:
					VarFound = VarFound+1
					Continue = 1
					TempDist.append('Fail')

			if VarFound==0:
				if 324 in array:
					Continue = 1
					TempDist.append(Incr)
					
			i = i+1
			Incr = 15*i
			array = None
			src = None

	
	DistMax.append([aItem[1], TempDist[0]])
	# Delete the working files
	os.remove('Analysis_Files/WorkingFile.shp')
	os.remove('Analysis_Files/WorkingClipFile.tif')


	############################################################

	# ############ WORKS Below, old method Keep for possible mod
	# This uses a buffer from water to points, does not allow for
	# detection of interference. Method fails due to too many items
	# to unpack from shapefile.

	# # Rename header Nest Id to work with loop, need to remove space
	# pts = pts.rename(columns = {'Nest Id': 'NEST_ID'})


	# # # Open the pixel quality shapefile and extract water only
	# src = gpd.read_file(aItem[-1])
	# src = src.rename(columns = {'Nest Id': 'NEST_ID'})

	# srcWater = src[src.Value == float(324)]
	# # Land = 328
	# # Water = 324

	# # Reset the index of the gpd to allow for loop with m
	# srcWater.reset_index(inplace = True)

	# # This is our temp distance data set. Allows for determination
	# # of best distance fit in the case of multiple water polys
	# TempL = []
	# m = 0
	# while m < len(srcWater):
	# 	# print srcWater.geometry[0]
	# 	Incr = 0
	# 	i = 1
	# 	Continue = 0
	# 	while Continue == 0:
	# 		Continue = 1
	# 		for index, row in pts.iterrows():
	# 			Geom = row.geometry
	# 			ID = int(row.NEST_ID)
	# 			WaterBuff = srcWater.geometry[m].buffer(Incr)
	# 			# for aValue in NestIds:
	# 			# 	if ID == int(aValue):
	# 			# 		# Answer was in bool format. Changed to string
	# 			# 		# A bit messy, but this works ok!
	# 			Temp = str(WaterBuff.contains(Geom))
	# 			StringAnswer = str('False' in Temp)
	# 			if StringAnswer == 'True':
	# 				Incr = 15*i
	# 				Continue = 0
	# 		i = i+1
	# 	TempL.append(Incr)
	# 	m = m+1
	# TempL.sort()
	# src = None
	# WaterBuff = None
	# srcWater = None
	# pts = None
	# DistMax.append([aItem[3], Segment, TempL[0]])
	# print DistMax
	# sys.exit()

	# Previous method, save for future needs
	##########################################################

DistDf = pd.DataFrame(DistMax, columns = ['Nest Id', 'Distance'])

DistDf.to_csv('Analysis_Files/NestDistance.csv')

sys.exit()

# Distance calulation END MUTE 2

############################# NDVI HABITAT IDENTIFY ####################

# Identify Habitat from NDVI

# print NDVI_DF

# Extract raw NDVI values to determine the range of values
NDVI_RAW = NDVI_DF['NDVI'].tolist()
NDVI_RAW.sort()

# Define Max and Min NDVI value
# MaxNDVI = float(NDVI_RAW[-1])
# MinNDVI = float(NDVI_RAW[0])

# Manual Def
MaxNDVI = float(0.13)
MinNDVI = float(0.1)

# print MaxNDVI
# print MinNDVI

for Scene in OptimalData:
	# Filter through optimal data and match the record to the 
	# correct NRS filter item which contains nest Ids.
	SegmentI = Scene[0]

	# Open the correct NDVI file. Define the landsat code in 
	# the accepted scene.
	LandsatCode = Scene[4]
	for Path in NDVI_Folders:
		print SegmentI
		# Split the folder path info to extract the base 
		# landsat image name. This will include the info
		# for the Sensor, Path/Row, and date needed for
		# matching the correct NDVI clip
		SplitPath = Path.split('/')
		# Again define the segment, this time from the NDVI Folders.
		# This ensures the correct NDVI folder is searched
		SegmentP=SplitPath[1][7::]
		if SegmentP==Segment:
			# Open the correct landsat image
			NDVI_Files = listdir(Path+'/'+'TIF_Files')
			for aFile in NDVI_Files:
				SplitName = aFile.split('_')
				BaseFile = SplitName[0]+'_'+SplitName[1]+'_'+SplitName[2]
				# the \n is in the orriginal dataset
				# I simply do not have it removed here and must
				# add it to ensure a match is made
				if str(BaseFile+'\n')==str(LandsatCode):
					CorrectNDVIFile = Path+'/'+'TIF_Files'+'/'+aFile

	 		driver = gdal.GetDriverByName('GTiff')
			ds = gdal.Open(CorrectNDVIFile)
			band =  ds.GetRasterBand(1)
			lista = np.array(band.ReadAsArray())
			# ds = None
			# lista[np.where( lista < 200 )] = 1
			lista[np.where((MinNDVI <= lista) & (lista <= MaxNDVI)) ] = 1
			lista[np.where((MinNDVI > lista) & (lista > MaxNDVI)) ] = 2

			# create new file
			file2 = driver.Create( 'THISONE.tif', ds.RasterXSize , ds.RasterYSize , 1)
			file2.GetRasterBand(1).WriteArray(lista)

			# spatial ref system
			proj = ds.GetProjection()
			georef = ds.GetGeoTransform()
			file2.SetProjection(proj)
			file2.SetGeoTransform(georef)
			file2.FlushCache()
			ds = None
			file2 = None
			src = rasterio.open('THISONE.tif')
			pyplot.imshow(src.read(1))
			pyplot.show()
			src.close()


		# Write a clipped B1 tif file
		# with rasterio.open(NDVITIF+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+SegCount+'_'+Band1+'_'+'Clip.tif', 'w', **out_meta_B1) as dst:
		#     dst.write(out_image_B1)
	# src = rasterio.open(CorrectNDVIFile)
	# sys.exit()


# lista[np.where( lista < 200 )] = 1
# lista[np.where((200 < lista) & (lista < 400)) ] = 2
# lista[np.where((400 < lista) & (lista < 600)) ] = 3
# lista[np.where((600 < lista) & (lista < 800)) ] = 4
# lista[np.where( lista > 800 )] = 5




















