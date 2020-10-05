

import sys
sys.path.append("TrainingData/PythonFiles")

import Modules
from Modules import*



######################################### EXTRACT NDVI AND DISTANCE PARAMS ################################
# Define the year of interest

PrevYear = 2019
Years = [2014, 2015, 2016, 2017, 2018, 2019, 2020]
# We are looking for the meta file in MasterData under training that
# contains a combination of all years that have been trained, excluding
# our year of interest. This will provide a threshold for this specific 
# analysis that excludes its own information.

AllNDVIMeta = 'TrainingData/MasterData/TrainData/CombinationsNDVI_Meta.txt'
AllDistanceMeta = 'TrainingData/MasterData/TrainData/CombinationsDistance_Meta.txt'
ShapeFiles = 'TrainingData/Shapefile/'
Meta_Dir = 'TrainingData/Analysis_Files'
# Define a list of years that must be searched. This list will be all years
# minus the one of interest
LO = list(set(Years)-set([PrevYear]))

# Open the Meta files and select records that match the LO list

# Open datasets
NDVI = []
with open(AllNDVIMeta, 'r') as Doc:
	for aItem in Doc:
		# Split the data
		Sp1 = aItem.split(']')
		Filter1 = []
		for Data in Sp1:
			Temp = []
			Mod1 = Data.replace(']', '')
			Mod2 = Mod1.replace('[', '')
			Mod3 = Mod2.split(',')
			# remove blank spaces
			for aValue in Mod3:
				if aValue!='':
					Temp.append(aValue)
			# check if the data is there.
			# given the data format, some blank records existed
			# therefore, the lenght of some items may be 0
			if len(Temp)>0:
				NDVI.append(Temp)
Doc.close()


Distance = []
with open(AllDistanceMeta, 'r') as Doc:
	for aItem in Doc:
		# Split the data
		Sp1 = aItem.split(']')
		Filter1 = []
		for Data in Sp1:
			Temp = []
			Mod1 = Data.replace(']', '')
			Mod2 = Mod1.replace('[', '')
			Mod3 = Mod2.split(',')
			# remove blank spaces
			for aValue in Mod3:
				if aValue!='':
					Temp.append(aValue)
			# check if the data is there.
			# given the data format, some blank records existed
			# therefore, the lenght of some items may be 0
			if len(Temp)>0:
				Distance.append(Temp)
Doc.close()


# Find Records
for aItem in Distance:
	if len(aItem[3::])==len(LO):
		Convert = []
		YearList = aItem[3::]
		for aValue in YearList:
			Convert.append(int(aValue))
		CompList = set(Convert)-set(LO)
		if len(CompList)==0:
			# We have a match because it is identical
			FullDistanceTrainRecord=aItem

for aItem in NDVI:
	if len(aItem[3::])==len(LO):
		Convert = []
		YearList = aItem[3::]
		for aValue in YearList:
			Convert.append(int(aValue))
		CompList = set(Convert)-set(LO)
		if len(CompList)==0:
			# We have a match because it is identical
			FullNDVITrainRecord=aItem

# These records contain all the needed info for predicting habitat
# print FullDistanceTrainRecord
# print FullNDVITrainRecord

# PREDICTION PARAMETERS
MaxNDVI = float(FullNDVITrainRecord[0])
MinNDVI = float(FullNDVITrainRecord[1])
# MaxDistance = float(FullDistanceTrainRecord[0])
# MOD
# Was cutting out middle of bars. This will at least help when at lake.
MaxDistance = float(FullDistanceTrainRecord[0])


########################## WRITE FOLDERS FOR PREDICTION DATA #########################

dir = 'ParsedSegmentData'
if not os.path.exists(dir):
    os.makedirs(dir)
dir = 'Analysis_Files'
if not os.path.exists(dir):
    os.makedirs(dir)

PSD = 'ParsedSegmentData/'
SF_Dir = 'TrainingData/Shapefile/'

# LOCS
LakeSak = [3, 4, 5, 6, 7, 8, 9, 10, 11]
NorthDakota = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]
SouthDakota = [25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43]
Proj = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12,13,14,15,16]
RiverProj = [12,13,14,15,16]
################### Query

IDs = []
with fiona.open(SF_Dir+'RiverFull.shp') as File:
	for feature in File:
		IDs.append(feature['properties']['Id'])

File.close()

SegmentSections = [LakeSak, NorthDakota, SouthDakota]
Options = ['Lake Sakakawea', 'North Dakota (River Segments)', 'South Dakota']

# SKIP QUERY HERE


################## MAIN SEGS

# c = 1
# print 'Available Segments:'
# for aItem in Options:
# 	print '%i) %s' % (c, aItem)
# 	c=c+1

# query = raw_input('\nSelect Segmemt/s for download (Separate by comma) or A for all: ')

# Selection_Index = query.split(',')


# ID_Selection = []
# for aItem in Selection_Index:
# 	try:
# 		aItem = int(aItem)
# 		# Index for lists start at 0, started at 1 in selection index
# 		# must subtract by 1 to equal true index of selection.
# 		I_Val = aItem-1
# 		for aValue in SegmentSections[I_Val]:
# 			ID_Selection.append(aValue)
# 	except:
# 		try:
# 			aItem = str(aItem)
# 			if aItem == 'A':
# 				print '\nAll Selected...'
# 		except:
# 			print 'Fail'
# 			sys.exit()

# if aItem=='A':
# 	for aVar in IDs:
# 		ID_Selection.append(aVar)

# ID_Selection = LakeSak
# ID_Selection = []
# for x in range(1,44):
# 	ID_Selection.append(x)
ID_Selection = Proj

################################ QUERY END

# Run write segs

RunSegs = 'Yes'
# Will write individual shapefiles for each segment
# with fiona.open(SF_Dir+'LakeParsed.shp') as File:
if RunSegs == 'Yes':
	with fiona.open(SF_Dir+'RiverFull.shp') as File:

		# with fiona.open(SF_Dir+'ParsedLake.shp') as File:
			# schema = File.schema
		meta = File.meta
		# print meta
		# sys.exit()
		for feature in File:
			# ID = feature['properties']['NAME']
			ID = feature['properties']['Id']
			for aItem in ID_Selection:
				aItem = int(aItem)
				if aItem==int(ID):
					# Write folders for all the data
					dir = PSD+'Segment'+str(ID)
					if not os.path.exists(dir):
						os.makedirs(dir)
					TempDir = PSD+'Segment'+str(ID)+'/'
					# Write the shapefile
					# with fiona.open(TempDir+str(Number)+'SegShp', 'w', crs=from_epsg(26914), driver='ESRI Shapefile', schema=schema) as output:
					with fiona.open(TempDir+'Seg'+str(ID)+'Shp', 'w', **meta) as output:
							output.write(feature)
	File.close()

	# Write file to specify the specific parameters for this analysis
	with open('Analysis_Files/Analysis_Meta.txt', 'w') as Doc:
		for ID in ID_Selection:
			# 5/3/2020 Desktop Mod on ID, changed to have correct values.
			# Required the addition of str() around ID
			Doc.write('ID:')
			if ID != ID_Selection[-1]:
				Doc.write(str(ID)+',')
			if ID == ID_Selection[-1]:
				Doc.write(str(ID))
		Doc.write('-')

# By this point, the segment folders have been written

####################################### PARSE IMAGE FILES ##########################

# All image files for the given year need to be clipped and modified

# Use the selectscene file
# import SelectScene
# from SelectScene import*
# SelectScenes('Predict')


# Pull the appropriate image file names based upon P/R and date
Data = []
# Use from training data section
with open(Meta_Dir+'/Analysis_Meta.txt', 'r') as Doc:
	for aItem in Doc:
		Data.append(aItem)
Doc.close()

# Parse the data and relate segment ID with PR values
Sp1 = Data[0].split('-')
IDs = Sp1[0].split(',')
PRs = Sp1[1].split(',')
Filtered = []

for aItem in PRs:
	Sp1 = aItem.split('.')
	for aValue in ID_Selection:
		if int(Sp1[0][3::]) == aValue:
			Filtered.append([int(Sp1[0][3::]), Sp1[1]])

# Sort through image files and determine which ones relate to the segment id

# List the available images in the image folder
dir_path = str(os.path.dirname(os.path.realpath(__file__)))
TifDir = dir_path+'/TrainingData/TIF_Code_Images/'
SF_Dir = dir_path+'/TrainingData/Shapefile/'

Images = os.listdir(TifDir)
# print Images

for aItem in Filtered:
	Temp = []
	# Find appropriate images for path, row, date
	for ImageFolder in Images:
		Sp1 = ImageFolder.split('_')
		PRIF = Sp1[1]
		Year = int(Sp1[-1][0:4])
		if PrevYear==Year:
			if PRIF==aItem[-1]:
				Temp.append(ImageFolder)
	aItem.append(Temp)

# At this point, the segment, pr, and associated images are in the Filtered list
# Start filtering through the filtered dataset and pull the images to write
# habitat figures and data
PLD = 'ParsedSegmentData/'
Cloud = 'Yes'
################################## CLOUDS AND DISTURBANCE DATA #####################
if Cloud=='Yes':
	ProgBarLimit = len(Filtered)
	# Coverage data for the analysis
	CoverAll = []
	# ProgBarLimit = 2
	print '\nCloud Segment Analysis...\n'
	for r in tqdm.tqdm(range(ProgBarLimit)):
		aSeg = Filtered[r][0]
		SegCount = str(Filtered[r][0])
		# SegCount = 31
		# SegCount = 0
		SHPfl = 'Segment%s/Seg%sShp/Seg%sShp.shp' % (SegCount, SegCount, SegCount)
		SHPflFolder = 'Segment%s/Seg%sShp/' % (SegCount, SegCount)
		# Loop through the lake segments and write data
		###########################################################################
		#							Calculate Segement Cloud Coverage
		CPR = Filtered[r][1]

		# Identify folder paths within image output
		# LANDSAT_PATH = glob('TIF_Code_Images*/*')
		LANDSAT_PATH = []
		for aItem in Filtered[r][-1]:
			LANDSAT_PATH.append(TifDir+aItem)
		# List of folders containing files
		Paths = []
		# Files within each folder
		Files = []
		# Define folder paths for NDVI out for multiple dates
		Cloud_Out_Paths = []
		# Define sat codes
		SatCodes = []
		# Define dates
		Dates = []
		# Select available image files with correct path and row
		PRCodes = []
		# they will all be processed
		for aItem in LANDSAT_PATH:
			# Select the path and row for the image collection to 
			# compare to the CorrectPR object
			S1 = aItem.split('\\')
			S = S1[-1].split('/')

			SCode = str(S[-1])
			SCode = SCode.split('_')
			# Satellite
			SAT = SCode[0]
			# Date
			Date = SCode[-1]
			# Path and Row
			PR = SCode[1]
			# Compare the PR code to the CorrectPR
			# If the item matches, then the image file
			# for that PR will be selected for the NDVI analysis
			if PR==CPR:
				Paths.append(aItem)
				L = listdir(aItem)
				Files.append(L)
				SatCodes.append(SAT)
				Dates.append(Date)
				PRCodes.append(PR)

		# Write output paths for cloud assessment
		# General output for cloud data
		dir = 'ParsedSegmentData/Segment%s/Seg%sCLOUD' % (SegCount, SegCount)
		if not os.path.exists(dir):
			os.makedirs(dir)
		# Shapefile Directory
		dir = 'ParsedSegmentData/Segment%s/Seg%sCLOUD/SHP_Files' % (SegCount, SegCount)
		if not os.path.exists(dir):
			os.makedirs(dir)
		# tif clip Directory
		dir = 'ParsedSegmentData/Segment%s/Seg%sCLOUD/TIF_Files' % (SegCount, SegCount)
		if not os.path.exists(dir):
			os.makedirs(dir)
		# Tif clips/shp files
		dir = 'ParsedSegmentData/Segment%s/Seg%sCLOUD/PNG_Files' % (SegCount, SegCount)
		if not os.path.exists(dir):
			os.makedirs(dir)

		# Define the output folders for each data type 
		CloudPNG = 'ParsedSegmentData/Segment%s/Seg%sCLOUD/PNG_Files/' % (SegCount, SegCount)
		CloudSHP = 'ParsedSegmentData/Segment%s/Seg%sCLOUD/SHP_Files/' % (SegCount, SegCount)
		CloudTIF = 'ParsedSegmentData/Segment%s/Seg%sCLOUD/TIF_Files/' % (SegCount, SegCount)
		CloudGEN = 'ParsedSegmentData/Segment%s/Seg%sCLOUD/' % (SegCount, SegCount)


		# Write a general metadata file for the cloud covereage
		with open(CloudGEN+'Image_Meta.txt', 'w') as Doc:
			Doc.close()

		# Write cloud coverage data for each image that will be used
		# Cloud quality data is included in landsat image folders ending
		# with pixel_qa.tif. This will be done over all image dates in the
		# paths provided for this segment.
		Index = 0
		for Path in Paths:
			# Temp list for the coverage variables
			TempCov = []
			# Find the cloud pixel quality tif image in select path
			FileEnd = 'pixel_qa'
			CloudCover = glob(Path+'*/*%s.TIF' % FileEnd)
			CloudCover = CloudCover[0]
			# Write bound file based on parsed lake data

			with fiona.open(SF_Dir+'RiverFull.shp', "r") as shapefile:
				for feature in shapefile:
					if str(feature['properties']['Id'])==SegCount:
						source_schema = shapefile.schema
						crs = shapefile.crs
						source_properties = shapefile['properties']
						# Identify the location of the bound file for the DEM segment
						BoundFile = PLD+'Segment'+SegCount+'/'+'Seg'+SegCount+'Shp/'+SegCount+'Bounds.shp'
						with fiona.open(PLD+'Segment'+SegCount+'/'+'Seg'+SegCount+'Shp/'+SegCount+'Bounds.shp', 'w', 
							driver = shapefile.driver, schema = source_schema, crs = crs) as NewFile:
							# NewFile['geometry'] = feature['geometry']
							# NewFile['properties'] = source_properties
							# for f in feature:
							NewFile.write(feature)
							NewFile.close()
			shapefile.close()

			# Use a temp open file to write an extent shape 
			# this will be used to clip all tif images and must be converted
			# to match the crs of the tif images downloaded
			Temp = CloudCover
			# Use the extent of the BoundFile for the DEM data
			crop_extent = gpd.read_file(BoundFile)
			Temp = rasterio.open(Temp)
			crop_extent = crop_extent.to_crs(Temp.crs)

			# Save the crop extent as the segment extent
			Save = crop_extent.to_file(PLD+SHPflFolder+str(SegCount)+'Extent'+'.shp')

			# Open the extent file to clip the pixel_qa band by geometry
			import rasterio.mask
			with fiona.open(PLD+SHPflFolder+str(SegCount)+'Extent'+'.shp', "r") as shapefile:
				shapes = [feature["geometry"] for feature in shapefile]
			shapefile.close()

			with rasterio.open(CloudCover) as src:
				out_image_c, out_transformc = rasterio.mask.mask(src, shapes, crop=True)
				out_meta_c = src.meta
				# Define the parameters of the tif clip
				out_meta_c.update({"driver": "GTiff",
		                 "height": out_image_c.shape[1],
		                 "width": out_image_c.shape[2],
		                 "transform": out_transformc})
			src.close()

			# Write tif clip of cloud coverage data
			with rasterio.open(CloudTIF+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+str(SegCount)+'Clip.tif', 'w', **out_meta_c) as dst:
			    dst.write(out_image_c)
			dst.close()

			# FPlot is for future use in the script. This will be used in the final
			# png graphing output.
			FPlot = out_image_c

			from osgeo import osr
			gdal.UseExceptions()
			#  get raster datasource
			src_ds = gdal.Open(CloudTIF+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+str(SegCount)+'Clip.tif')
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
			dst_layername = 'CVar'
			drv = ogr.GetDriverByName("ESRI Shapefile")
			dst_ds = drv.CreateDataSource(CloudSHP+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+str(SegCount)+'pixel_qa' + ".shp" )
			dst_layer = dst_ds.CreateLayer(dst_layername, srs = srs)

			field_defn = ogr.FieldDefn("Value", ogr.OFTReal)
			dst_layer.CreateField(field_defn)
			gdal.Polygonize( srcband, None, dst_layer, 0, callback=None )

			src_ds = None
			srcband = None
			srs = None
			drv = None
			dst_ds = None
			dst_layer = None
			field_defn = None
			dst_layer = None
			########################## Output PNG for data ################################

			# Use gdal to open raster then extract unique values in raster clip with numpy
			# ds = gdal.Open(CloudTIF+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+str(SegCount)+'Clip.tif')
			# band =  ds.GetRasterBand(1)
			# array = np.array(band.ReadAsArray())
			# ds = None
			
			################## AREA CALCS
			# # Write lists of values from the pixel_qa array
			# # Here we will extract key variables and thier occurance in the 
			# # array
			# values, Freq = np.unique(array, return_counts=True)
			# values = values.tolist()
			# print values
			# sys.exit()
			# Freq = Freq.tolist()
			# # Drop no data value is present
			# i = 0
			# Loc = []
			# if 0 in values:
			# 	for aItem in values:
			# 		if aItem==0:
			# 			Loc.append(i)
			# 		i = i+1
			# if len(Loc)!=0:
			# 	values = values[1::]
			# 	Freq = Freq[1::]	
			# # This will convert values to scientific notation to simplify graphing
			# # this is being added due to large gaps in numbers in the data set
			def format_e(n):
				a = '%E' % n
				return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]

			# # This will reformat the 2 lists. ConValues will now be in sci notation
			# # ConFreq into float
			# ConValues = []
			# ConFreq = []
			# i = 0
			# while i < len(values):
			# 	Vfloat = values[i]
			# 	TempV = format_e(Vfloat)
			# 	Ffloat = Freq[i]
			# 	TempF = Ffloat
			# 	ConValues.append(TempV)
			# 	ConFreq.append(TempF)
			# 	i = i+1

			# # print ConValues
			# # print ConFreq
			# # sys.exit()

			# # Calculate coverage area for each pixel type
			# # Coverage will contain the [pixel value, area]
			# Coverage = []
			# # CovVal = []
			# i = 0
			# TotalArea = 0
			# while i < len(values):
			# 	# Skip 1 for calculations
			# 	Pixels = Freq[i]
			# 	# Area of a pixel is 900 square meters 30*30
			# 	Area = Pixels*900
			# 	# Do no include 1 in total area calc
			# 	if values[i]!=1:
			# 		TotalArea = TotalArea+Area
			# 	# TotalArea = TotalArea+Area
			# 	InsertItem = [values[i], Area]
			# 	Coverage.append(InsertItem)
			# 	i = i+1
			############# AREA CALC MOD
			# This will calculate based on polygon. Previous version was on 
			# number of pixels, estimating area single pix to be 900^m2. This 
			# uses polygon area

			# MOD 9/26/2020
			# Define IDs
			File = CloudSHP+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+str(SegCount)+'pixel_qa' + ".shp"
			TotalArea = 0
			valuesAll = []
			with fiona.open(File, "r") as src:
				for Polygon in src:
					# Define the id for a given poly
					ID = int(Polygon['properties']['Value'])
					Geom = Polygon['geometry']['coordinates'][0]
					Geom = tuple(Geom)
					# make a polygon
					from shapely.geometry import Polygon
					Po = Polygon(Geom)
					# calculate area (meters)
					Area = Po.area
					if ID!=1:
						TotalArea = TotalArea+Area
					valuesAll.append(ID)
			src.close()

			values = list(collections.Counter(valuesAll).keys())
			# Sort data to fix the labeling issue
			values.sort()



			# Open the interference shapefile for clip
			Coverage = []
			# Iterate through Values
			for Val in values:
				TempArea = 0
				with fiona.open(File, "r") as src:
					ESPG = src.crs
					for Polygon in src:
						# Define the id for a given poly
						ID = int(Polygon['properties']['Value'])
						if ID == Val:
							# Omit the data we do not want to clip. Such as land, water, saturation
							# and the image edge.
							Geom = Polygon['geometry']['coordinates'][0]
							Geom = tuple(Geom)
							# make a polygon
							from shapely.geometry import Polygon
							Po = Polygon(Geom)
							# calculate area (meters)
							Area = Po.area
							TempArea = TempArea+Area

				Coverage.append([Val, TempArea])
				src = None



			# Calculate the coverage of each datatype. This will vary from each 
			# satellite used. We will use this section to conduct a cloud coverage
			# assesment for the chunk of data we are actually using for the image.
			# Do not include 1 in the percent coverage, it does not contain
			# relevant data. It is the edges of the tif clip which fill the 
			# map figure. Basically no data, but is technically an item that is 
			# callable.
			for Value in Coverage:
				if Value[0]!=1:
					PercentCoverage = float(Value[-1])/float(TotalArea)
					Value.append(PercentCoverage)

			if SatCodes[Index] == 'LC08':
				Clear = 322
				Water = 324
			if SatCodes[Index] == 'LT05':
				Clear = 66
				Water = 68
			# Define Labels
			labels = []
			sizes = []
			i = 1
			CovI = 0
			Matched = False
			# Add data to temp coverage set. This will be used later for NDVI.
			# If complete interference is occuring, NDVI cannot be calculated.
			TempCov.append(SegCount)
			# TempCov.append(str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index]))
			# Add check variables to coverage. If a specific type of variable
			# is not present, it needs to be represented as a zero in the 
			# coverage final output file.
			CC = 0
			WW = 0
			IN = 0
			# ND = 0
			TotIN = 0

			for aValue in Coverage:
				# print aValue
				if aValue[0]==Clear:
					Percent = float(Coverage[CovI][-1])*100
					PercentRound = str(round(Percent, 2))
					# Add data to temp coverage set
					CC = copy.deepcopy(PercentRound)
					labels.append(('Clear '+'%s'+'%%') % PercentRound)
					Matched = True
				if aValue[0]==Water:
					Percent = float(Coverage[CovI][-1])*100
					PercentRound = str(round(Percent, 2))
					# Add data to temp coverage set
					WW = copy.deepcopy(PercentRound)
					labels.append(('Water '+'%s'+'%%') % PercentRound)
					Matched = True
				# Pass the image edge int value of 1
				if aValue[0]==1:
					Matched = True
				if Matched==False:
					Percent = float(Coverage[CovI][-1])*100
					PercentRound = str(round(Percent, 2))
					# Add data to temp coverage set
					IN = copy.deepcopy(PercentRound)
					TotIN = round(Percent, 2)+TotIN
					labels.append(('%i_INT '+'%s'+'%%') % (aValue[0], PercentRound))
					i = i+1
				Matched = False
				CovI = CovI+1
				# Do not append the size of this, otherwise it will be added to pie
				if aValue[0]!=1:
					sizes.append(aValue[-1])
			# Add the Clear land, water, and interference data to the temp dataset	
			TempCov.append(CC)
			TempCov.append(WW)
			TempCov.append(TotIN)
			# TempCov Format: [Segment, Landsat info, Land, Water, Interference]

			#	Data for text box on plot ax2
			# Clip image properties statements: *****************
			# Define cover of clear pixels
			TotalInt = 0
			for aItem in Coverage:
				if aItem[0]==Water:
					# TotalClear = TotalClear+aItem[-1]
					pass
				if aItem[0]==Clear:
					# TotalClear = TotalClear+aItem[-1]
					pass
				# Do not include values of 1, which are the edges, to be included in
				# area of interference cal below. It is still technically an item,
				# but contains no data.
				# if aItem[0]!=Clear and aItem[0]!=Water and aItem[0]!=1:
				# 	TotalInt = TotalInt+aItem[-1]
			######################## Not possible for this file format
			######################## Not possible for this file format
			
			# Area string
			Area = 'CLIP_AREA: %s^2 Meters' % str(format_e(TotalArea))
			# Total calculated interference
			TotalInt = str(TotIN)
			interference = 'ESTIMATED_CLIP_INTERFERENCE: %s' % (TotalInt)

			Cinput = 'Clip Image Properties:'+'\n'+interference+'\n'+Area

			# Set colors and labels for the ax1 and ax3 figures

			# Change numpy array values for clean plotting
			# Open the complete data set for the conversion. Insert an
			# integer to replace the values in the NP array with
			# ReplaceInt = 1
			# Given the nature of our problem with labeling being off, this is now
			# set to 0. It prevents the data in the map from being offset and/or
			# incorrect
			ReplaceInt = 0
			for aItem in Coverage:
				aItem.append(ReplaceInt)
				ReplaceInt=ReplaceInt+1
			# Use the appended list item to match values in the np array
			# then replace them with smaller integers
			for aItem in Coverage:
				FPlot = np.where(FPlot==aItem[0], aItem[-1], FPlot)

			Colors = []
			# Random color map colors, will use 2 separate color map schemes so they
			# do not match. Every other one will be from a separate scheme.
			# The two schemes
			jet = cm.get_cmap('jet', len(Coverage)+1)
			binary = cm.get_cmap('binary', len(Coverage)+1)

			# Use the inserted value to determine color
			bounds = []
			for aItem in Coverage:
				# print aItem
				Var = int(1+aItem[-1])
				# Check if index is even
				if aItem[-1]%2==0:
					if aItem[0]!=Water and aItem[0]!=Clear:
						RandomColor = jet(Var)
				# Check if index is odd
				if aItem[-1]%2!=0:
					if aItem[0]!=Water and aItem[0]!=Clear:
						RandomColor = binary(Var)
				if aItem[0]==Water:
					RandomColor = 'blue'
				if aItem[0]==Clear:
					RandomColor = 'green'
				# Add color for the "no data" which is just the area outside of the
				# clipped tif data
				if aItem[0]==1:
					RandomColor = 'white'

				bounds.append(aItem[-1])
				Colors.append(RandomColor)
			# print bounds
			# print Colors 
			# sys.exit()
			# Copy the colors as is for the pie chart
			# ColorsPie = copy.deepcopy(Colors)
			ColorsPie = copy.deepcopy(Colors)
			# Insert a value for nodat on the ax1 plot
			# Colors.insert(0, 'white')
			# Colors.insert(1, 'white')
			# print Colors
			# print len(Colors)

			# Define bounds on the ax1 plot
			Colors = LinearSegmentedColormap.from_list('Custom cmap', Colors)
			bounds = np.linspace(0, len(ColorsPie), len(ColorsPie)+1)
			# bounds = np.linspace(0, len(Colors), len(Colors)+2)
			norm = colors.BoundaryNorm(bounds, Colors.N)

			# Plot the data over a subplot2grid setup
			gridsize = (5,2)
			fig = plt.figure(figsize=(12,9))
			ax1 = plt.subplot2grid(gridsize, (0, 0), colspan= 2, rowspan= 3)
			ax2 = plt.subplot2grid(gridsize, (3, 0), colspan = 1, rowspan=2)
			ax3 = plt.subplot2grid(gridsize, (3, 1), colspan = 1, rowspan = 2)

			# Plot the raster image on ax1
			ax1.set_title(('Pixel Quality (pixel_qa) Values Segment %s' % SegCount), fontsize='14')
			out_image_c = np.rollaxis(FPlot[0], 0, 1)
			out_image_c = ax1.imshow(out_image_c, cmap=Colors, norm=norm)

			# Modify color bar to fit size of graph
			aspect = 20
			pad_fraction = 0.5
			divider = make_axes_locatable(ax1)
			width = axes_size.AxesY(ax1, aspect=1./aspect)
			pad = axes_size.Fraction(pad_fraction, width)
			cax = divider.append_axes("right", size=width, pad=pad)

			boundary_means = [np.mean([norm.boundaries[ii], norm.boundaries[ii - 1]])
			                  for ii in range(1, len(norm.boundaries))]
			# Modify the labels for ax1 colorbar
			a1L = []
			for aItem in labels:
				aItem=aItem.split(' ')
				a1L.append(aItem[0])
			# Insert a no data label
			# Add additional label for no data (values of 1)
			a1L.insert(0, 'No Data')
			# a1L.insert(1, 'No Data')
			# print a1L
			# print len(a1L)
			# # print Colors
			# sys.exit()

			# Plot the colorbar and set ticklabels
			plt.colorbar(out_image_c, cax=cax, ticks=boundary_means).set_ticklabels(a1L)
			
			# Insert sat, path/row, and date information
			Input = SatCodes[Index] + '\n' + PR + '\n' + Dates[Index]
			plt.text(0,1.0025, Input, fontsize='7')

			# Plot the image information
			# textstr = Finput+Cinput
			textstr = Cinput
			props = dict(boxstyle='round', facecolor='wheat')
			# Hide the empty graph shhh...
			ax2.spines['bottom'].set_color('white')
			ax2.spines['top'].set_color('white') 
			ax2.spines['right'].set_color('white')
			ax2.spines['left'].set_color('white')
			ax2.tick_params(axis='x', colors='white')
			ax2.tick_params(axis='y', colors='white')
			# Plot the text
			ax2.text(0, 1, textstr, fontsize=10,
	    		verticalalignment='top', bbox=props)

			# Plot the pie of percentage of cover
			# Select all color values outside first to accomodate for added white
			# since there is no "No Data" info due to how the images are processed,
			# that area is technically part of the array data. So it must be accounted
			# for for the plot to split it for spatial analysis and must be removed
			# from the total area and disturbance calculations. 9/21/20
			patches, texts = ax3.pie(sizes, startangle=90, colors=ColorsPie[1::])

			ax3.legend(patches, labels, loc="upper right", fontsize='6')
			ax3.axis('equal')
			ax3.set_title('Percent of Coverage (pixel_qa)', fontsize='8')

			plt.tight_layout(pad=4, w_pad=4, h_pad=4)
			plt.savefig((CloudPNG+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+SegCount+'CloudCover.png'), dpi=None, 
			    facecolor='w', edgecolor='b', orientation='portrait', papertype=None, 
			    format=None, transparent=False, bbox_inches=None, pad_inches=0.1, 
			    frameon=None)

			# plt.show()
			# print str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+str(SegCount)
			plt.close()
			# sys.exit()

			# Define the insertable item to be appended to the meta file.
			# Insert the satellite data into the coverage meta insertable item

			# Extract estimated water and land cover:
			# Coverage data format: [pixel val, area, % cover, index val]
			WaterArea = 0
			LandArea = 0
			for aItem in Coverage:
				if aItem[0]==Water:
					WaterArea = aItem[1]
				if aItem[0]==Clear:
					LandArea = aItem[1]
			
			SatCode = str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])
			InsertItem = SatCode+'^'+TotalInt+'^'+str(TotalArea)+'^'+str(WaterArea)+'^'+str(LandArea)+'*'
			
			# Append the information from the analysis to the image meta file in the 
			# cloud folder for that given segment. The format of the output will be 
			# satellite information, segment location, and clip coverage.
			with open(CloudGEN+'Image_Meta.txt', 'a') as Doc:
				Doc.write(InsertItem)
			Index = Index+1
			Doc.close()
			TempCov.append(SatCode)
			CoverAll.append(TempCov)
		# sys.exit()

	# Write the list as a pandas dataframe, then use the pandas 
	# module to write that to a csv file.
	df = pd.DataFrame(CoverAll)
	df.to_csv(Meta_Dir+'SegmentImageInfo.csv')


	####################################################################################
						# Write shapefiles from the cloud raster data

	CloudFolderPaths = []

	for Segment in IDs:
		# Define the path for cloud shapefiles
		Path = (PLD+'Segment%s/'+'Seg%sCLOUD/'+'SHP_Files') % (Segment, Segment)
		CloudFolderPaths.append(Path)

	Temp = []
	CloudShapesO = []

	# Write list of available shapefiles in each folder. If multiple dates
	# are being processed, there will be multiple shp files available.
	for Folder in CloudFolderPaths:
		Files = glob(Folder+'*/*.shp')
		CloudShapesO.append(Files)

	# Open cloud cover shapefiles
	for Collection in CloudShapesO:
		for SHP in Collection:
			if SHP[-12::]== 'pixel_qa.shp':
				
				# Split up the shapefile name, the satellite and date info is included.
				# We use this to determine which values are associated with clear coverage
				Sp1 = SHP.split('/')
				Sp2 = Sp1[-1].split('\\')
				Sp3 = Sp2[-1].split('_')
				Sat = Sp3[0]
				# Define folder path from splits
				TID = Sp1[1][7::]
				FolderLoc = Sp1[0]+'/'+Sp1[1]+'/'+Sp2[0]+'/'+Sp2[1]+'/'
				Name = FolderLoc+Sp3[0]+'_'+Sp3[1]+'_'+Sp3[2]+'_'+'%sINT.shp' % (TID)
				# Define our clear values
				if Sat=='LC08':
					Water = 324
					Clear = 322
				if Sat=='LT05':
					Clear = 66
					Water = 68
				with fiona.open(SHP, 'r') as source:
					meta = source.meta
					with fiona.open(Name, 'w', **meta) as output:
						for feature in source:
							if feature['properties']['Value'] == Clear:
								output.write(feature)
							if feature['properties']['Value'] == Water:
							 	output.write(feature)
				source.close()
				output.close()

############################# END OF CLOUD/DISTURBANCE COVER ####################

###################################### NDVI Analysis #####################################
# This is not the prediction section, it simply provides outputs explaining
# the conditions and figures supporting it.

# Open the coverage quality meta file. This will contain information
# explaining the intereference within a given segment.

ErrorLog = []

MetaCovRaw = []

with open(Meta_Dir+'SegmentImageInfo.csv', 'r') as File:
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

ProgBarLimit = len(Filtered)
HabitatMetaData = []
print '\nNDVI Segment Analysis...\n'
NDVIRUN = 'yes'
if NDVIRUN == 'yes':
	for r in tqdm.tqdm(range(ProgBarLimit)):
		try:
			# r =1
			aSeg = Filtered[r][0]
			SegCount = str(Filtered[r][0])

			SHP_Extent_File = 'Segment%s/Seg%sShp/%sExtent.shp' % (SegCount, SegCount, SegCount)
			# Write output paths for cloud assessment
			# General output for cloud data
			dir = 'ParsedSegmentData/Segment%s/Seg%sNDVI' % (SegCount, SegCount)
			if not os.path.exists(dir):
				os.makedirs(dir)
			# tif clip Directory
			dir = 'ParsedSegmentData/Segment%s/Seg%sNDVI/TIF_Files' % (SegCount, SegCount)
			if not os.path.exists(dir):
				os.makedirs(dir)
			# PNG files
			dir = 'ParsedSegmentData/Segment%s/Seg%sNDVI/PNG_Files' % (SegCount, SegCount)
			if not os.path.exists(dir):
				os.makedirs(dir)
			dir = 'ParsedSegmentData/Segment%s/Seg%sNDVI/SHP_Files' % (SegCount, SegCount)
			if not os.path.exists(dir):
				os.makedirs(dir)
			# Define the output folders for each data type 
			NDVIPNG = 'ParsedSegmentData/Segment%s/Seg%sNDVI/PNG_Files/' % (SegCount, SegCount)
			NDVITIF = 'ParsedSegmentData/Segment%s/Seg%sNDVI/TIF_Files/' % (SegCount, SegCount)
			NDVISHP = 'ParsedSegmentData/Segment%s/Seg%sNDVI/SHP_Files/' % (SegCount, SegCount)
			# try:

			CPR = Filtered[r][1]
		 	# Identify folder paths within image output
			# LANDSAT_PATH = glob('TIF_Code_Images*/*')
			LANDSAT_PATH = []
			for aItem in Filtered[r][-1]:
				LANDSAT_PATH.append(TifDir+aItem)
			# List of folders containing files
			Paths = []
			# Files within each folder
			Files = []
			# Define folder paths for NDVI out for multiple dates
			Cloud_Out_Paths = []
			# Define sat codes
			SatCodes = []
			# Define dates
			Dates = []
			# Select available image files with correct path and row
			PRCodes = []
			# they will all be processed
			for aItem in LANDSAT_PATH:
				# Select the path and row for the image collection to 
				# compare to the CorrectPR object
				S1 = aItem.split('\\')
				S = S1[-1].split('/')

				SCode = str(S[-1])
				SCode = SCode.split('_')
				# Satellite
				SAT = SCode[0]
				# Date
				Date = SCode[-1]
				# Path and Row
				PR = SCode[1]
				# Compare the PR code to the CorrectPR
				# If the item matches, then the image file
				# for that PR will be selected for the NDVI analysis
				if PR==CPR:
					Paths.append(aItem)
					L = listdir(aItem)
					Files.append(L)
					SatCodes.append(SAT)
					Dates.append(Date)
					PRCodes.append(PR)


			Index = 0
			for Path in Paths:
				# We are not calculating NDVI now. Select the NDVI band
				NDVIBand = 'ndvi'

				# ****FILTER BY COVERAGE****
				# Only process records with favorable coverage conditions.
				# Based upon cloud cover, currently set at 10%, see below.
				SkipRecord = False
				# Parse the path to extract date from image
				Mod1 = Path.split('\\')
				Mod2 = Mod1[-1].split('_')
				DateCodeCheck = Mod2[-1]
				# MetaCov Format: [Segment, Clear land, Water, Interference, Landsat Code, Date Code]
				for aRow in MetaCov:
					if aRow[0]==SegCount:
						# MetaCov aRow[-1] was converted to int to drop the
						# new line text (\n). Both of these variables below are
						# missing leading zeros.
						if str(aRow[-1])==str(DateCodeCheck):
							if float(aRow[3])>10.0:
								SkipRecord = True
								Index = Index+1

				# WORK WORK
				if SkipRecord==False:

					for aItem in glob(Path+'*/*%s.TIF' % NDVIBand):
						NDVIImage = aItem
						TempNDVIim = aItem
					import rasterio.mask
					with fiona.open(PSD+SHP_Extent_File, "r") as shapefile:
						shapes = [feature["geometry"] for feature in shapefile]

					with rasterio.open(NDVIImage) as src:
						out_image_B1, out_transform1 = rasterio.mask.mask(src, shapes, crop=True)
						out_meta_B1 = src.meta
						# Define the parameters of the tif clip
						out_meta_B1.update({"driver": "GTiff",
				                 "height": out_image_B1.shape[1],
				                 "width": out_image_B1.shape[2],
				                 "transform": out_transform1})
					# Write a clipped B1 tif file
					with rasterio.open(NDVITIF+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+SegCount+'_'+'ndvi.tif', 'w', **out_meta_B1) as dst:
					    dst.write(out_image_B1)
						# print out_meta_B4
						# print out_transform4
					dst.close()
					################################################
					# Define the area of data to be removed
					# Open the shapefile containing no data shapes
					# WE DO NOT WANT 352 and 480, these seem to be bad calcs

					# Must write a new clip file without the 352 and 480 interference
					FolderLoc = 'ParsedSegmentData/Segment%s/Seg%sCLOUD/SHP_Files/' % (SegCount, SegCount)
					File = str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+SegCount+'pixel_qa.shp'
					
					
					# # Open the interference shapefile for clip
					# with fiona.open(FolderLoc+File, "r") as shapefile:
					# 	shapes = [feature["geometry"] for feature in shapefile]

					# MOD 9/26/2020
					ClipKeep = 0.0
					TotalAreaOfClip = 0
					# Open the interference shapefile for clip
					ShapeFileDF = pd.DataFrame(columns=['id', 'area', 'geometry'])
					WaterDF = pd.DataFrame(columns=['id', 'area', 'geometry'])
					with fiona.open(FolderLoc+File, "r") as src:
						ESPG = src.crs
						for Polygon in src:
							# Define the id for a given poly
							ID = int(Polygon['properties']['Value'])
							# Omit the data we do not want to clip. Such as land, water, saturation
							# and the image edge.
							Geom = Polygon['geometry']['coordinates'][0]
							Geom = tuple(Geom)
							# make a polygon
							from shapely.geometry import Polygon
							Po = Polygon(Geom)
							# calculate area (meters)
							Area = Po.area
							if ID != 1:
								TotalAreaOfClip = TotalAreaOfClip+Area
							if ID != 352 and ID != 480 and ID != 324 and ID != 322 and ID != 1:
								# Define insert item
								InsertItem = [ID, Area, Po]
								# Start a df so you can transpose it to the
								# correct format
								TempDF = pd.DataFrame(InsertItem)
								# Rotate the data to match our desired format
								TempDF = TempDF.transpose()
								# Insert the column names
								TempDF.columns = ['id', 'area', 'geometry']
								# Add to main df
								ShapeFileDF = ShapeFileDF.append(TempDF, ignore_index = True)
							# Check the area of image data which is not going to be clipped
							if ID == 352 or ID == 480:	
								ClipKeep=ClipKeep+Area
							# Write water shapefile
							if ID == 324:
								InsertItem = [ID, Area, Po]
								# Start a df so you can transpose it to the
								# correct format
								TempDF = pd.DataFrame(InsertItem)
								# Rotate the data to match our desired format
								TempDF = TempDF.transpose()
								# Insert the column names
								TempDF.columns = ['id', 'area', 'geometry']
								WaterDF = WaterDF.append(TempDF, ignore_index = True)

					src = None

					######################### USE WATER SHAPE TO CLIP SHORE #################
					################################# MAX DISTANCE USED ####################
					# Write water shapefile
					gdf = gpd.GeoDataFrame(WaterDF, geometry='geometry')
					gdf.crs = ESPG
					gdf.to_file(FolderLoc+"WaterData.shp")
					gdf = None
					# Add the buffer to it
					ShapeFileDF = pd.DataFrame(columns=['id', 'geometry'])
					with fiona.open(FolderLoc+"WaterData.shp") as src:
						ESPG = src.crs
						# print '\nAdding Buffers...\n'
						ProgBarLimit = len(src)
						for i in tqdm.tqdm(range(ProgBarLimit)):
							Polygon = src[i]
							# Define the id for a given poly
							ID = int(Polygon['properties']['id'])
							Geom = Polygon['geometry']['coordinates'][0]
							Geom = tuple(Geom)
							from shapely.geometry import Polygon
							Po = Polygon(Geom)
							PRP = 'No'
							for aValue in RiverProj:
								if aValue==int(aSeg):
									PRP = 'Yes'
							if PRP == 'Yes':
								BFP = Po.buffer(MaxDistance*3)
							if PRP == 'No':
								BFP = Po.buffer(MaxDistance)
							InsertItem = [ID, BFP]
							TempDF = pd.DataFrame(InsertItem)
							# Rotate the data to match our desired format
							TempDF = TempDF.transpose()
							# Insert the column names
							TempDF.columns = ['id', 'geometry']
							ShapeFileDF = ShapeFileDF.append(TempDF, ignore_index = True)

					# Convert the dataframe to a geodataframe
					gdf = gpd.GeoDataFrame(ShapeFileDF, geometry='geometry')
					# No refrence system is in this file yet, we must set it to the correct one
					# before changing it.
					gdf.crs = ESPG
					# Write a shapefile from the geodataframe object
					gdf.to_file(FolderLoc+"WaterBuff.shp")
					src.close()

					# Water buffer is now written
					# rewrite the tif image file to only look for pixels near water
					with fiona.open(FolderLoc+"WaterBuff.shp", "r") as shapefile:
						shapes = [feature["geometry"] for feature in shapefile]
					shapefile.close()
					# Use the data mask
					with rasterio.open(NDVITIF+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+SegCount+'_'+'ndvi.tif', 'r') as src2:
						out_image_B2, out_transform2 = rasterio.mask.mask(src2, shapes, crop=True)
						out_meta_B2 = src2.meta
						out_meta_B2.update({"driver": "GTiff",
								"height": out_image_B2.shape[1],
								"width": out_image_B2.shape[2],
								"transform": out_transform2})
					src.close()
					with rasterio.open(NDVITIF+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+SegCount+'_'+'ndvi.tif', 'w', **out_meta_B2) as dst:
					    dst.write(out_image_B2)
					dst.close()

					# Remove all the water shapefiles now
					FilesInShp = list(os.listdir(FolderLoc))
					DelFiles = []
					for aFile in FilesInShp:
						# Sp1 = aFile.split('_')
						Sp2 = aFile.split('.')
						Type = Sp2[0]
						if Type == 'WaterBuff':
							DelFiles.append(FolderLoc+aFile)	
						if Type =='WaterData':
							DelFiles.append(FolderLoc+aFile)

					for aFile in DelFiles:
						os.remove(aFile)

					###################### New tif with distance buffer clip ##########################
					############################## MADE ABOVE ###########################
					# Copied code			
					# with fiona.open(PLD+SHPflFolder+str(SegCount)+'Extent'+'.shp', "r") as shapefile:
					# 	shapes = [feature["geometry"] for feature in shapefile]
					# shapefile.close()

					# with rasterio.open(CloudCover) as src:
					# 	out_image_c, out_transformc = rasterio.mask.mask(src, shapes, crop=True)
					# 	out_meta_c = src.meta
					# 	# Define the parameters of the tif clip
					# 	out_meta_c.update({"driver": "GTiff",
				 #                 "height": out_image_c.shape[1],
				 #                 "width": out_image_c.shape[2],
				 #                 "transform": out_transformc})
					# src.close()

					# # Write tif clip of cloud coverage data
					# with rasterio.open(CloudTIF+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+str(SegCount)+'Clip.tif', 'w', **out_meta_c) as dst:
					#     dst.write(out_image_c)
					# dst.close()

					# Determine interference clips
					Passed = 0
					if len(ShapeFileDF)>0:
						gdf = gpd.GeoDataFrame(ShapeFileDF, geometry='geometry')
						gdf.crs = ESPG
						gdf.to_file(FolderLoc+"ClipData.shp")
						gdf = None
						Passed = 1

					ShapeFileDF = None
					# If there are disturbances to clip, we will do this
					if Passed == 1:
						# Open the clip file
						with fiona.open(FolderLoc+"ClipData.shp", "r") as shapefile:
							shapes = [feature["geometry"] for feature in shapefile]
						# Use the data mask
						with rasterio.open(NDVITIF+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+SegCount+'_'+'ndvi.tif', 'r') as src2:
							out_image_B2, out_transform2 = rasterio.mask.mask(src2, shapes, crop=True)
							out_meta_B2 = shapefile.meta
						ndvi = np.array(out_image_B2)
						src2.close()
						out_image_B2 = None
						out_transform2 = None
						out_meta_B2 = None
						shapefile.close()

					if Passed == 0:
						from PIL import Image
						ndviIm = Image.open(NDVITIF+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+SegCount+'_'+'ndvi.tif')
						ndvi = np.asarray(ndviIm)
						ndviIm = None

					# from PIL import Image
					# ndviIm = Image.open(NDVITIF+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+SegCount+'_'+'ndvi.tif')
					# ndvi = np.asarray(ndviIm)
					# ndviIm = None
					ndviA = np.array(ndvi)*0.0001
					#######################################################################
					#							NDVI DATA EXTRACT
					# with rasterio.open(NDVITIF+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+SegCount+'_'+'ndvi.tif', 'r') as src2:
					# 	out_image_B2, out_transform2 = rasterio.mask.mask(src2, shapes, crop=True)
					# 	out_meta_B2 = src.meta
					# Open the NDVI tif and read values
					# ndvi = np.array(out_image_B2)
					from PIL import Image
					ndviIm = Image.open(NDVITIF+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+SegCount+'_'+'ndvi.tif')
					ndvi = np.asarray(ndviIm)
					ndviIm = None
					ndviA = np.array(ndvi)*0.0001


					# Write habitat data to tif
					###############################################
					# Habitat = ndviA[np.where((MinNDVI <= ndviA) & (ndviA <= MaxNDVI)) ] = 1
					CorrectFile = str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+SegCount+'_'+'ndvi.tif'
			 		driver = gdal.GetDriverByName('GTiff')
					ds = gdal.Open(NDVITIF+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+SegCount+'_'+'ndvi.tif')
					band =  ds.GetRasterBand(1)
					lista = np.array(band.ReadAsArray())*0.0001

					lista[np.where((MinNDVI <= lista) & (lista <= MaxNDVI)) ] = 1
					lista[np.where((MinNDVI > lista) & (lista > MaxNDVI)) ] = 2
					# create new file
					file2 = driver.Create(NDVITIF+'Working.tif', ds.RasterXSize , ds.RasterYSize , 1)
					file2.GetRasterBand(1).WriteArray(lista)

					# spatial ref system
					proj = ds.GetProjection()
					georef = ds.GetGeoTransform()
					file2.SetProjection(proj)
					file2.SetGeoTransform(georef)
					file2.FlushCache()
					ds = None
					file2 = None
					lista = None
					driver = None
					band = None
					proj = None

					# Write habitat to poly
					######################### Habitat Polygon

					gdal.UseExceptions()
					#  get raster datasource
					src_ds = gdal.Open(NDVITIF+'Working.tif')
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
					dst_layername = NDVISHP+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+SegCount+'_'+'ndvi.shp'
					drv = ogr.GetDriverByName("ESRI Shapefile")
					dst_ds = drv.CreateDataSource(dst_layername)
					dst_layer = dst_ds.CreateLayer(dst_layername, srs = srs)

					field_defn = ogr.FieldDefn("val", ogr.OFTReal)
					dst_layer.CreateField(field_defn)

					gdal.Polygonize(srcband, None, dst_layer, 0, callback=None )

					# In order to open this again after its creation, it needs to
					# be closed out completely.
					dst_layer = None
					srs = None
					dst_ds = None
					drv = None
					field_defn = None
					src_ds = None

					# Delete working tif file that was classified
					os.remove(NDVITIF+'Working.tif')

					# Write a habitat only file
					############################
					File = NDVISHP+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+str(SegCount)+'_'+'ndvi.shp'
					mainDat = str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+str(SegCount)+'_'
					src = None
					# Total Area of habitat
					AreaSum = 0
					# Filter out area as well
					ShapeFileDF = pd.DataFrame(columns=['id', 'area', 'geometry'])
					with fiona.open(File, "r") as src:
						ESPG = src.crs
						for Polygon in src:
							# Define the id for a given poly
							ID = int(Polygon['properties']['val'])
							if ID == 1:
								# Define geom
								Geom = Polygon['geometry']['coordinates'][0]
								Geom = tuple(Geom)
								# make a polygon
								from shapely.geometry import Polygon
								Po = Polygon(Geom)
								# calculate area (meters)
								Area = Po.area
								# This can remove single pixels or small patches 
								# from the final data set
								if Area>1000:
									AreaSum = AreaSum + Area
									# Define insert item
									InsertItem = [ID, Area, Po]
									# Start a df so you can transpose it to the
									# correct format
									TempDF = pd.DataFrame(InsertItem)
									# Rotate the data to match our desired format
									TempDF = TempDF.transpose()
									# Insert the column names
									TempDF.columns = ['id', 'area', 'geometry']
									# Add to main df
									ShapeFileDF = ShapeFileDF.append(TempDF, ignore_index = True)

					src = None
					gdf = gpd.GeoDataFrame(ShapeFileDF, geometry='geometry')
					gdf.crs = ESPG
					gdf.to_file(NDVISHP+mainDat+"Habitat.shp")
					gdf = None
					# Remove working shapefile
					# Define all NDVI files
					FilesInShp = list(os.listdir(NDVISHP))
					DelFiles = []
					for aFile in FilesInShp:
						Sp1 = aFile.split('_')
						Sp2 = Sp1[-1].split('.')
						Type = Sp2[0]
						if Type == 'ndvi':
							DelFiles.append(NDVISHP+aFile)
					# now remove them 
					# Done this way as there are multiple files to a shapefile
					for aFile in DelFiles:
						os.remove(aFile)

					# Define the estimated area of habitat
					HabArea = 0

					with fiona.open(NDVISHP+mainDat+"Habitat.shp") as src:
						for Polygon in src:
							# Define geom
							Geom = Polygon['geometry']['coordinates'][0]
							Geom = tuple(Geom)
							# make a polygon
							from shapely.geometry import Polygon
							Po = Polygon(Geom)
							# calculate area (meters)
							HabArea = Po.area+HabArea
					src.close()


					# BUFFER STUFF IF NEEDED

					# ShapeFileDF = pd.DataFrame(columns=['id', 'geometry'])
					# with fiona.open(Shapefile) as input:
					# 	print '\nAdding Buffers...\n'
					# 	ProgBarLimit = len(input)
					# 	for i in tqdm.tqdm(range(ProgBarLimit)):
					# 		Polygon = input[i]
					# 		# Define the id for a given poly
					# 		ID = int(Polygon['properties']['id'])
					# 		Geom = Polygon['geometry']['coordinates'][0]
					# 		Geom = tuple(Geom)
					# 		from shapely.geometry import Polygon
					# 		Po = Polygon(Geom)
					# 		BFP = Po.buffer(0.1)
					# 		InsertItem = [ID, BFP]
					# 		TempDF = pd.DataFrame(InsertItem)
					# 		# Rotate the data to match our desired format
					# 		TempDF = TempDF.transpose()
					# 		# Insert the column names
					# 		TempDF.columns = ['id', 'geometry']
					# 		ShapeFileDF = ShapeFileDF.append(TempDF, ignore_index = True)

					# # Convert the dataframe to a geodataframe
					# gdf = gpd.GeoDataFrame(ShapeFileDF, geometry='geometry')
					# # No refrence system is in this file yet, we must set it to the correct one
					# # before changing it.
					# gdf.crs = fiona.crs.from_epsg(26914)
					# # Write a shapefile from the geodataframe object
					# print 'Writting the shapefile: 567WithBuff_.1.shp'
					# gdf.to_file("567WithBuff_.1.shp")

					# #####################################################
					
					# Histogram bin data
					# Define output text information
					bins = [-1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.8, 0.9, 1]

					# sys.exit()
					MinVal = -1
					MaxVal = 1

					BinSelect = []

					for aValue in bins:
						if aValue>MinVal:
							if aValue<MaxVal:
								BinSelect.append(aValue)

					ImagePropertiesRaw = []
					ImageFile = 'Image_Meta.txt'
					CloudGEN = 'ParsedSegmentData/Segment%s/Seg%sCLOUD/' % (SegCount, SegCount)
					with open(CloudGEN+ImageFile) as Doc:
						for Item in Doc:
							ImagePropertiesRaw.append(Item)
					Doc.close()

					# Clean up the data, split by the asterisk *
					SplitData = ImagePropertiesRaw[0].split('*')
					# Keep all the data but the first and final blank spot left over from split
					SplitData = SplitData[0:-1]
					# Split each element in each sub list to write list of list
					ImageProperties = []
					for aItem in SplitData:
						aItem = aItem.split('^')
						ImageProperties.append(aItem)
					# New format: [Sat_PR_Date_Seg, Total Int, Total area, Water Area, Land Area]
					for aList in ImageProperties:
						# Split the data by _ to split the sat code and check the date
						Temp = aList[0].split('_')
						# TD is the temporary date, if it matches the one the loop
						# is on, it will select the clip meta data for this set
						TD = Temp[-1]
						if TD==Dates[Index]:
							# Identify it as an item
							SelectedData = aList

					# Calculate area left after clip, we only kept water and land
					TotalArea = float(SelectedData[-1])+float(SelectedData[-2])+float(ClipKeep)


					# CLIP NDVI Image properties statements: *****************
					# Build our second section of the text input for ax2
					# We are no longer removing saturation. So we need to calculate
					# what is being left. 
					###################### Clip not removed
					PercOfTotal = float(ClipKeep)/TotalAreaOfClip
					PercCorrected = round(PercOfTotal*100, 2)
					RemovedD = str(float(SelectedData[1])-PercCorrected)

					#############################
					SEG_INT = 'CLIP_DATA_INTERFERENCE: %s%%' % SelectedData[1]
					W_Area = 'FMASK_WATER_EST: %s^2m' % SelectedData[3]
					C_Area = 'FMASK_CLEAR/LAND_EST: %s^2m' % SelectedData[4]
					InitialArea = 'INITIAL_DATA_AREA: %s^2m' % SelectedData[2]
					RemovedData = 'REMOVED_DATA_INTERFERENCE: %s%%' % RemovedD
					RemainA = 'REMAINING_DATA_AREA: %i^2m' % TotalArea
					EstHabitatArea = 'ESTIMATED_HABITAT: %f^2m' % HabArea

					HabitatMetaData.append([aSeg, Dates[Index], SelectedData[1], SelectedData[3], SelectedData[4], RemovedD, HabArea])

					Cinput = 'NDVI Image Properties:\n'+SEG_INT+'\n'+W_Area+'\n'+C_Area+'\n'+InitialArea+'\n'+RemovedData+'\n'+RemainA+'\n'+EstHabitatArea

					# Plot the data over a subplot2grid setup
					gridsize = (5,2)
					fig = plt.figure(figsize=(12,9))
					ax1 = plt.subplot2grid(gridsize, (0, 0), colspan= 2, rowspan= 3)
					ax2 = plt.subplot2grid(gridsize, (3, 0), colspan = 1, rowspan=2)
					ax3 = plt.subplot2grid(gridsize, (3, 1), colspan = 1, rowspan = 2)

					# Plot the raster image on ax1
					Input = 'Satellite: '+SatCodes[Index]+'|| Path & Row: '+PR+ '|| Date' + Dates[Index]

					ax1.set_title(('NDVI Segment %s\n%s' % (SegCount, Input)), fontsize='11')
					# out_image_c = np.rollaxis(ndvi[0], 0, 1)
					# Convert data and drop bad data
					ndviA[ndviA == -0.9999] = np.nan
					out_image_c = np.rollaxis(ndviA, 0, 1)
					out_image_c = ax1.imshow(out_image_c, cmap='plasma', vmin=-1.0, vmax=1.0)

					# out_image_c = ax1.imshow(ndviIm, cmap='plasma', vmin=-1.0, vmax=1.0)
					# Modify color bar to fit size of graph
					aspect = 20
					pad_fraction = 0.5
					divider = make_axes_locatable(ax1)
					width = axes_size.AxesY(ax1, aspect=1./aspect)
					pad = axes_size.Fraction(pad_fraction, width)
					cax = divider.append_axes("right", size=width, pad=pad)

					# Plot the colorbar and set ticklabels
					plt.colorbar(out_image_c, cax=cax)
					# Insert sat, path/row, and date information
					# Input = SatCodes[Index] + '\n' + PR + '\n' + Dates[Index]
					# plt.text(1,1.025, Input, fontsize='7')

					# Plot the image information in ax2
					# textstr = Finput+'\n'+'\n'+Cinput
					textstr = Cinput
					props = dict(boxstyle='round', facecolor='wheat')
					# Hide the empty graph shhh...
					ax2.spines['bottom'].set_color('white')
					ax2.spines['top'].set_color('white') 
					ax2.spines['right'].set_color('white')
					ax2.spines['left'].set_color('white')
					ax2.tick_params(axis='x', colors='white')
					ax2.tick_params(axis='y', colors='white')
					# Plot the text
					ax2.text(0, 1, textstr, fontsize=10,
						verticalalignment='top', bbox=props)

					# Plot the histogram
					ax3.hist(ndviA, bins=BinSelect)
					# ax3.hist(ndviA)
					ax3.set_title('NDVI Value Histogram', fontsize=9)
					ax3.set_xlabel('NDVI Value')
					ax3.set_ylabel('Frequency')
					textstr = 'Max NDVI Value: %f\nMin NDVI Value: %f' % (MaxVal, MinVal)
					ax3.text(0.01,0.98, textstr, fontsize=6, verticalalignment='top', transform=ax3.transAxes)

					plt.tight_layout(pad=4, w_pad=4, h_pad=4)

					plt.savefig((NDVIPNG+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+SegCount+'NDVI_Fmasked.png'), dpi=None, 
					    facecolor='w', edgecolor='b', orientation='portrait', papertype=None, 
					    format=None, transparent=False, bbox_inches=None, pad_inches=0.1, 
					    frameon=None)

					Index  = Index+1
					# plt.show()
					plt.close()
					# sys.exit()
		except Exception as e:
			print e
			print '\n'
			print 'Continue....'
			ErrorLog.append([SegCount, e])

	print ErrorLog

	# Write file to explain errors
	with open('Analysis_Files/NDVI_Error.txt', 'w') as Doc:
		for Exc in ErrorLog:
			Doc.write(Exc)

	df = pd.DataFrame(HabitatMetaData)
	# Header = Segment, date, interference, water area, land area, removed percent, total area left 
	Header = ['Segment', 'Date', 'Interference', "Water_Area", 'Land_Area', 'Perc_Removed', 'Habitat_Area']
	df.to_csv('Analysis_Files/HabitatEstimate.csv', header=Header)


# Write a single habitat file
# Identify all habitat files for given days

# WORK ON THIS LATER


# ShapeFileDF = pd.DataFrame(columns=['id', 'area', 'geometry'])
# with fiona.open(File, "r") as src:
# 	ESPG = src.crs
# 	for Polygon in src:
# 		# Define the id for a given poly
# 		ID = int(Polygon['properties']['val'])
# 		if ID == 1:
# 			# Define geom
# 			Geom = Polygon['geometry']['coordinates'][0]
# 			Geom = tuple(Geom)
# 			# make a polygon
# 			from shapely.geometry import Polygon
# 			Po = Polygon(Geom)
# 			# calculate area (meters)
# 			Area = Po.area
# 			AreaSum = AreaSum + Area
# 			# Define insert item
# 			InsertItem = [ID, Area, Po]
# 			# Start a df so you can transpose it to the
# 			# correct format
# 			TempDF = pd.DataFrame(InsertItem)
# 			# Rotate the data to match our desired format
# 			TempDF = TempDF.transpose()
# 			# Insert the column names
# 			TempDF.columns = ['id', 'area', 'geometry']
# 			# Add to main df
# 			ShapeFileDF = ShapeFileDF.append(TempDF, ignore_index = True)

# src = None
# gdf = gpd.GeoDataFrame(ShapeFileDF, geometry='geometry')
# gdf.crs = ESPG
# gdf.to_file(NDVISHP+mainDat+"Habitat.shp")
# gdf = None


