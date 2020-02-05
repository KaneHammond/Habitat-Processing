import Modules
from Modules import *

# Identify basic file paths
PLD = 'ParsedLakeData/'

SHP_Loc = 'Shapefile/'

# List folders within directory, this is the segment folder location
SegmentFolders = listdir(PLD)

# # List sub dir for NDVI data
# SubDir1 = []
# for aItem in SegmentFolders:
# 	InsertItem = listdir(PLD+aItem)
# 	SubDir1.append(InsertItem[1])
	
# # List the ndvi image folders within each segment 
# # folder (Within the sub dirs)
# NDVI_Folders = []
# Index = 0
# for aItem in SubDir1:
# 	TempDir = PLD+SegmentFolders[Index]+'/'+aItem+'/'
# 	FolderItems = listdir(TempDir)
# 	for aItem in FolderItems:
# 		# Define folder path to image
# 		PathInsert = TempDir+aItem
# 		# Define singel tif in folder
# 		FolderItem = listdir(PathInsert)
# 		# Combine the folder path with file name
# 		InsertItem = PathInsert+'/'+FolderItem[0]
# 		NDVI_Folders.append(InsertItem)
# 	Index=Index+1

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
			# print pt['properties']['Id']
			# sys.exit()
			point = shape(pt['geometry'])
			# SHP = shape(FileUse['geometry'])
			# print SHP
			# print point
			# print shape(Geom)
			if point.within(shape(Geom)):
				NSR.append([str(pt['properties']['Id']), Segment])
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

##################################### NDVI SECTION ########################

# Set for single date... Open our csv dataset to determine
# which images will be used. We are not using segments with 
# greater than 10% intereference
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

# Define the date to pull the data from
Year = 2014
Month = 5
Date = str(Year)+'-'+str(Month)
DateFilter = datetime.datetime.strptime(Date, '%Y-%m')

# Write a datetime item at the end of 
# each element to aid in filtering the data
for aItem in MetaCov:
	I_Date = str(aItem[-1])
	Date = datetime.datetime.strptime(I_Date, '%Y%m%d')
	aItem.append(Date)

# Filter images by given date and write path/row data list
Month_Year_Filter = []
PRs = []
for aItem in MetaCov:
	Mod = aItem[4].split('_')
	PRs.append(Mod[1])
	if aItem[-1].year==Year:
		if aItem[-1].month==Month:
			Month_Year_Filter.append(aItem)

# Extract unique PR values
PRs = list(set(PRs))

# Pull best fitting scenes by the date filter and PRs.
# This will filter the month of interest and select scene clips 
# with less than 10% interference. They will be counted in the end
# and the date with the highest coverage for the datasets will be used
for aItem in PRs:
	Temp = []
	# Filter the prs in date filter
	for Scene in Month_Year_Filter:
		Mod = Scene[4].split('_')
		PR=Mod[1]
		# Select by PR value
		if PR==aItem:
			Temp.append(Scene)
	# Extract date variables
	Dates = []
	for Selection in Temp:
		Dates.append(Selection[-1])
	Dates = list(set(Dates))
	# Determine which date has the clearest coverage
	Dates.sort()
	Coverage = []
	for aDate in Dates:
		Temp = []
		# Filter the selected scenes for the given date by 
		# interference
		for Scene in Month_Year_Filter:
			Segment = Scene[0]
			Interference = Scene[3]
			Date = Scene[-1]
			if Date == aDate:
				if float(Interference)<10.0:
					Temp.append(Scene)
		Cov = float(len(Temp))/float(len(TempMod))
		# Write the list item of info
		Coverage.append([Cov, Temp])

# Select the highest coverage, it will be the last list item
def SortList(L):
	return(sorted(L, key=lambda x: x[0]))
Coverage = SortList(Coverage)
# Last item will be the best fit. This selection ignores the coverage
# value and selects only the information for segment image files.
# NOTE: Not all segments are represented in this method, as they
# are overlooked due to poor data quality.
OptimalData = Coverage[-1][1]

# Define a total DF to input all nest info data for NDVI
NDVI_DF = pd.DataFrame(columns=['Segment', 'Nest_ID', 'NDVI'])

# Call the segment and nest point data for extracting NDVI values
for aItem in NSR_Filtered:
	# Idenfify the segment in the NRS Filter, each segment is defined
	# based upon it's nest content. Areas without nests are ignored.
	# This is only the data training section, so only locations with nests
	# can be used to train NDVI values. Again, not all segments will be used
	# if data quality prevents it.
	Segment = aItem[0]
	for Scene in OptimalData:
		# Filter through optimal data and match the record to the 
		# correct NRS filter item which contains nest Ids.
		SegmentI = Scene[0]
		if Segment==SegmentI:
			# NestIds are the rest of the NRS Filter item, the first item 
			# is the Segment ID
			NestIds = aItem[1::]
			# Open out nest converted file, this is in the correct
			# crs format.
			pts = gpd.read_file(SHP_Loc+'NestsConverted.shp')
			# Convert the Id variables to string to match the list used
			# for filtering nests within a given segment, this is given
			# in NRS filter.
			pts['Id'] = pts['Id'].astype('str')
			# Mask only the points within that segment. This removes all 
			# other nests and focuses only on those within the correct
			# NDVI clip bounds.
			mask = ~pts['Id'].isin(NestIds)
			pts = pts.loc[~mask]
			# pts = pts['geometry']
			# Write columns containing x and y coordinate data.
			# This specific method of data extraction requires this
			# for looping the data.
			pts['x'] = pts['geometry'].x
			pts['y'] = pts['geometry'].y
			pts = pts[['x', 'y', 'Id', 'geometry']]
			coords = [(x,y) for x, y in zip(pts.x, pts.y)]
			# Open the correct NDVI file. Define the landsat code in 
			# the accepted scene.
			LandsatCode = Scene[4]
			for Path in NDVI_Folders:
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
			# Open the filtered NDVI clip file
			src = rasterio.open(CorrectNDVIFile)

			# Loop through the selected nests and extract the NDVI 
			# values at those locations
			pts['Raster Value'] = [x[0] for x in src.sample(coords)]
			# Initially failed below. The code above was as follows:
			# pts['Raster Value'] = [x for x in src.sample(coords)]
			# It was modified to allow it to pass and complete the data
			# extraction. Leaving this here incase it needs to be used.
			# pts['Rater Value'] = probes.apply(lambda x: x['Raster Value'][0], axis=1)
			src.close()

			# Format a temporary dataframe
			Tempdf = pd.DataFrame()
			# (columns=['Segment', 'Nest_ID', 'NDVI'])
			# pts[['x', 'y', 'Id', 'geometry']]
			Tempdf['Nest_ID'] = pts['Id']
			Tempdf['NDVI'] = pts['Raster Value']
			Tempdf['Segment'] = Segment

			# append the temporary dataframe to the complete one, append
			# doesn't happen in place, must define the df as the append
			NDVI_DF = NDVI_DF.append(Tempdf, sort=False, ignore_index=True)























