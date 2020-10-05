import sys
sys.path.append("TrainingData/PythonFiles")

import Modules
from Modules import*

### Our files that have just been processed are in the Parsed SegmentData folder.
# Here we can further proc this data. 

SegsRan = os.listdir('ParsedSegmentData/')

# List the folders in each seg

Folders = []
for aSegFol in SegsRan:
	Info = os.listdir('ParsedSegmentData/'+aSegFol)
	Folders.append(Info)

### 
# dir = 'MasterData/TrainData/Combinations'
# if not os.path.exists(dir):
#     os.makedirs(dir)
# MainOut = 'MasterData/TrainData/Combinations'

# Pull all habitat shapefiles and cloud data for figures
DatesOfHabitat = []
DatesOfCloud = []
Years = []
HabitatSHP = []
CloudFilesOnly = []
for FileFold in Folders:
	# Define the segment
	Seg = FileFold[0][3:-5]
	# Define the NDVI files, this is where the habitat 
	# data is.
	Loc = ('ParsedSegmentData/'+'Segment%s/'+FileFold[1]+'/SHP_Files') % Seg
	NDVIFILES = os.listdir(Loc)
	LocCloud = ('ParsedSegmentData/'+'Segment%s/'+'Seg'+Seg+'Cloud/SHP_Files') % Seg
	CLOUDFILES = os.listdir(LocCloud)
	# Pull date data from NDVI files for compiling all dates into single
	# file
	for aFile in NDVIFILES:
		Sp1 = aFile.split('_')
		DC = Sp1[2]
		Years.append(Sp1[2][0:4])
		if Sp1[-1]=='Habitat.shp':
			# Only one shapefile per date. Multiple supporting files.
			# This makes it to where only one date is pulled
			DatesOfHabitat.append(DC)
			# add file to list for combining later
			HabitatSHP.append(Loc+'/'+aFile)
	for aFile in CLOUDFILES:
		try:
			Sp1 = aFile.split('_')
			DC = Sp1[2]
			if Sp1[-1]=='qa.shp':
				DatesOfCloud.append(DC)
				CloudFilesOnly.append(LocCloud+'/'+aFile)
		except:
			pass


YearFilt = collections.Counter(Years).keys()
DatesFilt = collections.Counter(DatesOfHabitat).keys()

CDatesFilt = collections.Counter(DatesOfCloud).keys()
# # Write an output folder specific to the year of interest
# dir = '%s_FullShapes' % YearFilt[0]
# if not os.path.exists(dir):
#     os.makedirs(dir)
# MainOut = '%s_FullShapes/' % YearFilt[0]

# # Filter through the dates and combine appropriate shapefiles
# # Combine all NDVI/Habitat data
# for aDate in DatesFilt:
# 	FilesOfInterest = []
# 	for aFile in HabitatSHP:
# 		Sp1 = aFile.split('_')
# 		DateCode = Sp1[-3]
# 		if DateCode==aDate:
# 			FilesOfInterest.append(aFile)

# 	# Open the shapefiles and append them to a df
# 	ShapeFileDF = pd.DataFrame(columns=['id', 'area', 'geometry'])
# 	for aFile in FilesOfInterest:
# 		Sp1 = aFile.split('_')
# 		SegID = Sp1[-2]
# 		with fiona.open(aFile, "r") as src:
# 			ESPG = src.crs
# 			for Polygon in src:

# 				# Define the id for a given poly
# 				# ID = int(Polygon['properties']['id'])
# 				# Omit the data we do not want to clip. Such as land, water, saturation
# 				# and the image edge.
# 				Geom = Polygon['geometry']['coordinates'][0]
# 				Geom = tuple(Geom)
# 				# make a polygon
# 				from shapely.geometry import Polygon
# 				Po = Polygon(Geom)
# 				# calculate area (meters)
# 				Area = Po.area
# 				# Define insert item
# 				InsertItem = [SegID, Area, Po]
# 				# Start a df so you can transpose it to the
# 				# correct format
# 				TempDF = pd.DataFrame(InsertItem)
# 				# Rotate the data to match our desired format
# 				TempDF = TempDF.transpose()
# 				# Insert the column names
# 				TempDF.columns = ['id', 'area', 'geometry']
# 				# Add to main df
# 				ShapeFileDF = ShapeFileDF.append(TempDF, ignore_index = True)


# 		src = None
# 	# Write shapefile
# 	gdf = gpd.GeoDataFrame(ShapeFileDF, geometry='geometry')
# 	gdf.crs = ESPG
# 	gdf.to_file(MainOut+"%s_HabitatData.shp" % aDate)
# 	gdf = None
# 	ShapeFileDF = None

# # Combine all Cloud data
# for aDate in CDatesFilt:
# 	FilesOfInterest = []
# 	for aFile in CloudFilesOnly:
# 		Sp1 = aFile.split('_')
# 		DateCode = Sp1[-3]
# 		if DateCode==aDate:
# 			FilesOfInterest.append(aFile)

# 	# Open the shapefiles and append them to a df
# 	ShapeFileDF = pd.DataFrame(columns=['id', 'Value', 'area', 'geometry'])
# 	for aFile in FilesOfInterest:
# 		Sp1 = aFile.split('_')
# 		SegID = Sp1[-2]
# 		with fiona.open(aFile, "r") as src:
# 			ESPG = src.crs
# 			for Polygon in src:
# 				Val = Polygon['properties']['Value']
# 				# Define the id for a given poly
# 				# ID = int(Polygon['properties']['id'])
# 				# Omit the data we do not want to clip. Such as land, water, saturation
# 				# and the image edge.
# 				Geom = Polygon['geometry']['coordinates'][0]
# 				Geom = tuple(Geom)
# 				# make a polygon
# 				from shapely.geometry import Polygon
# 				Po = Polygon(Geom)
# 				# calculate area (meters)
# 				Area = Po.area
# 				# Define insert item
# 				InsertItem = [SegID, Val, Area, Po]
# 				# Start a df so you can transpose it to the
# 				# correct format
# 				TempDF = pd.DataFrame(InsertItem)
# 				# Rotate the data to match our desired format
# 				TempDF = TempDF.transpose()
# 				# Insert the column names
# 				TempDF.columns = ['id', 'Val', 'area', 'geometry']
# 				# Add to main df
# 				ShapeFileDF = ShapeFileDF.append(TempDF, ignore_index = True)


# 		src = None
# 	# Write shapefile
# 	gdf = gpd.GeoDataFrame(ShapeFileDF, geometry='geometry')
# 	gdf.crs = ESPG
# 	gdf.to_file(MainOut+"%s_IntData.shp" % aDate)
# 	gdf = None

############################## Write final figures for complete output
# return StringExp, BitVal
def ProcBit(Value):
	BitVal = '{0:016b}'.format(Value)
	# Define the basic meanings
	Fill = BitVal[-1]
	TerrainOcclusion = BitVal[-2]
	RadSat = BitVal[-4:-2]
	Cloud = BitVal[-5]
	CloudConfidence = BitVal[-7:-5]
	CloudShadowCon = BitVal[-9:-7]
	SnowIce = BitVal[-11:-9]
	CirrusCon = BitVal[-13:-11]

	BinaryValues = [Fill, TerrainOcclusion, RadSat, Cloud, CloudConfidence, CloudShadowCon, SnowIce, CirrusCon]
	######################### Explain Binary return
	# Build the associated defs for the binary values
	# Check fill value
	if Fill == '0':
		FillExp = 'This Condition Does Not Exist'
	if Fill == '1':
		FillExp = 'This Condition Exists'
	# Check for TerrainOcclusion
	if TerrainOcclusion=='0':
		TOExp = 'This Condition Does Not Exist'
	if TerrainOcclusion == '1':
		TOExp = 'This Condition Exists'		
	# Check for radiomet saturation
	if RadSat == '00':
		RSExp = 'No Bands Contain Saturation'
	if RadSat == '01':
		RSExp = '1-2 Bands Contain Saturation'
	if RadSat == '10':
		RSExp = '3-4 Bands Contain Saturation'
	if RadSat == '11':
		RSExp = '5 or More Bands Contain Saturation'
	# Check for cloud existance
	if Cloud == '0':
		CExp = 'This Condition Does Not Exist'
	if Cloud =='1':
		CExp = 'This Condition Exists'
	# Define cloud confidence
	if CloudConfidence == '00':
		CCExp = 'This Condition Does Not Exist'
	if CloudConfidence == '01':
		CCExp = 'Low Probability: Algorithm Has Low Confidence\nThat This Condition Exists (0-33 Percent)'
	if CloudConfidence == '10':
		CCExp = 'Medium Probability: Algorithm Has Medium Confidence\nThat This Condition Exists (34-66 Percent)'
	if CloudConfidence == '11':
		CCExp = 'High Probability: Algorithm Has High Confidence\nThat This Condition Exists (67-100 Percent)'
	# Define cloud shadow conf
	if CloudShadowCon == '00':
		CSExp = 'This Condition Does Not Exist'
	if CloudShadowCon == '01':
		CSExp = 'Low Probability: Algorithm Has Low Confidence\nThat This Condition Exists (0-33 Percent)'
	if CloudShadowCon == '10':
		CSExp = 'Medium Probability: Algorithm Has Medium Confidence\nThat This Condition Exists (34-66 Percent)'
	if CloudShadowCon == '11':
		CSExp = 'High Probability: Algorithm Has High Confidence\nThat This Condition Exists (67-100 Percent)'
	# Define Snow Ice confidence
	if SnowIce == '00':
		SIExp = 'This Condition Does Not Exist'
	if SnowIce == '01':
		SIExp = 'Low Probability: Algorithm Has Low Confidence\nThat This Condition Exists (0-33 Percent)'
	if SnowIce == '10':
		SIExp = 'Medium Probability: Algorithm Has Medium Confidence\nThat This Condition Exists (34-66 Percent)'
	if SnowIce == '11':
		SIExp = 'High Probability: Algorithm Has High Confidence\nThat This Condition Exists (67-100 Percent)'
	# Define cirrus con
	if CirrusCon == '00':
		CRCExp = 'This Condition Does Not Exist'
	if CirrusCon == '01':
		CRCExp = 'Low Probability: Algorithm Has Low Confidence\nThat This Condition Exists (0-33 Percent)'
	if CirrusCon == '10':
		CRCExp = 'Medium Probability: Algorithm Has Medium\nConfidence That This Condition Exists (34-66 Percent)'
	if CirrusCon == '11':
		CRCExp = 'High Probability: Algorithm Has High\nConfidence That This Condition Exists (67-100 Percent)'
	############## Defs define for binary values
	Exps = [FillExp, TOExp, RSExp, CExp, CCExp, CSExp, SIExp, CRCExp]
	TestedVars = ['Fill', 'Terrain Occlusion', 'Radiometric Saturation', 'Cloud', 'Cloud Confidence', 'Cloud Shadow Confidence', 'Snow/Ice Confidence', 'Cirrus Confidence']
	# Build the string for this specific value
	# Check which are present
	FilteredExps = []
	FilteredTV = []
	i = 0
	while i < len(BinaryValues):
		if int(BinaryValues[i])>0:
			FilteredExps.append(Exps[i])
			FilteredTV.append(TestedVars[i])
		i = i+1


	StringExp = ''
	i = 0
	while i < len(FilteredExps):
		if i == 0:
			# First loop
			StringExp = FilteredTV[i]+': '+FilteredExps[i]+'\n'
		if i>0:
			StringExp = StringExp + '\n'+ FilteredTV[i]+': '+FilteredExps[i]+'\n'
		i = i+1
	# Return the StringExp and Binary result
	return StringExp, BitVal


# This will be done for each segment.
# Only on dates habitat was computed.

for aFile in HabitatSHP:
	# Define the segment
	Sp1 = aFile.split('_')
	SegmentHabitat = Sp1[-2]
	DateCode = Sp1[-3]
	# Select the appropriate cloud segment file
	for CloudFile in CloudFilesOnly:
		Sp1 = CloudFile.split('_')
		SegmentCloud = Sp1[-2][0:-5]
		DateCodeCloud = Sp1[-3]
		if DateCodeCloud == DateCode:
			if SegmentCloud==SegmentHabitat:
				CorrectCloudFile = CloudFile
				break
	# The correct Habitat and Cloud file are now selected. Now the files
	# need to be pulled for the output figure. Interferences will be defined
	# here as well using the ProcBit def. This def converts the int values
	# to 16-bit which have associated defs for the interferences.

	Locs = gpd.read_file(aFile)
	Locs.plot()
	plt.show()
	sys.exit()


	# Cinput = 'NDVI Image Properties:\n'+SEG_INT+'\n'+W_Area+'\n'+C_Area+'\n'+InitialArea+'\n'+RemovedData+'\n'+RemainA+'\n'+EstHabitatArea

	# Plot the data over a subplot2grid setup
	# gridsize = (5,2)
	# fig = plt.figure(figsize=(12,9))
	# ax1 = plt.subplot2grid(gridsize, (0, 0), colspan= 2, rowspan= 3)
	# ax2 = plt.subplot2grid(gridsize, (3, 0), colspan = 1, rowspan=2)
	# ax3 = plt.subplot2grid(gridsize, (3, 1), colspan = 1, rowspan = 2)

	# # Plot the raster image on ax1
	# # Input = 'Satellite: '+SatCodes[Index]+'|| Path & Row: '+PR+ '|| Date' + Dates[Index]

	# # ax1.set_title(('NDVI Segment %s\n%s' % (SegCount, Input)), fontsize='11')
	# # out_image_c = np.rollaxis(ndvi[0], 0, 1)
	# # Convert data and drop bad data
	# ndviA[ndviA == -0.9999] = np.nan
	# out_image_c = np.rollaxis(ndviA, 0, 1)
	# out_image_c = ax1.imshow(out_image_c, cmap='plasma', vmin=-1.0, vmax=1.0)

	# # out_image_c = ax1.imshow(ndviIm, cmap='plasma', vmin=-1.0, vmax=1.0)
	# # Modify color bar to fit size of graph
	# aspect = 20
	# pad_fraction = 0.5
	# divider = make_axes_locatable(ax1)
	# width = axes_size.AxesY(ax1, aspect=1./aspect)
	# pad = axes_size.Fraction(pad_fraction, width)
	# cax = divider.append_axes("right", size=width, pad=pad)

	# # Plot the colorbar and set ticklabels
	# plt.colorbar(out_image_c, cax=cax)
	# # Insert sat, path/row, and date information
	# # Input = SatCodes[Index] + '\n' + PR + '\n' + Dates[Index]
	# # plt.text(1,1.025, Input, fontsize='7')

	# # Plot the image information in ax2
	# # textstr = Finput+'\n'+'\n'+Cinput
	# textstr = Cinput
	# props = dict(boxstyle='round', facecolor='wheat')
	# # Hide the empty graph shhh...
	# ax2.spines['bottom'].set_color('white')
	# ax2.spines['top'].set_color('white') 
	# ax2.spines['right'].set_color('white')
	# ax2.spines['left'].set_color('white')
	# ax2.tick_params(axis='x', colors='white')
	# ax2.tick_params(axis='y', colors='white')
	# # Plot the text
	# ax2.text(0, 1, textstr, fontsize=10,
	# 	verticalalignment='top', bbox=props)

	# # Plot the histogram
	# ax3.hist(ndviA, bins=BinSelect)
	# # ax3.hist(ndviA)
	# ax3.set_title('NDVI Value Histogram', fontsize=9)
	# ax3.set_xlabel('NDVI Value')
	# ax3.set_ylabel('Frequency')
	# textstr = 'Max NDVI Value: %f\nMin NDVI Value: %f' % (MaxVal, MinVal)
	# ax3.text(0.01,0.98, textstr, fontsize=6, verticalalignment='top', transform=ax3.transAxes)

	# plt.tight_layout(pad=4, w_pad=4, h_pad=4)

	# # plt.savefig((NDVIPNG+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+SegCount+'NDVI_Fmasked.png'), dpi=None, 
	# #     facecolor='w', edgecolor='b', orientation='portrait', papertype=None, 
	# #     format=None, transparent=False, bbox_inches=None, pad_inches=0.1, 
	# #     frameon=None)

	# Index  = Index+1
	# # plt.show()
	# plt.close()
	# sys.exit()
	