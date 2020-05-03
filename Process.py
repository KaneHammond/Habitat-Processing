
# Imports for program
import Modules
from Modules import *
# print 'Imports Complete...'
# Define folder paths
SF_Dir = 'Shapefile/'
PLD = 'ParsedLakeData/'
Meta_Dir = 'Analysis_Files/'
ImageMetaFiles = Meta_Dir+'Image_Meta/'



###############################################################################

# Write index to determine how many segments are present, as well as IDs
# This will be based off of how many segments were chosen in the first section.
ID = []
SegmentFolders = os.listdir(PLD)

for Folder in SegmentFolders:
	ID.append(Folder[7::])

# Read Analysis_Meta.txt to pull defined path and row for each area
# within the analysis. This will be indexed by the segment number
Content = []
with open(Meta_Dir+'Analysis_Meta.txt', 'r') as Doc:
	for item in Doc:
		Content.append(item)

# Clean up the data
Content = Content[0].split('-')
Content = Content[0:-1]
# Split the path row information from the text, this section will also include the 
# id associated with each P/R pair
PR_Raw = Content[-1]
PR_Raw = PR_Raw.split(',')
PRs = []
for aItem in PR_Raw:
	aItem = aItem.replace('PR:', '')
	aItem = aItem.split('.')
	PR_Split = aItem[-1].split('/')
	PRs.append([str(aItem[0]), PR_Split[0], PR_Split[1]])

# Coverage Data list. This will be used to write a csv metafile which
# contains all coverage information for each clip and date
CoverAll = []

ProgBarLimit = len(ID)
# ProgBarLimit = 2
print '\nCloud Segment Analysis...\n'
for r in tqdm.tqdm(range(ProgBarLimit)):
	aSeg = PRs[r]
	SegCount = aSeg[0]
	# SegCount = 31
	# SegCount = 0
	SHPfl = 'Segment%s/Seg%sShp/Seg%sShp.shp' % (SegCount, SegCount, SegCount)
	SHPflFolder = 'Segment%s/Seg%sShp/' % (SegCount, SegCount)
	# Loop through the lake segments and write data
	###########################################################################
	#							Calculate Segement Cloud Coverage and NDVI
	CPR = aSeg[1]+aSeg[2]

	# Identify folder paths within image output
	LANDSAT_PATH = glob('TIF_Code_Images*/*')
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
		S = aItem.split('\\')
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
	dir = 'ParsedLakeData/Segment%s/Seg%sCLOUD' % (SegCount, SegCount)
	if not os.path.exists(dir):
		os.makedirs(dir)
	# Shapefile Directory
	dir = 'ParsedLakeData/Segment%s/Seg%sCLOUD/SHP_Files' % (SegCount, SegCount)
	if not os.path.exists(dir):
		os.makedirs(dir)
	# tif clip Directory
	dir = 'ParsedLakeData/Segment%s/Seg%sCLOUD/TIF_Files' % (SegCount, SegCount)
	if not os.path.exists(dir):
		os.makedirs(dir)
	# Tif clips/shp files
	dir = 'ParsedLakeData/Segment%s/Seg%sCLOUD/PNG_Files' % (SegCount, SegCount)
	if not os.path.exists(dir):
		os.makedirs(dir)

	# Define the output folders for each data type 
	CloudPNG = 'ParsedLakeData/Segment%s/Seg%sCLOUD/PNG_Files/' % (SegCount, SegCount)
	CloudSHP = 'ParsedLakeData/Segment%s/Seg%sCLOUD/SHP_Files/' % (SegCount, SegCount)
	CloudTIF = 'ParsedLakeData/Segment%s/Seg%sCLOUD/TIF_Files/' % (SegCount, SegCount)
	CloudGEN = 'ParsedLakeData/Segment%s/Seg%sCLOUD/' % (SegCount, SegCount)


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

		##################### DEM GRID BASED

		# # Open the shapefile and determine which grid the section
		# # belongs to. Then use the full grid to clip the image data.
		# # This will be the new extent. For this purpose, we are only
		# # looking at the phase. Expansion will require more specific
		# # filters. Currently we only have the JamesRiver data.
		# # Split between 2 phases, 4 and 5.
		# with fiona.open(PLD+SHPflFolder+'Seg'+str(SegCount)+'Shp'+'.shp', "r") as shapefile:
		# 	phase = [feature["properties"]["Phase"] for feature in shapefile]
				
		# shapefile.close()
		# # Filter through the 2 phases available for the lake location
		# if int(phase[0])==4:
		# 	GridFile = 'IndexLiDARJamesRiverPh4QL3.shp'

		# if int(phase[0])==5:
		# 	GridFile = 'IndexLiDARJamesRiverPh5QL3.shp'

		# # Open the appropriate gridded file to pull the full
		# # polygon. This will be the correct bounds for the
		# # specific DEM file.

		# with fiona.open(SF_Dir+GridFile, "r") as shapefile:
		# 	for feature in shapefile:
		# 		if str(feature['properties']['TileNum'])==SegCount:
		# 			source_schema = shapefile.schema
		# 			crs = shapefile.crs
		# 			source_properties = shapefile['properties']
		# 			# Identify the location of the bound file for the DEM segment
		# 			BoundFile = PLD+'Segment'+SegCount+'/'+'Seg'+SegCount+'Shp/'+SegCount+'Bounds.shp'
		# 			with fiona.open(PLD+'Segment'+SegCount+'/'+'Seg'+SegCount+'Shp/'+SegCount+'Bounds.shp', 'w', 
		# 				driver = shapefile.driver, schema = source_schema, crs = crs) as NewFile:
		# 				# NewFile['geometry'] = feature['geometry']
		# 				# NewFile['properties'] = source_properties
		# 				# for f in feature:
		# 				NewFile.write(feature)
		# 				NewFile.close()
		# shapefile.close()

		####################### END DEM GRID BASED

		# Write bound file based on parsed lake data

		with fiona.open(SF_Dir+'LakeParsed.shp', "r") as shapefile:
			for feature in shapefile:
				if str(feature['id'])==SegCount:
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
		ds = gdal.Open(CloudTIF+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+str(SegCount)+'Clip.tif')
		band =  ds.GetRasterBand(1)
		array = np.array(band.ReadAsArray())
		ds = None
		
		# Write lists of values from the pixel_qa array
		# Here we will extract key variables and thier occurance in the 
		# array
		values, Freq = np.unique(array, return_counts=True)
		values = values.tolist()
		Freq = Freq.tolist()
		# Drop the nodata value (0) **
		i = 0
		Loc = []
		if 0 in values:
			for aItem in values:
				if aItem==0:
					Loc.append(i)
				i = i+1
		if len(Loc)!=0:
			values = values[1::]
			Freq = Freq[1::]	

		# This will convert values to scientific notation to simplify graphing
		# this is being added due to large gaps in numbers in the data set
		def format_e(n):
			a = '%E' % n
			return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]

		# This will reformat the 2 lists. ConValues will now be in sci notation
		# ConFreq into float
		ConValues = []
		ConFreq = []
		i = 0
		while i < len(values):
			Vfloat = values[i]
			TempV = format_e(Vfloat)
			Ffloat = Freq[i]
			TempF = Ffloat
			ConValues.append(TempV)
			ConFreq.append(TempF)
			i = i+1

		# print ConValues
		# print ConFreq
		# sys.exit()

		# Calculate coverage area for each pixel type
		# Coverage will contain the [pixel value, area]
		Coverage = []
		# CovVal = []
		i = 0
		TotalArea = 0
		while i < len(values):
			Pixels = Freq[i]
			Area = Pixels*90
			TotalArea = TotalArea+Area
			InsertItem = [values[i], Area]
			Coverage.append(InsertItem)
			i = i+1

		# Calculate the coverage of each datatype. This will vary from each 
		# satellite used. We will use this section to conduct a cloud coverage
		# assesment for the chunk of data we are actually using for the image.

		for Value in Coverage:
			PercentCoverage = float(Value[-1]/float(TotalArea))
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

		for aValue in Coverage:
			if aValue[0]==Clear:
				Percent = float(Coverage[CovI][-1])*100
				PercentRound = str(round(Percent, 2))
				# Add data to temp coverage set
				CC = PercentRound
				labels.append(('Clear '+'%s'+'%%') % PercentRound)
				Matched = True
			if aValue[0]==Water:
				Percent = float(Coverage[CovI][-1])*100
				PercentRound = str(round(Percent, 2))
				# Add data to temp coverage set
				WW = PercentRound
				labels.append(('Water '+'%s'+'%%') % PercentRound)
				Matched = True
			if Matched==False:
				Percent = float(Coverage[CovI][-1])*100
				PercentRound = str(round(Percent, 2))
				# Add data to temp coverage set
				IN = PercentRound
				labels.append(('%i_INT '+'%s'+'%%') % (aValue[0], PercentRound))
				i = i+1
			Matched = False
			CovI = CovI+1
			sizes.append(aValue[-1])
		# Add the Clear land, water, and interference data to the temp dataset	
		TempCov.append(CC)
		TempCov.append(WW)
		TempCov.append(IN)
		# TempCov Format: [Segment, Landsat info, Land, Water, Interference]

		#	Data for text box on plot ax2

		# Pull meta data for the image file and write a text box explaining 
		# data quality for the entire image, and clipped segment.

		ImagePropertiesRaw = []
		ImageFile = str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index]+'_Properties.txt')
		# print ImageMetaFiles+ImageFile
		# sys.exit()
		with open(ImageMetaFiles+ImageFile) as Doc:
			for Item in Doc:
				ImagePropertiesRaw.append(Item)
		Doc.close()

		# Used * when writing the file to allow for splitting each item later
		ImagePropertiesRaw = ImagePropertiesRaw[0].split('*')
		# Restructure to list of lists and clean up string data
		# Data information:
		#	system:index - Will return sat_pr_date
		#	EARTH_SUN_DISTANCE - Will return distance to sun
		#	system:time_start - Will return time
		#	system:footprint - Will return type of data and the coordinates
		#	LEVEL1_PRODUCTION_DATE - Will return date of production
		#	SATELLITE - Will return sat 
		#	GEOMETRIC_RMSE_MODEL - Will return the RMSE model
		#	CLOUD_COVER = Give cloud cover for entire image
		#	GEOMETRIC_RMSE_MODEL_X - RMSE model for x
		#	GEOMETRIC_RMSE_MODEL_Y - RMSE model for y
		#	IMAGE_QUALITY_OLI - Image quality value
		#	ESPA_VERSION - 
		#	WRS_ROW - 
		# 	WRS_PATH - 
		#	system:assetsize - size of image
		#	LANDSAT_ID - Complete ID
		#	SENSING_TIME - date and time of image
		#	SR_APP_VERSION - Surface reflectance model
		Selection = ['LANDSAT_ID','CLOUD_COVER','SENSING_TIME','SR_APP_VERSION']
		# List for the properties to be oreded in
		ImageProperties = []
		# List to store the selected items before ordering
		Temp = []
		for aItem in ImagePropertiesRaw:
			aItem = aItem.replace(']', '')
			aItem = aItem.replace('[', '')
			aItem = aItem.replace("u'", '')
			aItem = aItem.replace("'", '')
			aItem = aItem.replace(" ", '')
			aItem = aItem.split(',')
			for selected in Selection:
				if selected==aItem[0]:
					Temp.append(aItem)
		# Order the selection by the order in the Selection input list
		for aItem in Selection:
			for selected in Temp:
				if selected[0]==aItem:
					ImageProperties.append(selected)

		# FULL Image properties statements: *****************
		Sat = ImageProperties[0][0]+':'+' '+ImageProperties[0][-1]
		Cloud = ImageProperties[1][0]+':'+' '+ImageProperties[1][-1]+'%'
		Time = ImageProperties[2][0]+':'+' '+ImageProperties[2][-1]
		SRver = ImageProperties[3][0]+':'+' '+ImageProperties[3][-1]
		Finput = 'Base Image Properties:\n'+Sat+'\n'+Cloud+'\n'+Time+'\n'+SRver
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
			if aItem[0]!=Clear and aItem[0]!=Water:
				TotalInt = TotalInt+aItem[-1]
		
		# Area string
		Area = 'CLIP_AREA: %s^2 Meters' % str(format_e(TotalArea))
		# Total calculated interference
		TotalInt = TotalInt*100
		TotalInt = str(round(TotalInt))
		interference = 'ESTIMATED_CLIP_INTERFERENCE: %s%%' % (TotalInt)

		Cinput = '\n\nClip Image Properties:'+'\n'+interference+'\n'+Area

		# Set colors and labels for the ax1 and ax3 figures

		# Change numpy array values for clean plotting
		# Open the complete data set for the conversion. Insert an
		# integer to replace the values in the NP array with
		ReplaceInt = 1
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
			bounds.append(aItem[-1])
			Colors.append(RandomColor)


		# Copy the colors as is for the pie chart
		ColorsPie = copy.deepcopy(Colors)
		# Insert a value for nodat on the ax1 plot
		Colors.insert(0, 'white')

		# Define bounds on the ax1 plot
		Colors = LinearSegmentedColormap.from_list('Custom cmap', Colors)
		bounds = np.linspace(0, len(ColorsPie), len(ColorsPie)+2)
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
		a1L.insert(0, 'No Data')

		# Plot the colorbar and set ticklabels
		plt.colorbar(out_image_c, cax=cax, ticks=boundary_means).set_ticklabels(a1L)
		
		# Insert sat, path/row, and date information
		Input = SatCodes[Index] + '\n' + PR + '\n' + Dates[Index]
		plt.text(0,1.0025, Input, fontsize='7')

		# Plot the image information
		textstr = Finput+Cinput
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
		patches, texts = ax3.pie(sizes, startangle=90, colors=ColorsPie)

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

# Write the list as a pandas dataframe, then use the pandas 
# module to write that to a csv file.
df = pd.DataFrame(CoverAll)
df.to_csv(Meta_Dir+'SegmentImageInfo.csv')


#####################################################################################
					# Write shapefiles from the cloud raster data

CloudFolderPaths = []

for Segment in ID:
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
# Close out of files
source.close()
output.close()

######################################################################################
#									NDVI SECTION

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

ProgBarLimit = len(ID)

print '\nNDVI Segment Analysis...\n'
for r in tqdm.tqdm(range(ProgBarLimit)):
	# if r == 4:
	# r = 4
	aSeg = PRs[r]
	# print aSeg
	# aSeg = PRs[6]
	SegCount = aSeg[0]
	# print SegCount
	SHP_Extent_File = 'Segment%s/Seg%sShp/%sExtent.shp' % (SegCount, SegCount, SegCount)
	# Write output paths for cloud assessment
	# General output for cloud data
	dir = 'ParsedLakeData/Segment%s/Seg%sNDVI' % (SegCount, SegCount)
	if not os.path.exists(dir):
		os.makedirs(dir)
	# tif clip Directory
	dir = 'ParsedLakeData/Segment%s/Seg%sNDVI/TIF_Files' % (SegCount, SegCount)
	if not os.path.exists(dir):
		os.makedirs(dir)
	# PNG files
	dir = 'ParsedLakeData/Segment%s/Seg%sNDVI/PNG_Files' % (SegCount, SegCount)
	if not os.path.exists(dir):
		os.makedirs(dir)

	# Define the output folders for each data type 
	NDVIPNG = 'ParsedLakeData/Segment%s/Seg%sNDVI/PNG_Files/' % (SegCount, SegCount)
	NDVITIF = 'ParsedLakeData/Segment%s/Seg%sNDVI/TIF_Files/' % (SegCount, SegCount)

	try:

		# Calculate Segement NDVI
		CPR = aSeg[1]+aSeg[2]

		# Identify folder paths within image output
		LANDSAT_PATH = glob('TIF_Code_Images*/*')
		# List of folders containing files
		Paths = []
		# Files within each folder
		Files = []
		# Define folder paths for NDVI out for multiple dates
		NDVI_Out_Paths = []
		# Define sat codes
		SatCodes = []
		# Define dates
		Dates = []
		# PR codes
		PRCodes = []
		for aItem in LANDSAT_PATH:
			# Select the path and row for the image collection to 
			# compare to the CorrectPR object
			S = aItem.split('\\')
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
			# Determine which bands to use based upon satellite
			if SatCodes[Index] == 'LC08':
				Band1 = 'B5' # NIR
				Band2 = 'B4' # Red
			if SatCodes[Index] == 'LT05':
				Band1 = 'B4' # NIR
				Band2 = 'B3' # Red

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
				# Select the B3 and B4 or B4 and B5 bands
				# Define the imagery for the path and row of the polygon
				for aItem in glob(Path+'*/*%s.TIF' % Band1):
					OriginalData = aItem
					B1 = aItem
				for aItem in glob(Path+'*/*%s.TIF' % Band2):
					B2 = aItem
					Temp = aItem

				import rasterio.mask
				with fiona.open(PLD+SHP_Extent_File, "r") as shapefile:
					shapes = [feature["geometry"] for feature in shapefile]

				with rasterio.open(B1) as src:
					out_image_B1, out_transform1 = rasterio.mask.mask(src, shapes, crop=True)
					out_meta_B1 = src.meta
					# Define the parameters of the tif clip
					out_meta_B1.update({"driver": "GTiff",
			                 "height": out_image_B1.shape[1],
			                 "width": out_image_B1.shape[2],
			                 "transform": out_transform1})
				# Write a clipped B1 tif file
				with rasterio.open(NDVITIF+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+SegCount+'_'+Band1+'_'+'Clip.tif', 'w', **out_meta_B1) as dst:
				    dst.write(out_image_B1)
					# print out_meta_B4
					# print out_transform4
				dst.close()

				with rasterio.open(B2) as src:
					out_image_B2, out_transform = rasterio.mask.mask(src, shapes, crop=True)
					out_meta_B2 = src.meta
					# Define the parameters of the tif clip
					out_meta_B2.update({"driver": "GTiff",
			                 "height": out_image_B2.shape[1],
			                 "width": out_image_B2.shape[2],
			                 "transform": out_transform})
				# Write a clipped b4 tif file
				with rasterio.open(NDVITIF+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+SegCount+'_'+Band2+'_'+'Clip.tif', 'w', **out_meta_B2) as dst:
				    dst.write(out_image_B2)
					# show(out_image_B3)
				dst.close()

				# Open the shapefile containing no data shapes
				FolderLoc = 'ParsedLakeData/Segment%s/Seg%sCLOUD/SHP_Files/' % (SegCount, SegCount)
				File = str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+SegCount+'INT.shp'
				
				# Open the interference shapefile for clip
				with fiona.open(FolderLoc+File, "r") as shapefile:
					shapes = [feature["geometry"] for feature in shapefile]
				shapefile.close()


				# Open the clipped bands for the NDVI calculation
				# We write the new tif files and open them for the calculation to aid in 
				# writing a new tif file for the ndvi calculation. Error was lifted when 
				# trying to save the ndvi tif image from the clips of band imagery.
				# print str((SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+SegCount+'_'+Band1+'_'+'Clip.tif')
				with rasterio.open(NDVITIF+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+SegCount+'_'+Band1+'_'+'Clip.tif', 'r') as src:
					out_image_B1, out_transform1 = rasterio.mask.mask(src, shapes, crop=True)
					out_meta_B1 = src.meta

				with rasterio.open(NDVITIF+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+SegCount+'_'+Band2+'_'+'Clip.tif', 'r') as src2:
					out_image_B2, out_transform2 = rasterio.mask.mask(src2, shapes, crop=True)
					out_meta_B2 = src.meta

				# Calculate NDVI
				# Allow division by zero
				np.seterr(divide='ignore', invalid='ignore')
				ndvi_upper = (out_image_B1.astype(float) - out_image_B2.astype(float))
				ndvi_lower = (out_image_B1.astype(float) + out_image_B2.astype(float))
				ndvi = (ndvi_upper / ndvi_lower)

				# Update the meta data the new 
				out_meta_B1.update({"driver": "GTiff",
				                 "height": out_image_B1.shape[1],
				                 "width": out_image_B1.shape[2],
				                 "transform": out_transform1,
				                 # Define dytpe as float to accomodate for the dec in calc
				                 "dtype": "float64"})

				# WORK HERE DEFINE VMIN VMAX

				# Write the calculated NDVI data to a tif image
				with rasterio.open(NDVITIF+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+SegCount+'_'+'ndvi.tif', 'w', **out_meta_B1) as dst:
				    dst.write(ndvi)

				# Close out of files
				dst.close()
				src.close()
				src2.close()
				# Delete the band clips previously written
				os.remove(NDVITIF+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+str(SegCount)+'_'+Band1+'_'+'Clip.tif')
				os.remove(NDVITIF+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+str(SegCount)+'_'+Band2+'_'+'Clip.tif')

				#######################################################################
				#							NDVI PNG OUTPUT

				# Histogram bin data
				# Define output text information
				bins = [-1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.8, 0.9, 1]

				# Get the max and min variables in the data, this will be used for the 
				# the bins.
				MinVal = np.nanmin(ndvi[0])
				MaxVal = np.nanmax(ndvi[0])
				BinSelect = []

				for aValue in bins:
					if aValue>MinVal:
						if aValue<MaxVal:
							BinSelect.append(aValue)

				# Figure Text section on ax2

				# Pull meta data for the base image file 

				ImagePropertiesRaw = []
				ImageFile = str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index]+'_Properties.txt')
				with open(ImageMetaFiles+ImageFile) as Doc:
					for Item in Doc:
						ImagePropertiesRaw.append(Item)
				Doc.close()
				# Used * when writing the file to allow for splitting each item later
				ImagePropertiesRaw = ImagePropertiesRaw[0].split('*')
				# Restructure to list of lists and clean up string data
				# Data information:
				#	system:index - Will return sat_pr_date
				#	EARTH_SUN_DISTANCE - Will return distance to sun
				#	system:time_start - Will return time
				#	system:footprint - Will return type of data and the coordinates
				#	LEVEL1_PRODUCTION_DATE - Will return date of production
				#	SATELLITE - Will return sat 
				#	GEOMETRIC_RMSE_MODEL - Will return the RMSE model
				#	CLOUD_COVER = Give cloud cover for entire image
				#	GEOMETRIC_RMSE_MODEL_X - RMSE model for x
				#	GEOMETRIC_RMSE_MODEL_Y - RMSE model for y
				#	IMAGE_QUALITY_OLI - Image quality value
				#	ESPA_VERSION - 
				#	WRS_ROW - 
				# 	WRS_PATH - 
				#	system:assetsize - size of image
				#	LANDSAT_ID - Complete ID
				#	SENSING_TIME - date and time of image
				#	SR_APP_VERSION - Surface reflectance model
				Selection = ['LANDSAT_ID','CLOUD_COVER','SENSING_TIME','SR_APP_VERSION']
				# List for the properties to be oreded in
				ImageProperties = []
				# List to store the selected items before ordering
				Temp = []
				for aItem in ImagePropertiesRaw:
					aItem = aItem.replace(']', '')
					aItem = aItem.replace('[', '')
					aItem = aItem.replace("u'", '')
					aItem = aItem.replace("'", '')
					aItem = aItem.replace(" ", '')
					aItem = aItem.split(',')
					for selected in Selection:
						if selected==aItem[0]:
							Temp.append(aItem)
				# Order the selection by the order in the Selection input list
				for aItem in Selection:
					for selected in Temp:
						if selected[0]==aItem:
							ImageProperties.append(selected)

				# FULL Image properties statements: *****************
				Sat = ImageProperties[0][0]+':'+' '+ImageProperties[0][-1]
				Cloud = ImageProperties[1][0]+':'+' '+ImageProperties[1][-1]+'%'
				Time = ImageProperties[2][0]+':'+' '+ImageProperties[2][-1]
				SRver = ImageProperties[3][0]+':'+' '+ImageProperties[3][-1]
				Finput = 'Base Image Properties:\n'+Sat+'\n'+Cloud+'\n'+Time+'\n'+SRver

				# Pull data for the clipped image interference coverage
				# The image meta file is written for each cloud coverage section
				# and contains the image name and estimated interference. At this piont,
				# the interference will have been clipped.
				# The data format for the initial input for the meta file is below:
				# InsertItem = SatCode+'^'+TotalInt+'^'+str(TotalArea)+'^'+str(WaterArea)+'^'+str(LandArea)+'*'
				
				ImagePropertiesRaw = []
				ImageFile = 'Image_Meta.txt'
				CloudGEN = 'ParsedLakeData/Segment%s/Seg%sCLOUD/' % (SegCount, SegCount)
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
				TotalArea = int(SelectedData[-1])+int(SelectedData[-2])

				# CLIP NDVI Image properties statements: *****************
				# Build our second section of the text input for ax2
				SEG_INT = 'CLIP_DATA_REMOVED: %s%%' % SelectedData[1]
				W_Area = 'FMASK_WATER_EST: %s^2m' % SelectedData[3]
				C_Area = 'FMASK_CLEAR/LAND_EST: %s^2m' % SelectedData[4]
				InitialArea = 'INITIAL_DATA_AREA: %s^2m' % SelectedData[2]
				RemainA = 'REMAINING_DATA_AREA: %i^2m' % TotalArea

				Cinput = 'NDVI Image Properties:\n'+SEG_INT+'\n'+W_Area+'\n'+C_Area+'\n'+InitialArea+'\n'+RemainA

				# Plot the data over a subplot2grid setup
				gridsize = (5,2)
				fig = plt.figure(figsize=(12,9))
				ax1 = plt.subplot2grid(gridsize, (0, 0), colspan= 2, rowspan= 3)
				ax2 = plt.subplot2grid(gridsize, (3, 0), colspan = 1, rowspan=2)
				ax3 = plt.subplot2grid(gridsize, (3, 1), colspan = 1, rowspan = 2)

				# Plot the raster image on ax1
				ax1.set_title(('NDVI Segment %s' % SegCount), fontsize='14')
				out_image_c = np.rollaxis(ndvi[0], 0, 1)
				out_image_c = ax1.imshow(out_image_c, cmap='plasma', vmin=-1.0, vmax=1.0)

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
				Input = SatCodes[Index] + '\n' + PR + '\n' + Dates[Index]
				plt.text(0,1.025, Input, fontsize='7')

				# Plot the image information in ax2
				textstr = Finput+'\n'+'\n'+Cinput
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
				ax3.hist(ndvi[0], bins=BinSelect)
				ax3.set_title('NDVI Value Histogram', fontsize=9)
				ax3.set_xlabel('NDVI Value')
				ax3.set_ylabel('Frequency')
				# textstr = 'Max NDVI Value: %f\nMin NDVI Value: %f' % (MaxVal, MinVal)
				# ax3.text(0.01,0.98, textstr, fontsize=6, verticalalignment='top', transform=ax3.transAxes)

				plt.tight_layout(pad=4, w_pad=4, h_pad=4)

				plt.savefig((NDVIPNG+str(SatCodes[Index])+'_'+str(PRCodes[Index])+'_'+str(Dates[Index])+'_'+SegCount+'NDVI_Fmasked.png'), dpi=None, 
				    facecolor='w', edgecolor='b', orientation='portrait', papertype=None, 
				    format=None, transparent=False, bbox_inches=None, pad_inches=0.1, 
				    frameon=None)

				Index  = Index+1
				# plt.show()
				plt.close()
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

	#						Calculate Segment NDVI END



