# Imports for program
import Modules
from Modules import *

#########################################################################
# This section will write segment folders in prep for conducting the    #
# analysis.
#########################################################################

dir = 'ParsedLakeData'
if not os.path.exists(dir):
    os.makedirs(dir)

dir = 'Analysis_Files'
if not os.path.exists(dir):
    os.makedirs(dir)


# New location for parsed lake segments
PLD = 'ParsedLakeData/'
SF_Dir = 'Shapefile/'

# Will write the file with a coordinate system
# LakeSeg = gpd.GeoDataFrame.from_file(SF_Dir+'LakeSeg.shp')
# f = gpd.GeoDataFrame.from_file(SF_Dir+'ParsedLake.shp')
# LakeSeg = LakeSeg.to_crs(f.crs)
# LakeSeg.to_file('LakeSeg.shp')
# sys.exit()

# This will be used to select which segments are used in the analysis



# IDs = []
# with fiona.open(SF_Dir+'LakeSeg.shp') as File:
# 	for feature in File:
# 		IDs.append(feature['id'])

# For DEM based analysis
IDs = []
with fiona.open(SF_Dir+'DEM_Grid_Cut.shp') as File:
	for feature in File:
		IDs.append(feature['properties']['NAME'])

################################## QUERY

# This is a query for larger segments, not optimal for DEM grid
c = 1
print 'Available Segments:'
for aItem in IDs:
	print '%i) %s' % (c, aItem)
	c=c+1

query = raw_input('\nSelect Segmemt/s for download (Separate by comma): ')

Selection_Index = query.split(',')

ID_Selection = []
for aItem in Selection_Index:
	try:
		aItem = int(aItem)
	except:
		print 'Fail'
		sys.exit()
	# Index for lists start at 0, started at 1 in selection index
	# must subtract by 1 to equal true index of selection.
	I_Val = aItem-1
	ID_Selection.append(IDs[I_Val])

################################ QUERY END

# Full Run

# Will write individual shapefiles for each segment
with fiona.open(SF_Dir+'DEM_Grid_Cut.shp') as File:
	# with fiona.open(SF_Dir+'ParsedLake.shp') as File:
		# schema = File.schema
	meta = File.meta
	# print meta
	# sys.exit()
	for feature in File:
		ID = feature['properties']['NAME']
		for aItem in ID_Selection:
			aItem = int(aItem)
			if aItem==int(ID):
				# Write folders for all the data
				dir = 'ParsedLakeData/Segment'+str(ID)
				if not os.path.exists(dir):
					os.makedirs(dir)
				TempDir = 'ParsedLakeData/Segment'+str(ID)+'/'
				# Write the shapefile
				# with fiona.open(TempDir+str(Number)+'SegShp', 'w', crs=from_epsg(26914), driver='ESRI Shapefile', schema=schema) as output:
				with fiona.open(TempDir+'Seg'+str(ID)+'Shp', 'w', **meta) as output:
						output.write(feature)
File.close()

# Write file to specify the specific parameters for this analysis
with open('Analysis_Files/Analysis_Meta.txt', 'w') as Doc:
	for ID in ID_Selection:
		Doc.write('ID:')
		if ID != ID_Selection[-1]:
			Doc.write(ID+',')
		if ID == ID_Selection[-1]:
			Doc.write(ID)
	Doc.write('-')

##################################### Partial Run
# Use if all data is downloaded, choose first 5 records/datasets
# IDs = IDs[0:5]
# ID_Selection = IDs


# Will write individual shapefiles for each segment
with fiona.open(SF_Dir+'DEM_Grid_Cut.shp') as File:
	# with fiona.open(SF_Dir+'ParsedLake.shp') as File:
		# schema = File.schema
	meta = File.meta
	# print meta
	# sys.exit()
	for feature in File:
		ID = feature['properties']['NAME']
		for aItem in ID_Selection:
			aItem = int(aItem)
			if aItem==int(ID):
				# Write folders for all the data
				dir = 'ParsedLakeData/Segment'+str(ID)
				if not os.path.exists(dir):
					os.makedirs(dir)
				TempDir = 'ParsedLakeData/Segment'+str(ID)+'/'
				# Write the shapefile
				# with fiona.open(TempDir+str(Number)+'SegShp', 'w', crs=from_epsg(26914), driver='ESRI Shapefile', schema=schema) as output:
				with fiona.open(TempDir+'Seg'+str(ID)+'Shp', 'w', **meta) as output:
						output.write(feature)
File.close()

# Write file to specify the specific parameters for this analysis
with open('Analysis_Files/Analysis_Meta.txt', 'w') as Doc:
	for ID in ID_Selection:
		Doc.write('ID:')
		if ID != ID_Selection[-1]:
			Doc.write(ID+',')
		if ID == ID_Selection[-1]:
			Doc.write(ID)
	Doc.write('-')







