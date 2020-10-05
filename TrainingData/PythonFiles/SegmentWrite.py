# Imports for program
import Modules
from Modules import *

#########################################################################
# This section will write segment folders in prep for conducting the    #
# analysis.
#########################################################################

dir = 'ParsedSegmentData'
if not os.path.exists(dir):
    os.makedirs(dir)

dir = 'Analysis_Files'
if not os.path.exists(dir):
    os.makedirs(dir)

# New location for parsed lake segments
PSD = 'ParsedSegmentData/'
SF_Dir = 'Shapefile/'

# This will be used to select which segments are used in the analysis

# Switch to full River
IDs = []
with fiona.open(SF_Dir+'RiverFull.shp') as File:
	for feature in File:
		IDs.append(feature['properties']['Id'])
File.close()
# Define available options

LakeSak = [3, 4, 5, 6, 7, 8, 9, 10, 11]
NorthDakota = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]
SouthDakota = [25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43]

SegmentSections = [LakeSak, NorthDakota, SouthDakota]
Options = ['Lake Sakakawea', 'North Dakota (River Segments)', 'South Dakota']

################## MAIN SEGS

c = 1
print 'Available Segments:'
for aItem in Options:
	print '%i) %s' % (c, aItem)
	c=c+1

query = raw_input('\nSelect Segmemt/s for download (Separate by comma) or A for all: ')

Selection_Index = query.split(',')

ID_Selection = []
for aItem in Selection_Index:
	try:
		aItem = int(aItem)
		# Index for lists start at 0, started at 1 in selection index
		# must subtract by 1 to equal true index of selection.
		I_Val = aItem-1
		for aValue in SegmentSections[I_Val]:
			ID_Selection.append(aValue)
	except:
		try:
			aItem = str(aItem)
			if aItem == 'A':
				print '\nAll Selected...'
		except:
			print 'Fail'
			sys.exit()

if aItem=='A':
	for aVar in IDs:
		ID_Selection.append(aVar)

################################ QUERY END

# Will write individual shapefiles for each segment
# with fiona.open(SF_Dir+'LakeParsed.shp') as File:
# with fiona.open(SF_Dir+'RiverFull.shp') as File:
# Use a simplified version for training
with fiona.open(SF_Dir+'RiverSimp.shp') as File:

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




