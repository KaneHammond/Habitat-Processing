# Imports for program
import Modules
from Modules import *

# Add folder paths
sys.path.append("Shapefile")
SF_Dir = 'Shapefile/'

sys.path.append('ParsedLakeData')
Par_Dat = 'ParsedLakeData/'

#################################################################
# Selecting section
RawID = []
with open('Analysis_Files/Analysis_Meta.txt', 'r') as Doc:
	for aItem in Doc:
		RawID.append(aItem)
RawID = RawID[0].split(',')

# Filter through list items and clean them up
IDs = []
for aItem in RawID:
	if aItem[-1]=='-':
		Mod = aItem[0:-1]
		Temp = Mod.split(':')
		# Var = int(Temp[-1])
		Var = Temp[-1]
		IDs.append(Var)
	if aItem[-1]!='-':
		Temp = aItem.split(':')
		# Var = int(Temp[-1])
		Var = Temp[-1]
		IDs.append(Var)

Folders_In_Seg = list(os.listdir(Par_Dat))

# Filter segment folders based upon which segments were chosen in previous
# section.
SelectFolders = []
ShapeFolders = []
for aItem in Folders_In_Seg:
	# The name of each folder will be Segment....
	# we need the numbers after the folder name to 
	# determine its ID. Since some numbers may be larger, we use
	# the length of the word 'segment' to determine location of id variables.
	TempID = int(aItem[7::])
	for ID in IDs:
		ID = int(ID)
		if ID==TempID:
			SelectFolders.append(aItem)
			SHPfold = os.listdir(Par_Dat+aItem)
			ShapeFolders.append(Par_Dat+aItem+'/'+SHPfold[0])
# Identify shapefiles for each segment chosen.
# This is written to follow Segment write, so only one shapefile
# is supposed to be present at the moment.
ShapeFiles = []
for Folder in ShapeFolders:
	# print Folder
	SHPfl = glob(Folder+'/'+'*.shp')
	SHPfl = SHPfl[0].replace('\\', '/')
	ShapeFiles.append(SHPfl)



#######################################################################

PR = []

PreviousScene = []

# for Shape in ShapeFiles:
ProgBarLimit = len(ShapeFiles)
print '\nFiltering Scenes...\n'
for i in tqdm.tqdm(range(ProgBarLimit)):
	Shape = ShapeFiles[i]

	###################### MOD

	# WORK WORK 

	# You can use within here instead of intersect. It will ensure the possible shapes
	# chosed actually surround the polygon,

	#Read Rivermile shapefile to check coordinate system
	file1 = gpd.read_file(Shape)
	#Read WRS2 shapefile to check coordinate system
	file2 = gpd.read_file(SF_Dir+"WRS2_RiverFull.shp")
	# file2 = gpd.read_file(SF_Dir+'conus_ard_grid.shp')
	#Reproject file1 to file2's crs, assuming it's available (WANT epsg:4326)
	file1 = file1.to_crs(file2.crs)

	file1.to_file(Shape)

	file1 = None
	file2 = None

	PR1 = []
	with fiona.open(Shape, 'r') as File2: 
		for Ext in File2:
			ID = int(Ext['properties']['Id'])
			with fiona.open(SF_Dir+'WRS2_RiverFull.shp', 'r') as File:
				for feature in File:
					PRtemp = str(feature['properties']['PR'])
					Geom = feature['geometry']
					ChckExt = shape(Ext['geometry'])
					# print ChckExt.within(shape(Geom))
					if ChckExt.within(shape(Geom)) == True:
						PR1.append(PRtemp)
	File2.close()
	File.close()

	# Append this to the PR list
	PR.append([ID, PR1])

# ################################# Filter PRs ################################

# # PR provides a list of segments and their associated PR. In some cases,
# # more than one is available. The current purpose is to reduce the amount of 
# # data required to run the program. Therefore, this section will ensure the 
# # minimum number of image paths is used.

AllPRCodes = []

for selection in PR:
	for aItem in selection[-1]:
		AllPRCodes.append(aItem)

Keys = collections.Counter(AllPRCodes).keys()

# Stats = collections.Counter(AllPRCodes)


# Keys are organized in numerical order via frequency. Filter through them 
# to match segments with most common scenes. Filtering out those with multiple
# available.

Checked = []
Selection = []
for aItem in Keys:
	# aItem is the PR from the loop.
	for Seg in PR:
		# PR is the segment data with associated PRs it is within
		# chk is set to 0, if a segment considered to be 'passed'
		# is within the checked list, chk will be set to 1 to allow it 
		# to be skipped. A segment is only considered 'passed' when it 
		# has been matched with a single PR, idealy early through the loop process
		# as the first items in the Keys list are the most commonly encountered
		# image paths.
		chk = 0
		for Passed in Checked:
			if Passed==Seg[0]:
				chk = 1
		if chk==0:
			# chk is 0, so this segment has not been related to PR (a second time).
			# This time it is related, it will be by the PR values associated with
			# it. If aItem is found within PPRs (Passed PRs) from the initial
			# matching, this loop will utilize aItem as the most suitable PR for 
			# the segment as it is the most common throught the segments. Reducing
			# the total number of images needed to cover the study area.
			PPRs = Seg[-1]
			for val in PPRs:
				if val==aItem:
					Checked.append(Seg[0])
					Selection.append([Seg[0], aItem])



# Append the master file to illustrate which PR selection 
# goes with each segment

with open('Analysis_Files/Analysis_Meta.txt', 'a') as Doc:
	for Item in Selection:
		if Item != Selection[-1]:
			Doc.write('PR:%s.%s,' % (Item[0], Item[1]))
		if Item == Selection[-1]:
			Doc.write('PR:%s.%s-' % (Item[0], Item[1]))

