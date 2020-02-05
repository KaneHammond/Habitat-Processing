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

# for Shape in ShapeFiles:
ProgBarLimit = len(ShapeFiles)
print '\nSelecting Scenes...\n'
for i in tqdm.tqdm(range(ProgBarLimit)):
	Shape = ShapeFiles[i]
	# Identify the coordinate system of teh WRS2 data set
	wrs = gpd.GeoDataFrame.from_file(SF_Dir+'WRS2_descending.shp')

	#Read Rivermile shapefile to check coordinate system
	file1 = gpd.read_file(Shape)

	#Open with fiona to find bounds
	c = fiona.open(Shape)
	# This is the actual crs and we must define it before changing it
	# c.crs = fiona.crs.from_epsg(26914)

	#Identify bounds of Rivermile selection (Max x and y)
	bounds = c.bounds

	#Read WRS2 shapefile to check coordinate system
	file2 = gpd.read_file(SF_Dir+"WRS2_descending.shp")

	#Reproject file1 to file2's crs, assuming it's available
	file1 = file1.to_crs(file2.crs)

	wrs_intersection = wrs[wrs.intersects(file1.geometry[0])]

	paths, rows = wrs_intersection['PATH'].values, wrs_intersection['ROW'].values

	# #######################################################################

	# WILL SELECT SINGLE IMAGE PATH/ROW

	# # Write list for max and min x y for paths
	MxMn = []

	# Iterate through each Polygon of paths and rows intersecting the area
	for i, row in wrs_intersection.iterrows():
	    # Convert row to polygon, geometry attribut contains all
	    # coordinates
	    poly = Polygon(row.geometry)
	    # Append max and min to list
	    MxMn.append(poly.bounds)

	# identify max and min xy for the target shapefile. This will be used to determine if 
	# multiple images are required.

	SelectionMXMN = file1.bounds

	# Identify which scenes are needed. If multiple are
	# selected, filter through to check which is best fit.

	if len(MxMn)>1:
		# In list, Max x is most negative; Min x is least negative;
		# Max y is greatest number, Min y is second largest. (Logic valid for 
		# North America)

		# Convert bounds of Rivermile to list for order
		# temp = file3.bounds
		temp = file1.bounds
		List = temp.values.tolist()

		# Format (Max x, Min x, Min y, Max y)
		ListRMSort = sorted(List[0])

		# Loop through scenes to see which one fully encompasses area.
		# Also select which image is most centered on RM.

		# Differences calculates the difference between max and min 
		# coordinates between the scene and Rivermile selection
		Differences = []
		temp = []

		for aItem in MxMn:
			# Order coordinates from least to greatest 
			Order = sorted(aItem)
			i = 0
			# Loop calculation for each scene
			for var in Order:
				# Calculate the difference between variables
				# ListRMSort is the iteration of max and min
				# for the Rivermile selection.
				x = var-ListRMSort[i]
				i = i+1
				temp.append(x)
			Differences.append(temp)
			temp = []
		# Relationship between Differences and shape:
		# Valid if x for each var is => ( - , + , - , + )
		# If multiple scenes are valid, the one containing largest
		# valid difference is selected as it is most centered.

		Valid = []
		# 1 is no
		# 2 is yes
		temp = []
		for aList in Differences:
			i = 0
			temp = []
			for aItem in aList:
				# Identify if each variable for [0] and [2] are correct
				if i == 0 or i == 2:
					if aItem > 0:
						temp.append(1)
					if aItem < 0:
						temp.append(2)
				# Identify if each variable for [1] and [3] are correct
				if i == 1 or i == 3:
					if aItem < 0:
						temp.append(1)
					if aItem > 0:
						temp.append(2)
				i = i+1

			# Determine which scenes are valid
			# Use counter to determine keys (which numbers are present)
			y = collections.Counter(temp)
			Variables = y.keys()
			# If length is 1, all values are either correct or incorrect
			if len(Variables)==1:
				for aItem in Variables:
					if aItem == 1:
						Valid.append(1)
					if aItem == 2:
						Valid.append(2)
			# If lenght is greater than 1, then not all are true. 
			# Therefore, the set is not sufficient for the whole
			# study area.
			if len(Variables) > 1:
				Valid.append(1)
			temp = []

		# If value in Valid are 1, the scene does not fully cover the area.
		# If value in Valid is 2, the scene covers the study area entirely.

		# Select which scene covers area. If both fail, then both will be used.
		i = 0
		SelectScene = []
		for aItem in Valid:
			if aItem==2:
				SelectScene.append(i)
			if aItem!=2:
				pass
			i = i+1

		# If multiple pass, select which is most centered.

		# Differences between Max and Min values to Rivermile
		Mx = []
		My = []
		Mnx = []
		Mny = []
		i = 0
		temp = []
		# Which scene has highest distance for that variable
		HighestDist = []

		if len(SelectScene) > 1:
			# Add separate variables for diff to
			# individual list for ordering
			for aList in Differences: 
				Mx.append(aList[0])
				Mnx.append(aList[1])
				Mny.append(aList[2])
				My.append(aList[3])
			# Sort the lists in order
			# Most negative is furthest
			Mx = sorted(Mx)
			# Higher positive is further
			Mny = sorted(Mny)
			# Most negative is furthest
			Mnx = sorted(Mnx)
			# Higher positive is further
			My = sorted(My)

			for aList in Differences: 
				# Highest Diff Max x
				if aList[0]==Mx[0]:
					temp.append(i)
				if aList[0]!=Mx[0]:
					temp.append('Null')
				# Highest Diff Min x
				if aList[1]==Mnx[0]:
					temp.append(i)
				if aList[1]!=Mnx[0]:
					temp.append('Null')
				# Highest Diff Min y
				if aList[2]==Mny[-1]:
					temp.append(i)
				if aList[2]!=Mny[-1]:
					temp.append('Null')
				# Highest Diff Min y
				if aList[3]==My[-1]:
					temp.append(i)
				if aList[3]!=My[-1]:
					temp.append('Null')
				HighestDist.append(temp)
				temp = []
				i = i+1
			# Min and Max y can be the same. This is the Latitude of image.
			# Check if a set completely satisfies coverage or if multiple 
			# images are required.
			i = 0
			Reject = []
			for aList in HighestDist:
				for aItem in aList:
					if aItem=='Null':
						Reject.append(i)
				i = i+1
			# Count which are rejected based on index
			SceneR = collections.Counter(Reject)
			RejectedScenes = SceneR.keys()

			# Check if all scenes rejected
			DownloadAll = []
			if len(RejectedScenes)==len(wrs_intersection):
				DownloadAll.append(2)

			# Index and count all scenes for selection of best 
			i = 0
			Index = []

			for aItem in range(len(wrs_intersection)):
				Index.append(i)
				i=i+1

			res = list(set(Index)-set(RejectedScenes))
	# print ('path: %03d, row: %03d' % (row.PATH, row.ROW))
	Path = '%03d' % row.PATH
	Row = '%03d' % row.ROW
	PRITEM = Path+'/'+Row
	PR.append(PRITEM)

# Append the master file to illustrate which PR selection 
# goes with each segment

i = 0
with open('Analysis_Files/Analysis_Meta.txt', 'a') as Doc:
	for ID in IDs:
		if ID != IDs[-1]:
			Doc.write('PR:%s.%s,' % (ID, PR[i]))
		if ID == IDs[-1]:
			Doc.write('PR:%s.%s-' % (ID, PR[i]))
		i = i+1
