
# INFORMATION HERE BRUV!!!

import Modules
from Modules import *

CurrentLoc = os.path.abspath(__file__)
SplitCL = CurrentLoc.split('\\')
DesiredInfo = SplitCL[0:-2]

CompletePath = ''
i = 0
for aItem in DesiredInfo:
	if i>0:
		CompletePath = CompletePath+'\\'+aItem
	if i==0:
		CompletePath = aItem
	i = i+1


CSV_LOC = CompletePath+'/Nest_CSV/'

# Write a list of available files in Nest_CSV. This data 
# is not downloaded automatically, manual data download is 
# required. File name is not relevant. The information in the files are
# compared to deterine if the data is entered twice, or if an updated file 
# with more data has been added (useful mid season).

# ADD CHECK FOR SAME LISTS AND WRONG FILES!!!!!!*******************
# SAYS IT ABOVE BUT IT IS NOT YET IN HERE DON'T FORGET! 4/26/20

##################################### PULL DATA FROM CSV #############################
# Open the csv files in the folder and correct and extract the data into a list
# of lists.

Files = listdir(CSV_LOC)

# ****
################################### This method is not used at the moment!!!!!
##################################### 9/22/2020
############################################ Could be more stable
# This method fails in certain cases due to complications with added text from
# the field. May be more stable if decide to work out the issues, will pass them
# at the moment, but pandas is simpler to use.

# Data = []

# Header = '"Nest Id","Species","Year","Segment","River Mile","Site Name","X Coord","Y Coord","Z Coord","Loc Qual","Init Date","Fate","Fate Cause","Known Hch","Calc Hatch","Clutch Sz","Egg Hatch","Egg Unhch","Egg Unkn","Eggs Lost Early","Incidental","At Risk","Moved","Caged","Raised","N Comment"\n'

# for aFile in Files:
# 	# Open the csv file individually.
# 	with open(CSV_LOC+aFile, 'r') as Extract:
# 		Temp = []
# 		for aRow in Extract:
# 			# Identify the header, we need the format to stay the same.
# 			# Also, we do not need the header in the data.
# 			if str(aRow) == Header:
# 				pass
# 			if str(aRow) != Header:
# 				Temp.append(aRow)
# 				# Replace the empty data with Nan
# 				Mod = aRow.replace('""', 'Nan')
# 				# Drop the extra "" on objects
# 				Mod = Mod.replace('"', '')
# 				# Split the incoming data from string to list object
# 				Mod = Mod.split(',')
# 				# We only need a few sections of this data, select those 
# 				# objects.
# 				# Format: [Nest Id, Species Id, Year, X, Y, Z, Init Date, Calc Hatch Date]
# 				# Also convert these date items to datetime objects for anaysis and modification
# 				# Try dates first, if not present, insert null data
# 				# print Mod
# 				try:
# 					D1 = datetime.datetime.strptime(Mod[10], '%m/%d/%Y')
# 				except:
# 					# This means the D1 conversion failed
# 					D1 = 'Nan'
# 				try:
# 					D2 = datetime.datetime.strptime(Mod[14], '%m/%d/%Y')
# 				except:
# 					# This means the D2 conversion failed
# 					D2 = 'Nan'

# 				try:
# 					InsertItem = [Mod[0], Mod[1], int(Mod[2]), float(Mod[6]), float(Mod[7]), float(Mod[8]), D1, D2]
# 					Temp.append(InsertItem)
# 				except:
# 					pass
# 					# Assume bad data, can happen from comment format
# 	# Close out information and extract the data to a Data list
# 	Extract.close()
# 	Data.append(Temp)
# 	Temp = []

#****

# FORMAT OF CSV DATA (26 items each row)
# Nest Id, Species, Year, Segment, River Mile, Site Name,	
# X Coord, Y Coord, Z Coord, Loc Qual, Init Date, Fate,	
# Fate Cause, Known Hch, Calc Hatch, Clutch Sz,	Egg Hatch,
# Egg Unhch, Egg Unkn, Eggs Lost Early,	Incidental,
# At Risk, Moved, Caged, Raised, N Comment

# FORMAT OF Selected DATA (8 items each row)
# Format: [Nest Id, Species Id, Year, X, Y, Z, Init Date, Calc Hatch Date]

###################### GEOREFERENCE THE DATA USING PANDAS #######################

# Take lists and convert them to a pandas df

# Main dataframe
DF_M = pd.DataFrame(columns = ['Nest Id', 'Species', 'Year', 'X Coord', 'Y Coord', 'Z Coord', 'Init Date', 'Calc Hatch'])

# Pandas option

for aFile in Files:
	df = pd.read_csv(CSV_LOC+aFile)
	df = df.drop(['Segment', 'River Mile', 'Site Name', 'Loc Qual',
		'Fate', 'Fate Cause', 'Known Hch', 'Clutch Sz', 'Egg Hatch',
		'Egg Unhch', 'Egg Unkn', 'Eggs Lost Early', 'Incidental', 
		'At Risk', 'Moved', 'Caged', 'Raised', 'N Comment'], axis = 1)
	# print df
	# print df.head()
	DF_M = DF_M.append(df, ignore_index = True)
# Rename the columns to remove spaces
DF_M = DF_M.rename(columns = {'X Coord': 'X','Y Coord': 'Y','Z Coord': 'Z'})


# points from xy: Longitude, Latitude
try:
	GDF_M = gpd.GeoDataFrame(DF_M, geometry = gpd.points_from_xy(DF_M.X, DF_M.Y))
except:
	print 'GPD failure... Alternative Solution'
	# Use with fiona driver error here, but works in place of gpd if fails
	geometry = [Point(xy) for xy in zip(DF_M.X, DF_M.Y)]
	DF_M = DF_M.drop(['X', 'Y'], axis=1)
	GDF_M = GeoDataFrame(DF_M, geometry = geometry)	
	# https://gis.stackexchange.com/questions/174159/convert-a-pandas-dataframe-to-a-geodataframe



# Write it to standard WGS 84, this is what the data is in, it will
# be converted when opened next.
GDF_M.crs = {'init' : 'epsg:4326'}

# Write to a shapefile 
GDF_M.to_file('Shapefile/Nests.shp', driver = 'ESRI Shapefile')







