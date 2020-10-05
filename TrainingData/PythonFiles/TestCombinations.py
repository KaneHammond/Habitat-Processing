
import Modules
from Modules import*

# Define output loc
f1out = 'TrainData/Combinations/Analysis/'


# Pull in the data and prep it
# Open datasets

##########################TRAINED DATA
NDVI = []
with open('TrainData/Combinations/NDVI_Meta.txt', 'r') as Doc:
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
with open('TrainData/Combinations/Distance_Meta.txt', 'r') as Doc:
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

##################################NEST DATA
# Pull all relevant nesting data.

DistanceData = []

with open('Analysis_Files/NestDistance.csv') as File:
	for aRow in File:
		# Clean up the data before appending to item
		Sp1 = aRow.split(',')
		Val = Sp1[2].replace('\n', '')
		Nest = Sp1[1]
		try:
			# If converting to float fails, append fail
			Val = float(Val)
			DistanceData.append([Nest, Val])
		except:
			pass
# Pull the NDVI dataset from the previously written file. This file contains data
# from specific points in a raster image. Those points are nest locations.

# Open the NDVI dataset 

NDVIData = []

with open('Analysis_Files/NestNDVI.csv') as File:
	for aRow in File:
		# Clean up the data before appending to item
		Sp1 = aRow.split(',')
		Val = Sp1[2].replace('\n', '')
		Nest = Sp1[1]
		try:
			# If converting to float fails, append fail
			Val = float(Val)
			NDVIData.append([Nest, Val])
		except:
			pass
File.close()

# Remove the first record, it is the header
NDVIData = NDVIData[1::]

# Nest image info contains the nest and year data. Was previously used with
# training the data. Useful for this.
NestImageInfo = []

with open('Analysis_Files/NestImInfo.csv') as File:
	for aRow in File:
		# Clean up the data before appending to item
		Sp1 = aRow.split(',')
		Val = Sp1[2].replace('\n', '')
		Nest = Sp1[1]
		# Select year nest was recorded.
		Year = Val[-9:-5]
		# Only keeping the year from the landsat image name
		NestImageInfo.append([Nest, Year])
NestImageInfo = NestImageInfo[1::]

# Combine nest info with appropriate year
NDVIDataYear = []
DistanceDataYear = []
for aRec in NestImageInfo:
	Temp1 = []
	Temp2 = []
	# Format: nest id, year
	for Ndat in NDVIData:
		if Ndat[0]==aRec[0]:
			# Nest id matches
			# Add nest id
			Temp1.append(aRec[0])
			# Add year
			Temp1.append(aRec[1])
			# Add NDVI
			Temp1.append(Ndat[1])
	for Ddat in DistanceData:
		if Ddat[0]==aRec[0]:
			# Nest id matches
			# Add nest id
			Temp2.append(aRec[0])
			# Add year
			Temp2.append(aRec[1])
			# Add Distance
			Temp2.append(Ddat[1])

	# Append to total data lists if a record match was found
	if len(Temp1)>0:
		NDVIDataYear.append(Temp1)
	if len(Temp2)>0:
		DistanceDataYear.append(Temp2)

# Define the years of the study period. We want to test the number of nests 
# being placed within specific training scenarios.

YS = [2014, 2015, 2016, 2017, 2018, 2019]

###################### Test NDVI

YearlyAnalysisNDVI = []
TotalAnnualNestsNDVI = []

# Count how many nests are actually in a given year (with NDVI data)
for aYear in YS:
	TotalNests = 0
	for Nrec in NDVIDataYear:
		if int(Nrec[1])==aYear:
			TotalNests = TotalNests+1
	TotalAnnualNestsNDVI.append([aYear, TotalNests])

# Pull trained data which excludes each specific year.
# print NDVIDataYear[0]
# sys.exit()
ProgBarLimit = len(NDVI)
print '\nTesting NDVI Limitations...\n'
for kk in tqdm.tqdm(range(ProgBarLimit)):
# for aRecord in NDVI:
	aRecord = NDVI[kk]
	# temp will be all the data which can be tested for a single year
	Temp = []
	# aRecord is a training set for either a single year or combination.
	# Test which do not belong in that year.
	Years = aRecord[3::]
	for aYear in YS:
		FilteredNests = []
		# Select one of our analysis years.
		# Pull all nests for this year (NDVI only)
		YNDVI = []
		for Nrec in NDVIDataYear:
			if int(Nrec[1])==aYear:
				YNDVI.append(Nrec)

		# Now use this year to identify if the NDVI training record
		# does not include nest data from this year. If so, it can be used
		# to test nest data from aYear as it is not the training data.
		if aYear not in Years:
			# We can use this.
			# Test how many nests are within the data range
			MaxNDVI = float(aRecord[0])
			MinNDVI = float(aRecord[1])
			for aNest in YNDVI:
				if float(aNest[-1])<=MaxNDVI:
					if float(aNest[-1])>=MinNDVI:
						FilteredNests.append(aNest)
		if len(FilteredNests)>0:
			Temp.append(FilteredNests)

	# At this point, Temp is the years which can be tested with trained set
	# Define training set (Years)
	YearlyAnalysisNDVI.append([Years, Temp])

########################## Filter the NDVI data

NDVIDataCoverage = []

for aRecord in YearlyAnalysisNDVI:
	# aRecord is a dataset of a year that was tested against
	# Define the years which the training data is based upon. 
	YearsRaw = aRecord[0]
	YearsFix = []
	for aItem in YearsRaw:
		YearsFix.append(int(aItem))

	NumberOfY = len(YearsFix)

	# Analysis years are items within aRecord[1]
	AY = aRecord[1]

	# Define the key years being tested
	Temp = []
	for TrainYear in AY:
		for Nest in TrainYear:
			Temp.append(int(Nest[1]))
	KeyYears = collections.Counter(Temp).keys()

	# Count how many nest records past the fitering based upon previous filter
	# in limitation loop
	TempData = []
	for KeyYear in KeyYears:
		TempYear = []
		# Now that the key years are individually selected, filter through
		# the nesting data
		for TrainYear in AY:
			# Check if year matches
			YOI = int(TrainYear[0][1])
			# Check if year of interest YOI, is equal to key year
			if YOI==KeyYear:
				# TempYear.append([KeyYear, NumberOfY, len(TrainYear)])
				AppItem = [KeyYear, NumberOfY, len(TrainYear)]
		# Append the data to the overall list
		# DistanceDataCoverage.append(TempYear)
		NDVIDataCoverage.append(AppItem)

# print TotalAnnualNestsDistance
# print YS
GraphData = []
TempYearNum = []
for aItem in NDVIDataCoverage:
	# Define year of focus
	YOF = aItem[0]
	# Check annual nests
	for aRecord in TotalAnnualNestsNDVI:
		if aRecord[0]==YOF:
			Selection = aRecord
	# Calculate percentage pasing
	PercPass = aItem[-1]/float(Selection[1])
	# Write Graphing data
	# FORMAT: Year, number of years in train data, passed recs, perc of total recs
	GraphData.append([YOF, aItem[1], aItem[2], PercPass])
	TempYearNum.append(aItem[1])


###########################  Plot the NDVI

# Define key number of years in analysis
KeyNum = collections.Counter(TempYearNum).keys()
	
# Filter by years in analysis
FilteredGraph = []
Aves = []
STDs = []
Ys=[]

for aItem in KeyNum:
	TempFilt = []
	VarTemp = []
	# Pull records matching aItem, the number of years in analysis
	for aRecord in GraphData:
		if aRecord[1]==aItem:
			TempFilt.append(aRecord)
			VarTemp.append(aRecord[-1])
	FilteredGraph.append(TempFilt)
	# This is the plot data
	Aves.append(np.mean(VarTemp))
	STDs.append(np.std(VarTemp))
	Ys.append(aItem)

# Plot averages with error bar
fig = plt.figure(figsize=(12,12))

plt.scatter(Ys,Aves)
plt.errorbar(Ys, Aves, yerr=STDs, fmt='o')
plt.xlabel('Variables in Combination')
plt.ylabel('Average Ratio of Passed to Total Nests')
plt.title('Average Ratio of Passed Nests (NDVI)\nPer Combination Size')

# Plot the data info
i = 0
for xy in zip(Ys, Aves):
	plt.annotate('  Pass Ratio: %s\n  STD Error: %s' % (round(Aves[i], 2), round(STDs[i], 2)), xy=xy, textcoords='data')
	i = i+1
plt.savefig((f1out+'NDVI_Passed.png'), dpi=None, 
	facecolor='w', edgecolor='b', orientation='portrait', papertype=None, 
	format=None, transparent=False, bbox_inches=None, pad_inches=0.1, 
	frameon=None)
# plt.show()
plt.close()
# sys.exit()

########################################## Test Distance

# For counting passed nests
YearlyAnalysisDistance = []
# Counting annual nests in total
TotalAnnualNestsDistance = []

# Count how many nests are actually in a given year (with distance data)
for aYear in YS:
	TotalNests = 0
	for Nrec in DistanceDataYear:
		if int(Nrec[1])==aYear:
			TotalNests = TotalNests+1
	TotalAnnualNestsDistance.append([aYear, TotalNests])


# Pull trained data which excludes each specific year.
# print NDVIDataYear[0]
# sys.exit()
ProgBarLimit = len(Distance)
print '\nTesting Distance Limitations...\n'
for kk in tqdm.tqdm(range(ProgBarLimit)):
# for aRecord in NDVI:
	aRecord = Distance[kk]
	# temp will be all the data which can be tested for a single year
	Temp = []
	# aRecord is a training set for either a single year or combination.
	# Test which do not belong in that year.
	Years = aRecord[3::]
	for aYear in YS:
		FilteredNests = []
		# Select one of our analysis years.
		# Pull all nests for this year (NDVI only)
		YNDVI = []
		for Nrec in DistanceDataYear:
			if int(Nrec[1])==aYear:
				YNDVI.append(Nrec)

		# Now use this year to identify if the NDVI training record
		# does not include nest data from this year. If so, it can be used
		# to test nest data from aYear as it is not the training data.
		# Convert years to a list of integers to remove the issue with spaces
		YearsFix = []
		for aValue in Years:
			YearsFix.append(int(aValue))

		if aYear not in YearsFix:
			# We can use this.
			# Test how many nests are within the data range
			MaxDistance = float(aRecord[0])
			for aNest in YNDVI:
				if float(aNest[-1])<=MaxDistance:
					FilteredNests.append(aNest)
		if len(FilteredNests)>0:
			Temp.append(FilteredNests)

	# At this point, Temp is the years which can be tested with trained set
	# Define training set (Years)
	YearlyAnalysisDistance.append([Years, Temp])

######################### 
DistanceDataCoverage = []

for aRecord in YearlyAnalysisDistance:
	# aRecord is a dataset of a year that was tested against
	# Define the years which the training data is based upon. 
	YearsRaw = aRecord[0]
	YearsFix = []
	for aItem in YearsRaw:
		YearsFix.append(int(aItem))

	NumberOfY = len(YearsFix)

	# Analysis years are items within aRecord[1]
	AY = aRecord[1]

	# Define the key years being tested
	Temp = []
	for TrainYear in AY:
		for Nest in TrainYear:
			Temp.append(int(Nest[1]))
	KeyYears = collections.Counter(Temp).keys()

	# Count how many nest records past the fitering based upon previous filter
	# in limitation loop
	TempData = []
	for KeyYear in KeyYears:
		TempYear = []
		# Now that the key years are individually selected, filter through
		# the nesting data
		for TrainYear in AY:
			# Check if year matches
			YOI = int(TrainYear[0][1])
			# Check if year of interest YOI, is equal to key year
			if YOI==KeyYear:
				# TempYear.append([KeyYear, NumberOfY, len(TrainYear)])
				AppItem = [KeyYear, NumberOfY, len(TrainYear)]
		# Append the data to the overall list
		# DistanceDataCoverage.append(TempYear)
		DistanceDataCoverage.append(AppItem)

# print TotalAnnualNestsDistance
# print YS
GraphData = []
TempYearNum = []
for aItem in DistanceDataCoverage:
	# Define year of focus
	YOF = aItem[0]
	# Check annual nests
	for aRecord in TotalAnnualNestsDistance:
		if aRecord[0]==YOF:
			Selection = aRecord
	# Calculate percentage pasing
	PercPass = aItem[-1]/float(Selection[1])
	# Write Graphing data
	# FORMAT: Year, number of years in train data, passed recs, perc of total recs
	GraphData.append([YOF, aItem[1], aItem[2], PercPass])
	TempYearNum.append(aItem[1])

###########################  Plot the Distance Data

# Define key number of years in analysis
KeyNum = collections.Counter(TempYearNum).keys()
	
# Filter by years in analysis
FilteredGraph = []
Aves = []
STDs = []
Ys=[]

for aItem in KeyNum:
	TempFilt = []
	VarTemp = []
	# Pull records matching aItem, the number of years in analysis
	for aRecord in GraphData:
		if aRecord[1]==aItem:
			TempFilt.append(aRecord)
			VarTemp.append(aRecord[-1])
	FilteredGraph.append(TempFilt)
	# This is the plot data
	Aves.append(np.mean(VarTemp))
	STDs.append(np.std(VarTemp))
	Ys.append(aItem)

# Plot averages with error bar
fig = plt.figure(figsize=(12,12))

plt.scatter(Ys,Aves)
plt.errorbar(Ys, Aves, yerr=STDs, fmt='o')
plt.xlabel('Variables in Combination')
plt.ylabel('Average Ratio of Passed to Total Nests')
plt.title('Average Ratio of Passed Nests\nPer Combination Size')

# Plot the data info
i = 0
for xy in zip(Ys, Aves):
	plt.annotate('  Pass Ratio: %s\n  STD Error: %s' % (round(Aves[i], 2), round(STDs[i], 2)), xy=xy, textcoords='data')
	i = i+1
plt.savefig((f1out+'Distance_Passed.png'), dpi=None, 
	facecolor='w', edgecolor='b', orientation='portrait', papertype=None, 
	format=None, transparent=False, bbox_inches=None, pad_inches=0.1, 
	frameon=None)
# plt.show()
plt.close()



