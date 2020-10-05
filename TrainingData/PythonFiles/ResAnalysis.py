import Modules

from Modules import*


File = 'Analysis_Files/ResEl.csv'

Data = []

with open(File, 'r') as Doc:
	for aItem in Doc:
		# Split the data
		Sp1 = aItem.split(',')
		FloatsToConvert = Sp1[1::]
		Converted = []
		# Convert elevation data to float
		for aValue in FloatsToConvert:
			Converted.append(float(aValue))
		# insert the year into converted
		Converted.insert(0, int(Sp1[0]))
		# Append to data list
		Data.append(Converted)
Doc.close()

# Exclude all years for the study period
# Run for the years of the analysis
AnalysisYears = [2014, 2015, 2016, 2017, 2018, 2019]
FilterData = []

for Record in Data:
	Pass = 0
	for aYear in AnalysisYears:
		if Record[0]==aYear:
			Pass = 'NOOOOO'
	if Pass==0:
		FilterData.append(Record)

# Write def to compute the stds 

def CalculateAveStats(Data):
	# Separate data by month to compute the averages
	Jan = []
	Feb = []
	Mar = []
	Apr = []
	May = []
	Jun = []
	Jul = []
	Aug = []
	Sep = []
	Oct = []
	Nov = []
	Dec = []
	Months = [Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec]

	for aYear in Data:
		Jan.append(aYear[1])
		Feb.append(aYear[2])
		Mar.append(aYear[3])
		Apr.append(aYear[4])
		May.append(aYear[5])
		Jun.append(aYear[6])
		Jul.append(aYear[7])
		Aug.append(aYear[8])
		Sep.append(aYear[9])
		Oct.append(aYear[10])
		Nov.append(aYear[11])
		Dec.append(aYear[12])

	# Compute the monthly stats
	Stats = []
	for aMonth in Months:
		Mean = np.mean(aMonth)
		StandDev = np.std(aMonth)
		Stats.append([Mean, StandDev])
	# Format: Average, STD
	return Stats

# Run for the years of the analysis
AnalysisYears = [2014, 2015, 2016, 2017, 2018, 2019]


Stats = CalculateAveStats(FilterData)

# Yearly variation data
YVD = []

i = 0 
while i<len(AnalysisYears):
	DistanceFromAve = []
	Aves = []
	Stds = []
	for aValue in Stats:
		Aves.append(aValue[0])
		Stds.append(aValue[1])
	# monthly data for year
	Year = AnalysisYears[i]
	for aYear in Data:
		if aYear[0]==Year:
			MD = aYear
	MonthData = MD[1::]
	DistanceFromMean = []
	z = 0
	while z<len(MonthData):
		DM = MonthData[z]-float(Aves[z])
		STDval = DM/float(Stds[z])
		DistanceFromMean.append([DM, STDval, MonthData[z], Aves[z]])
		z = z+1
	YVD.append(DistanceFromMean)
	i = i+1


# Write a table of data for spring and summer
TableData = []
SprSum = []
i=0
while i<len(AnalysisYears):
	Year = AnalysisYears[i]
	# Select only spring and summer data
	SeasonData = YVD[i][2:8]
	SprSum.append(SeasonData)
	i = i+1

# Plot spring and summer data per year
i = 0
for aYear in AnalysisYears:
	Y1 = []
	Y2 = []
	Y3 = []
	Y4 = []
	for aRecord in SprSum[i]:
		# Month data
		Y1.append(aRecord[2])
		# average data
		Y2.append(aRecord[3])
		# Distance from mean
		Y3.append(aRecord[0])
		# Std from mean
		Y4.append(aRecord[1])

	# Months
	X = []
	for n in range(3,9):
		X.append(n)

	plt.plot(X, Y4)
	plt.xlabel('Month')
	plt.ylabel('Distance From Mean (STD)')
	plt.title('%s: Monthly Reservoir Distance From Mean (1968-2013)\n In Standard Deviations (STD)' % aYear)

	# plt.show()
	plt.close()
	# Index variable
	i = i+1


