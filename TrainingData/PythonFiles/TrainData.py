import Modules
from Modules import*

################################################
# Will train data based upon input nesting data.
# If only segments were selected, then train select
# will run and only train based upon those areas.

Setting = '855'


def TrainTheData(Years, Option):
	warnings.filterwarnings("ignore", category=RuntimeWarning)

	# File loc for either meta data or single year
	if Option =='Train All':
		# train all available years and conduct the analysis
		FileLocND = 'Analysis_Files/NestDistance.csv'
		FileLocIM = 'Analysis_Files/NestImInfo.csv'
		FileLocNDVI = 'Analysis_Files/NestNDVI.csv'
		dir = 'MasterData/TrainData/Combinations'
		if not os.path.exists(dir):
		    os.makedirs(dir)
		MainOut = 'MasterData/TrainData/Combinations'


	if Option =='Train Select':
		# train select available years and nests conduct the analysis
		FileLocND = 'Analysis_Files/NestDistance.csv'
		FileLocIM = 'Analysis_Files/NestImInfo.csv'
		FileLocNDVI = 'Analysis_Files/NestNDVI.csv'
		dir = 'TrainData/Combinations/'
		if not os.path.exists(dir):
		    os.makedirs(dir)
		MainOut = 'TrainData/Combinations/'

	# Define statistics to be used for specific years. The meta data will be saved
	# then tested separately.
	# Years = ['2014', '2015', '2016', '2017', '2018', '2019']

	# Pull the nest data to extract dates and nest info. This will be used to filter out
	# data later to select only NDVI data which will be relevant for testing each year
	# individually.

	NestImageInfo = []

	with open(FileLocIM) as File:
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
	# Open the distance data to conduct similar analysis to determine most probable 
	# distance from water.
	DistanceData = []

	with open(FileLocND) as File:
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

	with open(FileLocNDVI) as File:
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


	#################################### Defs #########################################

	BinMethods = ['auto', 'doane', 'scott', 'stone', 'rice', 'sturges', 'sqrt']

	def TestBins(Array):
		Freq = []
		for aBin in BinMethods:
			# https://numpy.org/devdocs/reference/generated/numpy.histogram_bin_edges.html
			var = np.histogram_bin_edges(Array, bins=aBin)
			# sys.exit()
			y, x = np.histogram(Array, bins = var)
			Freq.append(y)
		return Freq

	# The objective is to have the least amount of empty segments, with more than
	# 5 records on average per bin.
	def FilterBinOptions(Freqs):
		FreqLists = []
		# Index variable
		i = 0
		for aItem in Freqs:
			aItem = aItem.tolist()
			# Count the occurance of less than 5 in a bin
			Count = 0
			for aValue in aItem:
				if aValue>5:
					Count = Count+1
			# Calculate percentage of columns less than 5
			if Count>0:
				Percent = float(Count)/len(aItem)
			if Count == 0:
				# No low items found
				Percent = 0
			# Append to a final list
			FreqLists.append([Percent, Count, i])
			i = i+1
		# Sort by first variable and select the first in list, this is the best bin method
		# of the set
		SortedFreq = sorted(FreqLists, key=lambda x: x[0])
		# Define the index variable from the list of the first item, this will select the most
		# appropriate bin method.
		BinMethodOpt = BinMethods[SortedFreq[0][-1]]
		return BinMethodOpt

	def DetermineDistribution(DataArray, Size, BinMethodOpt):

		# Fails to import: geninvgauss
		# Fails to work: levy_stable has some issues, says it is in an experimental stage.
		# Waited a few minutes and it was not able to process. May just take several minutes.
		# Regardless, it was removed. 5/30/20

		dist_names = ['alpha', 'anglit', 'arcsine', 'argus', 'beta', 'betaprime',
			'bradford', 'burr', 'burr12', 'cauchy', 'chi', 'chi2', 'cosine', 'crystalball',
			'dgamma', 'dweibull', 'erlang', 'expon', 'exponnorm', 'exponweib', 'exponpow',
			'f', 'fatiguelife', 'fisk', 'foldcauchy', 'foldnorm', 'frechet_r', 'frechet_l',
			'genlogistic', 'gennorm', 'genpareto', 'genexpon', 'genextreme', 'gausshyper',
			'gamma', 'gengamma', 'genhalflogistic', 'gilbrat', 'gompertz',
			'gumbel_r', 'gumbel_l', 'halfcauchy', 'halflogistic', 'halfnorm', 'halfgennorm',
			'hypsecant', 'invgamma', 'invgauss', 'invweibull', 'johnsonsb', 'johnsonsu',
			'kappa4', 'kappa3', 'ksone', 'kstwobign', 'laplace', 'levy', 'levy_l', 
			'logistic', 'loggamma', 'loglaplace', 'lognorm', 'lomax',
			'maxwell', 'mielke', 'moyal', 'nakagami', 'ncx2', 'ncf', 'norm', 'norminvgauss',
			'pareto', 'pearson3', 'powerlaw', 'powerlognorm', 'powernorm', 'rdist', 'rayleigh',
			'rice', 'recipinvgauss', 'semicircular', 'skewnorm', 't', 'trapz', 'triang', 
			'truncexpon', 'truncnorm', 'tukeylambda', 'uniform', 'vonmises', 'vonmises_line',
			'wald', 'weibull_min', 'weibull_max', 'wrapcauchy']

		sse_statistics = []
		Params = []

		# 11 equi-distant bins of observed Data 
		# https://numpy.org/devdocs/reference/generated/numpy.linspace.html
		# ******************** Not my code 1
		# KANE NOTE: I do not like this. I also do not understand if this is best.
		# The distribution of this data is being based off a histogram, which has a
		# set number of bins. These may not adequately represent the data. Therefore,
		# the distribution this loop will suggest may not be optimal. A def in Test2.py
		# uses scipy tools to test several bin calculation methods to parse a histogram.
		# Just haven't been able to merge it with this script.

		# percentile_bins = np.linspace(0,100,11)
		# # # https://numpy.org/doc/1.18/reference/generated/numpy.percentile.html
		# percentile_cutoffs = np.percentile(DataArray, percentile_bins)
		# observed_frequency, bins = (np.histogram(DataArray, bins=percentile_cutoffs))
		# cum_observed_frequency = np.cumsum(observed_frequency)
		Parse_Bins = np.histogram_bin_edges(DataArray, bins=BinMethodOpt)
		y, x = np.histogram(DataArray, bins=Parse_Bins, density = True)
		# Clean up x
		x = (x + np.roll(x, -1))[:-1] / 2.0

		# ******************** Not my code 1

		# Loop through candidate distributions
		ProgBarLimit = len(dist_names)
		# print '\nDetermine Distribution...\n'
		for kk in tqdm.tqdm(range(ProgBarLimit)):
			# ******************** Not my code 2
			# for distribution in dist_names:
			distribution = dist_names[kk]
			# print distribution
			# Set up distribution and get fitted distribution parameters
			dist = getattr(scipy.stats, distribution)
			# Sources explaining .fit()
			# https://scikit-learn.org/stable/getting_started.html
			# https://stackoverflow.com/questions/45704226/what-does-fit-method-in-scikit-learn-do
			# https://stackoverflow.com/questions/34727919/fit-method-in-python-sklearn
			# https://scikit-learn.org/stable/tutorial/basic/tutorial.html
			param = dist.fit(DataArray)
			# print("{}\n{}\n".format(dist, param))
			Params.append(param)
			
			pdf = dist.pdf(x, *param)
			# Calculate fitted PDF error 
			sse = np.sum(np.power(y - pdf, 2.0))
			# Save the sse, lowest value is best fit
			sse_statistics.append(sse)




		#Sort by minimum ch-square statistics
		results = pd.DataFrame()
		results['Distribution'] = dist_names
		results['SSE'] = sse_statistics
		results['Mods'] = Params
		results.sort_values(['SSE'], inplace=True)
		return results
		# ******************** Not my code 2

	################################## Permutations and Combinations ######################

	# Compute the permutations. They are given as a list tuple items
	PermsTup = list(itertools.permutations(Years))
	# Convert these into lists
	PermsList = []

	for aTup in PermsTup:
		PermsList.append(list(aTup))

	# Filter through permutations to identifiy identical combinations. These only need
	# to be computed once!
	Ranges = []
	i = 0 
	# Define the number of variables in the list. Subtract 1, otherwise the final
	# combination and permutation loop will include all years available. Not required,
	# there is only one combination for this to begin with.
	# Max = len(PermsList[0])-1
	Max = len(PermsList[0])

	while i<Max:
		if i==0:
			Ranges.append([0])
		if i>0:
			Ranges.append([0, i])
		i = i+1

	# Pull all possible permutations for each range and filter out those which have the
	# the same combination. 

	# All permutation clips
	Filtered = []

	for aRange in Ranges:
		# Define a temp list to add permutation clips to
		TempPerms = []
		if len(aRange)==1:
			#This will be our first one.
			# Syntax to make a single variable work is different. There is no range.
			for perm in PermsList:
				# We are only looking at first variable here.
				TempPerms.append(int(perm[0]))
			# Now that all the first values have been identified, we can simply use
			# a counter import to select the key variables. In this case, they will be 
			# equal to the number of years available.
			KeyVars = collections.Counter(TempPerms).keys()
			# Append these key years to the Filtered list. I jsut realized this part
			# could have been much simpler. But I already wrote it.
			# Add to list
			Filtered.append(KeyVars)

		if len(aRange)==2:
			# We are now testing for identical pairs in a permutation.
			for perm in PermsList:
				# Must add 1 to the last variable to make this work. 
				SelectYears = perm[aRange[0]:aRange[-1]+1]
				# Append this section of the permutation to the temporary dataset
				TempPerms.append(SelectYears)

			# Filter out identical permutations 
			TempFilteri = []
			for aPermClip in TempPerms:
				if aPermClip not in TempFilteri:
					TempFilteri.append(aPermClip)

			# Remove identical combinations. This can be done in the same way as above
			# by converting the data to integer and ordering by value.
			TempFilterii = []
			for aItem in TempFilteri:
				# Convert the variables to integer.
				NewList = []
				for aString in aItem:
					# Convert to integer and append to list
					NewList.append(int(aString))
				# Sort this new list in place
				NewList.sort()
				# Append to list 
				TempFilterii.append(NewList)

			# Filter out identical items again
			TempFilteriii = []
			for aPermClip in TempFilterii:
				if aPermClip not in TempFilteriii:
					TempFilteriii.append(aPermClip)
			# append to main list
			Filtered.append(TempFilteriii)
			# Clear all previous lists
			TempFilteri = None
			TempFilterii = None
			TempFilteriii = None

	############################ Nests based on combinations ##########################

	# Extract all of the nest data combinations using our combination year data

	# Select nests based upon year
	NestsSelect = []
	# Start by selecting a combination of years from the filtered dataset
	for CombSet in Filtered:
		Tempi = []
		# Pull years, combinations from set
		for Values in CombSet:
			# List for the values in this specific Combination set
			Tempii = []
			# SINGLE YEARS HERE
			if type(Values)!=list:
				# It is a list of only years, first record in the set as they
				# are individual years.
				for aRec in NestImageInfo:
					if int(aRec[-1])==Values:
						Tempii.append(aRec[0])
				# This will append all nests with that specific year to the list
				Tempi.append(Tempii)
				Tempii = None
			# LISTS OF YEARS HERE
			if type(Values)==list:
				Tempii = []
				# Append all nests for each year into Tempii
				for aYear in Values:
					for aRec in NestImageInfo:
						if int(aRec[-1]) == aYear:
							Tempii.append(aRec[0])
				# Append all nests for the year selection to Tempi
				Tempi.append(Tempii)
				Tempii = None
		# Append all nest data to the selected nests list. This list will contain nest
		# ids for each combination of years that the program will iterate through.
		NestsSelect.append(Tempi)


	OANDVI = []
	OADistance = []

	# print NestsSelect
	# sys.exit()

	# Select specific set of records first based upon combination length
	# This will take some time...
	ProgBarLimit = len(NestsSelect)
	# ProgBarLimit = 2
	print '\nParsing Nest Data...\n'
	for zz in tqdm.tqdm(range(ProgBarLimit)):
	# for aRecSet in NestsSelect:
		aRecSet = NestsSelect[zz]
		# Define overall list for the data to be appended to.
		# Allowing it to match the format of the others
		TempNDVI = []
		TempDistance = []
		# Select list of nests
		for NestList in aRecSet:
			NDVI = []
			Distance = []
			# Now filter by individual nest
			for aRec in NestList:
				# print aRec
				# sys.exit()
				# filter and find the nests in other datasets
				for Vals in NDVIData:
					if Vals[0]==aRec:
						NDVI.append(Vals[1])
				for Vals2 in DistanceData:
					if Vals2[0]==aRec:
						Distance.append(Vals2[1])
			TempNDVI.append(NDVI)
			TempDistance.append(Distance)
		# Now append this data to a list
		OANDVI.append(TempNDVI)
		OADistance.append(TempDistance)

	##################################### TRAIN DATA ##################################

	# This section will define the parameters for each individual year, as well as
	# necesarry combinations.


	# YearlyStats = []

	# Loop through training data for each combination
	ProgBarLimit = len(OANDVI)
	# ProgBarLimit = 1
	OALimits = []
	print '\nDetermine NDVI Limitations...\n'
	for kk in tqdm.tqdm(range(ProgBarLimit)):
		# Define the NDVI values for specified nest sets
		DataSet = OANDVI[kk]
		# Define year combination list
		YCL = Filtered[kk]
		# Count index to determine what combination is being used
		# Important to know for output
		klk = 0
		# This will record all NDVI limitations for year combinations
		TempSet = []
		for NDVIdata in DataSet:
			Years = YCL[klk]
			# *****************NDVI Analysis
			# STEP 1
			# Define requirements for DetermineDistribution function
			DataArray = np.array(NDVIdata)
			Size = int(len(NDVIdata))
			# DataArray = np.array(Distdata)
			# Size = int(len(Distdata))

			Freqs = TestBins(np.array(NDVIdata))

			BinMethodOpt = FilterBinOptions(Freqs)

			# Initiate the function
			Output = DetermineDistribution(DataArray, Size, BinMethodOpt)

			# Write results to list
			ListofOut = Output.values.tolist()

			# Select the first record, it is the best fit based upon Chi-square test
			DistributionData = ListofOut[0]
			# STEP 2

			# CallFunc will call the specific graphing function needed for the 
			# optimal distribution identified by the DetermineDistribution function.
			# The first value in this list is the string name of the specific distribution.
			CallFunc = getattr(scipy.stats, DistributionData[0])

			Mod = str(DistributionData[-1])
			# Remove the parentheses
			Mod = Mod.replace('(', '')
			Mod = Mod.replace(')', '')
			# Split by comma
			Sp1 = Mod.split(',')

			Mod2 = []
			for aItem in Sp1:
				Mod2.append(float(aItem))

			# Define the distribution parameters
			p = tuple(Mod2)

			# This allows for the calculation of CDF using the identified distribution from above.
			# Sinc it is not a normal distribution, this is important. It is also the main objective.
			obj = CallFunc.cdf(DataArray, *p)
			# obj = CallFunc.cdf(0.1, *p)
			listobj = obj.tolist()

			# Define lim of distribution
			if Setting == '91':
				TopPercentile = CallFunc.ppf(0.9, *p)
				BottomPercentile = CallFunc.ppf(0.1, *p)
				PU = '90%'
				PL = '10%'
			if Setting == '855':
				TopPercentile = CallFunc.ppf(0.85, *p)
				BottomPercentile = CallFunc.ppf(0.05, *p)
				PU = '85%'
				PL = '5%'
			# Add this to yearly stats:
			TempSet.append([TopPercentile, BottomPercentile, DistributionData[0], Years])
			# print TopPercentile
			# print BottomPercentile

			# sys.exit()

			# Order this data together and show a graph of CDF with the optimal distribution function
			CombineData = []

			i = 0
			for aItem in listobj:
				# The order of listobj is the same as the NDVIdata, it is not ordered numerically
				# I want to combine these sets before ordering so I can graph them.
				CombineData.append([aItem, NDVIdata[i]])
				i = i+1

			def Sort(sub_li): 
				# Sort by second variable
				return(sorted(sub_li, key = lambda x: x[1])) 

			SortedDat = Sort(CombineData)

			# Define x and y
			x = []
			y = []

			for aItem in SortedDat:
				x.append(aItem[-1])
				y.append(aItem[0])  

			xa = np.array(x)
			ya = np.array(y)

			fig, ax = plt.subplots(1, 1)
			ax.grid(True)
			# Add in out texts defining out limits.
			ax.text(-0.01, 0.93, '%s' % ('90%'), verticalalignment='top', 
				horizontalalignment='right', transform=ax.transAxes, fontsize=9,
				fontweight='bold', color='r')
			ax.text(-0.01, 0.13, '%s' % ('10%'), verticalalignment='top', 
				horizontalalignment='right', transform=ax.transAxes, fontsize=9,
				fontweight='bold', color='r')
			# Add a text defining the limitations of the NDVI values
			ax.text(1, 0.05, 'NDVI Range: Max = %s Min = %s' % (round(TopPercentile,4), round(BottomPercentile,4)), verticalalignment='bottom', 
				horizontalalignment='right', transform=ax.transAxes, fontsize=9,
				fontweight='bold', color='black')

			ax.text(.5, 1.1, 'NDVI CDF: %s' % (DistributionData[0]), verticalalignment='center', 
				horizontalalignment='center', transform=ax.transAxes, fontsize=12,
				fontweight=None, color='black')
			ax.plot(xa, ya, label='CDF')

			ax.legend(loc='right')
			ax.set_title('Cumulative Distrbution of NDVI: Year(s) %s' % Years)
			ax.set_xlabel('NDVI')
			ax.set_ylabel('CDF Value')
			ax.set_yticks([0.1, 0.9], minor=True)
			# ax.set_yticks([0.90], minor=True)
			# ax.set_yticks([0.10], minor=True)
			# ax.yaxis.grid(True, which='major', color='g', linewidth=0.5)
			ax.yaxis.grid(True, which='minor', color='r', linewidth=1.25)
			plt.savefig((MainOut+str(Years)+'_NDVI_Trained.png'), dpi=None, 
				facecolor='w', edgecolor='b', orientation='portrait', papertype=None, 
				format=None, transparent=False, bbox_inches=None, pad_inches=0.1, 
				frameon=None)
			# plt.show()
			plt.close()
			klk = klk+1
		OALimits.append(TempSet)

	# Write the percentile data for each year to a txt doc
	with open(MainOut+'NDVI_Meta.txt', 'w') as Doc:
		for aValue in OALimits:
			Doc.write(str(aValue))
			Doc.write(',')
	Doc.close()


	##################################Distance Data

	# Copy and paste of NDVI loop!

	# Loop through training data for each combination
	ProgBarLimit = len(OADistance)
	# ProgBarLimit = 1
	OADLimits = []
	print '\nDetermine Distance Limitations...\n'
	for kk in tqdm.tqdm(range(ProgBarLimit)):
		# Define the NDVI values for specified nest sets
		DataSet = OADistance[kk]
		# Define year combination list
		YCL = Filtered[kk]
		# Count index to determine what combination is being used
		# Important to know for output
		klk = 0
		# This will record all NDVI limitations for year combinations
		TempSet = []
		for NDVIdata in DataSet:
			Years = YCL[klk]
			# *****************NDVI Analysis
			# STEP 1
			# Define requirements for DetermineDistribution function
			DataArray = np.array(NDVIdata)
			Size = int(len(NDVIdata))
			# DataArray = np.array(Distdata)
			# Size = int(len(Distdata))

			Freqs = TestBins(np.array(NDVIdata))

			BinMethodOpt = FilterBinOptions(Freqs)

			# Initiate the function
			Output = DetermineDistribution(DataArray, Size, BinMethodOpt)

			# Write results to list
			ListofOut = Output.values.tolist()

			# Select the first record, it is the best fit based upon Chi-square test
			DistributionData = ListofOut[0]
			# STEP 2

			# CallFunc will call the specific graphing function needed for the 
			# optimal distribution identified by the DetermineDistribution function.
			# The first value in this list is the string name of the specific distribution.
			CallFunc = getattr(scipy.stats, DistributionData[0])

			Mod = str(DistributionData[-1])
			# Remove the parentheses
			Mod = Mod.replace('(', '')
			Mod = Mod.replace(')', '')
			# Split by comma
			Sp1 = Mod.split(',')

			Mod2 = []
			for aItem in Sp1:
				Mod2.append(float(aItem))

			# Define the distribution parameters
			p = tuple(Mod2)

			# This allows for the calculation of CDF using the identified distribution from above.
			# Sinc it is not a normal distribution, this is important. It is also the main objective.
			obj = CallFunc.cdf(DataArray, *p)
			# obj = CallFunc.cdf(0.1, *p)
			listobj = obj.tolist()

			# Define lim of distribution
			if Setting == '91':
				TopPercentile = CallFunc.ppf(0.9, *p)
				BottomPercentile = CallFunc.ppf(0.1, *p)
				PU = '90%'
				PL = '10%'
			if Setting == '855':
				TopPercentile = CallFunc.ppf(0.85, *p)
				BottomPercentile = CallFunc.ppf(0.05, *p)
				PU = '85%'
				PL = '5%'
			# Add this to yearly stats:
			TempSet.append([TopPercentile, BottomPercentile, DistributionData[0], Years])
			# print TopPercentile
			# print BottomPercentile

			# sys.exit()

			# Order this data together and show a graph of CDF with the optimal distribution function
			CombineData = []

			i = 0
			for aItem in listobj:
				# The order of listobj is the same as the NDVIdata, it is not ordered numerically
				# I want to combine these sets before ordering so I can graph them.
				CombineData.append([aItem, NDVIdata[i]])
				i = i+1

			def Sort(sub_li): 
				# Sort by second variable
				return(sorted(sub_li, key = lambda x: x[1])) 

			SortedDat = Sort(CombineData)

			# Define x and y
			x = []
			y = []

			for aItem in SortedDat:
				x.append(aItem[-1])
				y.append(aItem[0])  

			xa = np.array(x)
			ya = np.array(y)

			fig, ax = plt.subplots(1, 1)
			ax.grid(True)
			# Add in out texts defining out limits.
			ax.text(-0.01, 0.93, '%s' % (PU), verticalalignment='top', 
				horizontalalignment='right', transform=ax.transAxes, fontsize=9,
				fontweight='bold', color='r')
			ax.text(-0.01, 0.13, '%s' % (PL), verticalalignment='top', 
				horizontalalignment='right', transform=ax.transAxes, fontsize=9,
				fontweight='bold', color='r')
			# Add a text defining the limitations of the NDVI values
			ax.text(1, 0.05, 'Distance Range: Max = %s Min = %s' % (round(TopPercentile,4), round(BottomPercentile,4)), verticalalignment='bottom', 
				horizontalalignment='right', transform=ax.transAxes, fontsize=9,
				fontweight='bold', color='black')

			ax.text(.5, 1.1, 'Distance CDF: %s' % (DistributionData[0]), verticalalignment='center', 
				horizontalalignment='center', transform=ax.transAxes, fontsize=12,
				fontweight=None, color='black')
			ax.plot(xa, ya, label='CDF')

			ax.legend(loc='right')
			ax.set_title('Cumulative Distrbution of Distance: Year(s) %s' % Years)
			ax.set_xlabel('Distance')
			ax.set_ylabel('CDF Value')
			ax.set_yticks([0.1, 0.9], minor=True)
			# ax.set_yticks([0.90], minor=True)
			# ax.set_yticks([0.10], minor=True)
			# ax.yaxis.grid(True, which='major', color='g', linewidth=0.5)
			ax.yaxis.grid(True, which='minor', color='r', linewidth=1.25)
			plt.savefig((MainOut+str(Years)+'_Distance_Trained.png'), dpi=None, 
				facecolor='w', edgecolor='b', orientation='portrait', papertype=None, 
				format=None, transparent=False, bbox_inches=None, pad_inches=0.1, 
				frameon=None)
			# plt.show()
			plt.close()
			klk = klk+1
		OADLimits.append(TempSet)

	# Write the percentile data for each year to a txt doc
	with open(MainOut+'Distance_Meta.txt', 'w') as Doc:
		for aValue in OADLimits:
			Doc.write(str(aValue))
			Doc.write(',')
	Doc.close()
	# Calculate train data on the fly. Each year will be tested with the
	# data before it.

