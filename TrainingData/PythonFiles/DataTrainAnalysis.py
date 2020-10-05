
import Modules
from Modules import*

def TrainDataAnalysis(Option, Years):

	# File loc for either meta data or single year
	if Option =='Train All':
		# train all available years and conduct the analysis
		FileLocND = 'MasterData/TrainData/CombinationsDistance_Meta.txt'
		# FileLocIM = 'MasterData/NestImInfoMaster.csv'
		FileLocNDVI = 'MasterData/TrainData/CombinationsNDVI_Meta.txt'
		f1out = 'MasterData/TrainData/Combinations'

		TrainSpecsOut = 'MasterData/'

	if Option =='Train Select':
		# train select available years and nests conduct the analysis
		FileLocND = 'Analysis_Files/CombinationsDistance_Meta.txt'
		FileLocIM = 'Analysis_Files/NestImInfo.csv'
		FileLocNDVI = 'Analysis_Files/CombinationsNDVI_Meta.txt'
		dir = 'TrainData/Combinations/Analysis'
		if not os.path.exists(dir):
		    os.makedirs(dir)
		f1out = 'TrainData/Combinations/Analysis'


	# Write a folder location for the train data analysis figures.

	# Open datasets
	NDVI = []
	with open(FileLocNDVI, 'r') as Doc:
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
	with open(FileLocND, 'r') as Doc:
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

	# Plot distributions

	NDVIdists = []
	for aItem in NDVI:
		NDVIdists.append(aItem[2])

	Vars = collections.Counter(NDVIdists).keys()
	Freqs = collections.Counter(NDVIdists)
	# print Freqs
	# print Vars
	FreqCount = []
	for aValue in Vars:
		FreqCount.append(Freqs[aValue])



	############################## Plotting Pie

	# Plot the data over a subplot2grid setup
	gridsize = (4,6)
	fig = plt.figure(figsize=(12,12))
	ax1 = plt.subplot2grid(gridsize, (0, 0), colspan= 2, rowspan= 2)
	ax2 = plt.subplot2grid(gridsize, (0, 4), colspan = 2, rowspan=2)
	ax3 = plt.subplot2grid(gridsize, (3, 0), colspan = 2, rowspan = 3)
	ax4 = plt.subplot2grid(gridsize, (3, 4), colspan = 2, rowspan = 3)

	# Plot 1

	wedges1, texts1 = ax1.pie(FreqCount)

	# bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
	# kw = dict(arrowprops=dict(arrowstyle="-"),
	#           bbox=bbox_props, zorder=0, va="center")
	kw = dict(arrowprops=dict(arrowstyle="-", color = 'k'), zorder=0, va="center", size = 7)

	for i, p in enumerate(wedges1):
	    ang = (p.theta2 - p.theta1)/2. + p.theta1
	    y = np.sin(np.deg2rad(ang))
	    x = np.cos(np.deg2rad(ang))
	    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
	    connectionstyle = "angle,angleA=0,angleB={}".format(ang)
	    kw["arrowprops"].update({"connectionstyle": connectionstyle})
	    ax1.annotate(Vars[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y),
	                horizontalalignment=horizontalalignment, **kw)
	# ax1.pie(Sizes18, labels=labels18)
	ax1.set_title('NDVI Distributions', size=12, pad = 15, loc = 'center')


	# Plot 2

	Ddists = []
	for aItem in Distance:
		Ddists.append(aItem[2])

	Vars2 = collections.Counter(Ddists).keys()
	Freqs2 = collections.Counter(Ddists)
	# print Freqs
	# print Vars
	FreqCount2 = []
	for aValue in Vars2:
		FreqCount2.append(Freqs2[aValue])

	wedges2, texts2 = ax2.pie(FreqCount2)

	# bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
	# kw = dict(arrowprops=dict(arrowstyle="-"),
	#           bbox=bbox_props, zorder=0, va="center")
	kw = dict(arrowprops=dict(arrowstyle="-", color = 'k'), zorder=0, va="center", size = 7)

	for i, p in enumerate(wedges2):
	    ang = (p.theta2 - p.theta1)/2. + p.theta1
	    y = np.sin(np.deg2rad(ang))
	    x = np.cos(np.deg2rad(ang))
	    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
	    connectionstyle = "angle,angleA=0,angleB={}".format(ang)
	    kw["arrowprops"].update({"connectionstyle": connectionstyle})
	    ax2.annotate(Vars2[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y),
	                horizontalalignment=horizontalalignment, **kw)

	ax2.set_title('Distance Distributions', size=12, pad = 15, loc = 'center')

	# PLOT 3

	# Add distributions per number of years in analysis
	NDVIYD = []
	NumYears = []


	for aRec in NDVI:
		Vals = len(aRec[3::])
		Dist = aRec[2]
		NDVIYD.append([Dist, Vals])
		NumYears.append(Vals)

	KeyNumbers = collections.Counter(NumYears).keys()

	DistsPerYear = []
	CountOfCombs = []
	for aKeyVar in KeyNumbers:
		Temp = []
		Count=0
		for Vals in NDVIYD:
			# Add dists to list for number of years associated with it
			# then count how many dists are associated with the number of years
			if Vals[-1]==aKeyVar:
				Temp.append(Vals[0])
				Count=Count+1
		KeyDists = collections.Counter(Temp).keys()
		DistsPerYear.append(len(KeyDists))
		CountOfCombs.append(Count)

	ax3.bar(KeyNumbers, DistsPerYear)
	ax3.set_xlabel('Years in Analysis', size=8)
	ax3.set_ylabel('Number of Distributions', size=8)
	ax3.set_title('Number of Distributions\n Per Length of Combination (NDVI)')

	# NOt using this, but I like it. Allows you to loop and write a string in
	# rows.
	# StrInsert = ''
	# i = 0
	# while i<len(CountOfCombs):
	# 	StrInsert = StrInsert+'\n'+str(KeyNumbers[i])+' Year(s); '+ str(CountOfCombs[i])+' Combination(s)'
	# 	i = i+1

	# ax3.text(6,0,StrInsert)

	# PLOT 4


	# Add distributions per number of years in analysis
	DistanceYD = []
	NumYears = []


	for aRec in Distance:
		Vals = len(aRec[3::])
		Dist = aRec[2]
		DistanceYD.append([Dist, Vals])
		NumYears.append(Vals)

	KeyNumbers = collections.Counter(NumYears).keys()

	DistsPerYear = []
	CountOfCombs = []
	for aKeyVar in KeyNumbers:
		Temp = []
		Count=0
		for Vals in NDVIYD:
			# Add dists to list for number of years associated with it
			# then count how many dists are associated with the number of years
			if Vals[-1]==aKeyVar:
				Temp.append(Vals[0])
				Count=Count+1
		KeyDists = collections.Counter(Temp).keys()
		DistsPerYear.append(len(KeyDists))
		CountOfCombs.append(Count)

	ax4.bar(KeyNumbers, DistsPerYear)
	ax4.set_xlabel('Years in Analysis', size=8)
	ax4.set_ylabel('Number of Distributions', size=8)
	ax4.set_title('Number of Distributions\n Per Length of Combination (Distance)')

	plt.savefig((f1out+'Distributions_Trained.png'), dpi=None, 
		facecolor='w', edgecolor='b', orientation='portrait', papertype=None, 
		format=None, transparent=False, bbox_inches=None, pad_inches=0.1, 
		frameon=None)

	# plt.show()
	plt.close()

	# Observe on plot fig, only looking at multiple years together
	X = []
	for aItem in NDVI:
		if len(aItem[3::]) > 1:
			X.append(len(aItem[3::]))
		if len(aItem[3::])==len(Years):
			FullNDVITrainRecord=aItem
	############# MAX NDVI
	# Plot averages with error bar
	fig = plt.figure(figsize=(12,12))

	# Calculate standard deviation for plot
	KeyYears = collections.Counter(X).keys()
	KeyY = collections.Counter(X)
	# print KeyY
	# sys.exit()
	Stds = []
	Aves = []
	for NumbYears in KeyYears:
		Temp = []
		for aRec in NDVI:
			if len(aRec[3::])==NumbYears:
				Temp.append(float(aRec[0]))
		# Calculate STD for that set
		STD = np.std(Temp)
		Ave = np.mean(Temp)
		Stds.append(STD)
		Aves.append(Ave)

	plt.scatter(KeyYears,Aves)
	plt.errorbar(KeyYears, Aves, yerr=Stds, fmt='o')
	plt.xlabel('Variables in Combination')
	plt.ylabel('Average Maximum NDVI')
	plt.title('Average Maximum NDVI Values\nPer Combination Size')
	plt.savefig((f1out+'NDVI_Trained_Max.png'), dpi=None, 
		facecolor='w', edgecolor='b', orientation='portrait', papertype=None, 
		format=None, transparent=False, bbox_inches=None, pad_inches=0.1, 
		frameon=None)
	# plt.show()
	plt.close()

	################# MIN NDVI
	# Plot averages with error bar
	fig = plt.figure(figsize=(12,12))

	# Calculate standard deviation for plot
	KeyYears = collections.Counter(X).keys()
	KeyY = collections.Counter(X)
	# print KeyY
	# sys.exit()
	Stds = []
	Aves = []
	for NumbYears in KeyYears:
		Temp = []
		for aRec in NDVI:
			if len(aRec[3::])==NumbYears:
				Temp.append(float(aRec[1]))
		# Calculate STD for that set
		STD = np.std(Temp)
		Ave = np.mean(Temp)
		Stds.append(STD)
		Aves.append(Ave)

	plt.scatter(KeyYears,Aves)
	plt.errorbar(KeyYears, Aves, yerr=Stds, fmt='o')
	plt.xlabel('Variables in Combination')
	plt.ylabel('Average Minimum NDVI')
	plt.title('Average Minimum NDVI Values\nPer Combination Size')
	plt.savefig((f1out+'NDVI_Trained_Min.png'), dpi=None, 
		facecolor='w', edgecolor='b', orientation='portrait', papertype=None, 
		format=None, transparent=False, bbox_inches=None, pad_inches=0.1, 
		frameon=None)
	# plt.show()
	plt.close()

	################ MAX DISTANCE

	X = []
	for aItem in Distance:
		if len(aItem[3::]) > 1:
			X.append(len(aItem[3::]))
		if len(aItem[3::])==len(Years):
			FullDistanceTrainRecord=aItem

	# Plot averages with error bar
	fig = plt.figure(figsize=(12,12))

	# Calculate standard deviation for plot
	KeyYears = collections.Counter(X).keys()
	KeyY = collections.Counter(X)

	Stds = []
	Aves = []
	for NumbYears in KeyYears:
		Temp = []
		for aRec in Distance:
			if len(aRec[3::])==NumbYears:
				Temp.append(float(aRec[0]))
		# Calculate STD for that set
		STD = np.std(Temp)
		Ave = np.mean(Temp)
		Stds.append(STD)
		Aves.append(Ave)

	plt.scatter(KeyYears,Aves)
	plt.errorbar(KeyYears, Aves, yerr=Stds, fmt='o')
	plt.xlabel('Variables in Combination')
	plt.ylabel('Average Maximum Distance (Meters)')
	plt.title('Average Maximum Distance Values\nPer Combination Size')
	plt.savefig((f1out+'Distance_Trained_Max.png'), dpi=None, 
		facecolor='w', edgecolor='b', orientation='portrait', papertype=None, 
		format=None, transparent=False, bbox_inches=None, pad_inches=0.1, 
		frameon=None)
	# plt.show()
	plt.close()


	###################
	Recs = [FullNDVITrainRecord, FullDistanceTrainRecord]
	# Write out training data to a text doc to predict habitat
	with open(TrainSpecsOut+'Trained_Meta.txt', 'w') as Doc:
		for aValue in Recs:
			Doc.write(str(aValue))
			Doc.write(',')
	Doc.close()



