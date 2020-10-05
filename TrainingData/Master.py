
# Combine and run files all at once
# Define location of python files
import sys
sys.path.append("PythonFiles")

import Modules
from Modules import*

############################ Check for master nest data

# # Download up to date imagery

# dir = 'MasterData'
# if not os.path.exists(dir):
#     os.makedirs(dir)

# MF = []
# try:
# 	MF = os.listdir('MasterData/')
# except:
# 	print '\nMaster Data Files Not Found...'

# if len(MF)==0:
# 	print "No master data is available."

# # Check to which extent the data has been processed.
# if len(MF)>0:
# 	Years = []
# 	with open('MasterData/NestImInfoMaster.csv', 'r') as Doc:
# 		i = 0
# 		for aItem in Doc:
# 			if i>0:
# 				# After skipping first record (header), extract the year
# 				# from the image file link name.
# 				Sp1 = aItem.split('_')
# 				ImageYear = Sp1[-1][0:4]
# 				Years.append(int(ImageYear))
# 			i = i+1
# 	Doc.close()
# 	# Determine key years in training data
# 	KY = collections.Counter(Years).keys()

# 	# Compare years in set to years from 2014
	
# 	CY = int(datetime.now().year)

# 	AllYears = []
# 	for x in range(2014, CY+1):
# 		AllYears.append(x)
	
# 	# Identify the difference in the list item
# 	Diff = list(set(AllYears)-set(KY))
# 	# If there is a year found here, then a data download and extract must 
# 	# be done. #WORK WORK WORK
import SegmentWrite
SegmentWrite
# Select scene
import SelectScene
from SelectScene import*
SelectScenes('Train')
# import Image_Request
# Image_Request
# import Download
# Download
# sys.exit()
# FIX Imagery prep

# Run the nest thing
# import Nest_Prep
# Nest_Prep


# sys.exit()
import ExtractRasterData
ExtractRasterData
# sys.exit()
# Start writting training logic

# Initiate training with given years

Years = ['2014', '2015', '2016', '2017', '2018', '2019', '2020']
Option = 'Train All'
# Option = 'Train Select'
import TrainData
from TrainData import*
# Year data that is given allows for development of combinations
TrainTheData(Years, Option)

# Run the train analysis
import DataTrainAnalysis
from DataTrainAnalysis import*
TrainDataAnalysis(Option, Years)




