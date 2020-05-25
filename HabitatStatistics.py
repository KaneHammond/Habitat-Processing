import Modules
from Modules import*

# Open the NDVI dataset 

NDVIdata = []

with open('Analysis_Files/NestNDVI.csv') as File:
	for aRow in File:
		# Clean up the data before appending to item
		Sp1 = aRow.split(',')
		Val = Sp1[2].replace('\n', '')
		try:
			# If converting to float fails, append fail
			Val = float(Val)
			NDVIdata.append(Val)
		except:
			NDVIdata.append('Fail')
File.close()

# Remove the first record, it is the header
NDVIdata = NDVIdata[1::]

plt.hist(NDVIdata, bins = 100)
plt.show()
