import sys
sys.path.append("PythonFiles")
import Modules
from Modules import*
import gzip
import shutil
import tarfile

############################## Unpack the download files #############################

# List Directory to find the correct folder used for the downloaded images
# This will include 'epsa'

File_Folders = list(os.listdir('Raw_Images/'))

for aItem in File_Folders:
	Sp1 = aItem.split('-')
	if Sp1[0] == 'espa':
		ParentDir = 'Raw_Images/'+aItem+'/'
		FilesList = list(os.listdir('Raw_Images/'+aItem))

# Unzip the files in the folder

ProgBarLimit = len(FilesList)
print '\nExtracting Zipped .gz files...\n'
Out = 'Raw_Images/'
for i in tqdm.tqdm(range(ProgBarLimit)):
	fname = ParentDir+FilesList[i]
	if fname.endswith("tar.gz"):
	    tar = tarfile.open(fname, "r:gz")
	    tar.extractall(ParentDir)
	    tar.close()
	elif fname.endswith("tar"):
	    tar = tarfile.open(fname, "r:")
	    tar.extractall(ParentDir)
	    tar.close()

############################## Rename and Move the tif Files #############################

# Write the organized tif file folder.

dir = 'TIF_Code_Images/'
if not os.path.exists(dir):
    os.makedirs(dir)
# Define output paths for this QL/Phase
f1out = 'TIF_Code_Images/'

# List all the Tif files
FilesList = list(os.listdir(ParentDir))

TifFiles = []

BaseCodes = []

print '\nRenaming and Moving .tif Files...\n'
ProgBarLimit = len(FilesList)
for i in tqdm.tqdm(range(ProgBarLimit)):
# for aFile in FilesList:
	aFile = FilesList[i]
	# print aFile
	Sp1 = aFile.split('.')
	if Sp1[-1]!= 'gz':
		if Sp1[-1]!='tar':
			Sp2 = aFile.split('_')
			BaseCode = Sp2[0]+'_'+Sp2[2]+'_'+Sp2[3]

			dir = 'TIF_Code_Images/'+BaseCode+'/'
			if not os.path.exists(dir):
			    os.makedirs(dir)
			# Define output paths for this QL/Phase
			f1out = 'TIF_Code_Images/'+BaseCode+'/'
			if Sp1[-1]=='tif':
				# TifFiles.append(aFile)
				# Split the image code and remove information to match
				# previous program version
				# Rename the file
				if Sp2[-1]=='qa.tif':
					Name = Sp2[0]+'_'+Sp2[2]+'_'+Sp2[3]+'_'+Sp2[7]+'_'+Sp2[8]
					os.rename(ParentDir+aFile, f1out+Name)
				if Sp2[-1]!='qa.tif':
					# print aFile
					# print Sp2
					Name = Sp2[0]+'_'+Sp2[2]+'_'+Sp2[3]+'_'+Sp2[8]
					os.rename(ParentDir+aFile, f1out+Name)


################################# Move the metadata files ########################

# Will do later, no need right now!				

