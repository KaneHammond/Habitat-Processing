# Habitat-Processing
Shore Bird Habitat Model

Order of scripts:
  1. SegmentWrite.py
  2. SelectScene.py
  3. Image_Request.py (Mod USGS Script)
  4. LandsatDownload.py (USGS Script)
  5. ImageryPrep.py
  6. Nest_Prep.py
  7. ExtractRasterData.py
  8. HabitatStatistics.py
  
# All scripts are run individually. Not set up for full automation at the moment.
  
The first step is to use the segment write file to determine the extent of the analysis.
From here, the appropriate scenes will be selected to match the extent of the analysis.
Starting by identifying the correct path and row of the Landsat 8 sensor. The Image request
will filter through a csv containing all Landsat 8 images. This csv is not written into to
program. Meaning it will not automatically be downloaded. It is a large file and I do not
intend to keep it this way. Once the correct images have been identified for the area and 
date (specified in Image_Request) the request will be sent to USGS. Once completed, you recieve 
an email that the request has been completed. The LandsatDownload.py will allow you to download this.
ImageyPrep.py will unpack and organize the images. The names are also changed to match with a previous
version of the program. Once the images have been organized, Nest_Prep.py should be run. This will
convert the csv nest data (Nest_CSV Folder) into a shapefile. Not provided as a data share agreement
for the data is required. ExtractRasterData.py will determine the NDVI value at nesting points, as
well as the distance to a water source (determined via pixel quality band). HabitatStatistics.py is
to analyze the data to train a model. Currently under development.

