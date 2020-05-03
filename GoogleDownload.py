import Modules
from Modules import *

sys.path.append('Analysis_Files')
Analysis_Files = 'Analysis_Files/'

dir = 'Analysis_Files/Image_Meta'
if not os.path.exists(dir):
    os.makedirs(dir)
Image_Meta_Dir = Analysis_Files+'/'+'Image_Meta'+'/'

# Season date logic. This section will allow for an indepth download of seasonal
# data for iteration. In order to select seasons, start/mid/end image dates must
# be defined. This will allow for a chronological assessment, as well as additional
# error checks for cloud coverage.

# Month range for season start range
# SeasonST1 = '5-1'
# SeasonST2 = '5-15'
# SeasonST3 = '5-15'
# SeasonST4 = '5-31'
# # Mid season range
# SeasonM1 = '6-1'
# SeasonM2 = '6-30'
# # Mid season range s
# SeasonM3 = '7-1'
# SeasonM4 = '7-31'
# # End season range
# SeasonE1 = '8-1'
# SeasonE2 = '8-31'

# SL = [[SeasonST1, SeasonST2], [SeasonM1, SeasonM2], [SeasonM3, SeasonM4], [SeasonE1, SeasonE2]]
SeasonST = '5-1'
SeasonE = '8-31'
SL = [[SeasonST, SeasonE]]

# Year range for the assessment
# Years = ['2000', '2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008',
    # '2009', '2010', '2011','2012', '2014', '2015', '2016', '2017', '2018', '2019']
# Years = ['2001', '2014']
Years = ['2014']

########################################################################

# Download single scene at different dates

# Define the download process for the loop
def filteredImageInCHIRPSToMapId(dateFrom, dateTo, Path, Row, Y):
    Y = int(Y)
    if Y < 2013:
        eeCollection = ee.ImageCollection("LANDSAT/LT05/C01/T1_SR")
        Band1 = 'B3'
        Band2 = 'B4'
    if Y > 2013:
        Band1 = 'B4'
        Band2 = 'B5'
        eeCollection = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR")
    Band3 = 'pixel_qa'

    eeFilterDate = ee.Filter.date(dateFrom, dateTo)
    eeCollection = eeCollection.filter(ee.Filter.eq('WRS_PATH', Path))
    eeCollection = eeCollection.filter(ee.Filter.eq('WRS_ROW', Row))
    # eeCollection = eeCollection.filter(ee.Filter.lt('CLOUD_COVER', 50))
    eeCollection = eeCollection.filter(eeFilterDate)

    Info = eeCollection.getInfo()

    i = 0
    while i < len(Info):
        ImageProperties = Info['features'][i]['properties']

        # Info = eeCollection.getInfo()['features'][0]['properties']
        FirstIm = [str(Info['features'][i]['id'])]
        FirstIm = ee.Image(FirstIm).select([Band1])

        SecondIm = [str(Info['features'][i]['id'])]
        SecondIm = ee.Image(SecondIm).select([Band2])

        ThirdIm = [str(Info['features'][i]['id'])]
        ThirdIm = ee.Image(ThirdIm).select([Band3])

        # Format a file name using landsat image name. This will help
        # filter image by date and P/R when pulling from Google Drive.

        Name = str(Info['features'][i]['id'])
        Split = Name.split('/')
        BaseFileName = Split[-1]

        # Suppress put section to write image properties without requesting 
        # images

        # out = batch.Export.image(FirstIm, description=BaseFileName+'_'+Band1)
        #     # region=ee.Feature(FirstIm.first()).geometry().bounds().getInfo()['coordinates'])
        # ## Process the image
        # process = batch.Task.start(out)
        # process

        # out = batch.Export.image(SecondIm, description=BaseFileName+'_'+Band2)
        #     # region=ee.Feature(FirstIm.first()).geometry().bounds().getInfo()['coordinates'])
        # ## Process the image
        # process = batch.Task.start(out)
        # process

        # out = batch.Export.image(ThirdIm, description=BaseFileName+'_'+Band3)
        #     # region=ee.Feature(FirstIm.first()).geometry().bounds().getInfo()['coordinates'])
        # ## Process the image
        # process = batch.Task.start(out)
        # process

        # Convert dictionary object to list
        list_key_value = [ [k,v] for k, v in ImageProperties.items() ]

        # Write image band properties file
        with open(Image_Meta_Dir+BaseFileName+'_Properties.txt', 'w') as Doc:
            for Object in list_key_value:
                if Object!= list_key_value[-1]:
                    Doc.write(str(Object)+'*')
                if Object== list_key_value[-1]:
                    Doc.write(str(Object))
        i = i+1

#################################################################################
#                   Identify the PR values from other modules
PR_For_Segs = []
All_PR = []
with open(Analysis_Files+'Analysis_Meta.txt', 'r') as Doc:
    for aItem in Doc:
        Temp = aItem
    # CLean up data, split by -
    Temp = Temp.split('-')
    # Select the PR data
    Temp = Temp[1]
    # Split it by , to separate individual PR values per segment
    Temp = Temp.split(',')
    # Each item is now an individual set for each segment
    for aItem in Temp:
        # Drop the PR lable
        aItem = aItem.replace('PR:', '')
        # Split by . to separate segment number from PR
        aItem = aItem.split('.')
        # Split our path row data
        PR_Split = aItem[-1].split('/')
        # Define our segment
        Sg = aItem[0]
        # Write an insertable item
        InsertItem = [int(PR_Split[0]), int(PR_Split[1])]
        ## Format [Seg #, Path, Row]
        All_PR.append(InsertItem)

# Filter through the Segs to ensure no doubles downloaded
for List in All_PR:
    if List not in PR_For_Segs:
        PR_For_Segs.append(List)

##################################################################################
#                   Run the Download
ProgBarLimit0 = len(PR_For_Segs)
ProgBarLimit1 = len(Years)
ProgBarLimit2 = len(SL)

ee.Initialize()

ErrorLog = []
print '\nRequesting Image Bands...\n'
for r in tqdm.tqdm(range(ProgBarLimit0)):
    # Select the correct PR
    Path = PR_For_Segs[r][0]
    Row = PR_For_Segs[r][1]
    for i in tqdm.tqdm(range(ProgBarLimit1)):
        Y = Years[i]
        for t in tqdm.tqdm(range(ProgBarLimit2)):
            Sep = SL[t]
            dateFrom = str(Y+'-'+Sep[0])
            dateTo = str(Y+'-'+Sep[1])
            # filteredImageInCHIRPSToMapId(dateFrom, dateTo, Path, Row, Y)
            try:
                filteredImageInCHIRPSToMapId(dateFrom, dateTo, Path, Row, Y)
            except Exception as e:
                Item = [Sep, dateFrom, dateTo, e]
                ErrorLog.append(Item)
                pass

# Write error log for the Google request
with open(Analysis_Files+'EarthEngine_Image_Request_Error.txt', 'w') as Doc:
    Doc.write(str(ErrorLog))

