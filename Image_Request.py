
import Modules
from Modules import*

import requests
import json
from getpass import getpass

# https://www.usgs.gov/land-resources/nli/landsat/bulk-metadata-service
# Link above is for the Landsat 8 meta data CSV. Not written into program
# for automatic download. This module is also a stand alone from the others.
# It was not written by me. Regardless, it is what I use after my requests 
# have been approved via email.


############################################## Define needed P/R values ######

# *&*&*&*&*&*&*&*& IF need to run request after this section, mute here.*&*&*&*&*&*&*&*&

# Open the file
RawData = []
with open('Analysis_Files/Analysis_Meta.txt', 'r') as File:
	for aRow in File:
		RawData.append(aRow)

File.close()

# print RawData

# Extract all p/r information

PR = []

for aList in RawData:
	Sp1 = aList.split('-')
	Sp2 = Sp1[1].split(',')

	for aItem in Sp2:
		aItem = aItem.split('.')
		# Sp4 = aItem[1].split('/')
		PR.append(aItem[1])


SelectPR = list(set(PR))

# print SelectPR
# sys.exit

# SelectPR.append('097/011')

Years = ['2014', '2015', '2016', '2017', '2018', '2019', '2020']
# Years = ['2015']

# ################################# Filter LS8 image codes #####################
# # FORMAT:
# # ['PANCHROMATIC_LINES', 'NADIR_OFFNADIR', 'sunAzimuth', 'REFLECTIVE_SAMPLES', 
# # 'upperLeftCornerLongitude', 'cloudCover', 'MAP_PROJECTION_L1', 'cartURL', 
# # 'sunElevation', 'path', 'BPF_NAME_TIRS', 'THERMAL_LINES', 
# # 'GROUND_CONTROL_POINTS_MODEL', 'row', 'imageQuality1', 'REFLECTIVE_LINES', 
# # 'ELLIPSOID', 'GEOMETRIC_RMSE_MODEL', 'browseURL', 'browseAvailable', 
# # 'dayOrNight', 'CPF_NAME', 'DATA_TYPE_L1', 'THERMAL_SAMPLES', 
# # 'upperRightCornerLatitude', 'lowerLeftCornerLatitude', 'sceneStartTime', 
# # 'dateUpdated', 'sensor', 'PANCHROMATIC_SAMPLES', 'GROUND_CONTROL_POINTS_VERSION', 
# # 'LANDSAT_PRODUCT_ID', 'acquisitionDate', 'upperRightCornerLongitude', 
# # 'PROCESSING_SOFTWARE_VERSION', 'GRID_CELL_SIZE_REFLECTIVE', 
# # 'lowerRightCornerLongitude', 'lowerRightCornerLatitude', 'sceneCenterLongitude',
# # 'COLLECTION_CATEGORY', 'GRID_CELL_SIZE_PANCHROMATIC', 'BPF_NAME_OLI', 
# # 'sceneCenterLatitude', 'CLOUD_COVER_LAND', 'lowerLeftCornerLongitude', 
# # 'GEOMETRIC_RMSE_MODEL_X', 'GEOMETRIC_RMSE_MODEL_Y', 'sceneStopTime', 
# # 'upperLeftCornerLatitude', 'UTM_ZONE', 'DATE_L1_GENERATED', 
# # 'GRID_CELL_SIZE_THERMAL', 'DATUM', 'COLLECTION_NUMBER', 'sceneID', 
# # 'RLUT_FILE_NAME', 'TIRS_SSM_MODEL', 'ROLL_ANGLE', 'receivingStation']

RequiredFiles = []

i = 0

for df in pd.read_csv('LANDSAT_8_C1.csv', chunksize=100000):
	ps = pd.DataFrame(df)
	AllCodes = ps['LANDSAT_PRODUCT_ID'].tolist()
	# l1 = ps['browseURL'].tolist()
	# l2 = ps['browseAvailable'].tolist()
	# l3 = ps['cartURL'].tolist()
	# l4 = ps['path'].tolist()
	# print l1[0]
	# print l2[0]
	# print l3[0]
	# print l4[0]
	# sys.exit()
	for LCODE in AllCodes:
		# Filter through the data to determine if it is needed.
		# Path Row information
		# LCODE = ps['LANDSAT_PRODUCT_ID'].tolist()
		Sp1 = LCODE.split('_')
		PRextract = str(Sp1[2])
		for aCode in SelectPR:
			if PRextract==aCode:
				# Filter by year in list
				Yearextract = Sp1[3][0:4]
				# print Yearextract
				# sys.exit()
				for aYear in Years:
					if aYear==Yearextract:
						# Filter by month and day
						Me = int(Sp1[3][4:6])
						De = int(Sp1[3][6::])
						if Me <= 8:
							if Me >= 4:
								# Month is acceptable
								if Me!=4:
									if Me!=8:
										RequiredFiles.append(LCODE)
								# If border month, filter by day as well
								if Me == 4:
									if De>15:
										RequiredFiles.append(LCODE)
								if Me==8:
									if De<15:
										RequiredFiles.append(LCODE)
	i = i+1
  # This takes some time and the actual size of file is not known.
  # This will print after every chunk (100000 records). This is global data
  # so the csv is large. Not currently written for automatic download. This will
  # require a little thought. We do not want it do download every time it runs. Only
  # when new data should be available. Also, in order to get new scenes, it will need
  # to download the global image data for filtering. Could be easily be improved upon.
  # Will require more knowledge on the EROS USGS server for advanced filtering options.
	print '\n'
	print i
	print len(AllCodes)
	print AllCodes[0]
	print '\n'

df = pd.DataFrame(RequiredFiles)
df.to_csv('BaseImages.csv')

# *&*&*&*&*&*&*&*& Mute to run request again. This section filters 
# through all sat data but writes to file to reduce time of a second attempt. *&*&*&*&*&*&*&*&

# sys.exit()
########################## Open the Base image csv

RawData = []

with open('BaseImages.csv', 'r') as File:
	for aRow in File:
		RawData.append(aRow)
		# sys.exit()
File.close()
RawData = RawData[1::]

AllImages = []
for aItem in RawData:
	Sp1 = aItem.split(',')
	Sp2 = Sp1[1].split('_')
	Q = Sp2[-1].replace('\n', '')
	if Q=='T1':
		Code = Sp1[1].replace('\n', '')
		AllImages.append(Code)

# Manually remove restricted data if applicable. This is returned upon checking
# the request prior to order. If ordered without checking, the order will fail
# and provide no reason as to why. Other than stating it was a bad
# request.
Remove = ["LC08_L1TP_030030_20160609_20170223_01_T1",
"LC08_L1TP_032028_20160607_20170223_01_T1",
"LC08_L1TP_032029_20160607_20170223_01_T1",
"LC08_L1TP_032027_20160607_20170223_01_T1",
"LC08_L1TP_034027_20160605_20170223_01_T1",
"LC08_L1TP_029030_20160602_20170223_01_T1",
"LC08_L1TP_031029_20160531_20170223_01_T1",
"LC08_L1TP_031030_20160531_20170223_01_T1"]


# We have 329 collections 2014-2020 FullR
Limit = len(AllImages)
# Limit = 50

FilteredImages = []
# i = 0
# while i<Limit:
#   FilteredImages.append(AllImages[i])
#   i = i+1

for aItem in AllImages:
  c = 0
  for aIm in Remove:
    if aItem==aIm:
      c = 1
  if c == 0:
    FilteredImages.append(aItem)

AllImages = None

# sys.exit()

###############################################

# Consider checking the version here, online is listed as v1: https://espa.cr.usgs.gov/static/docs/examples/api_demo.py
host = 'https://espa.cr.usgs.gov/api/v1/'
# host = 'https://espa.cr.usgs.gov/api/v0/'
username = 'EnterUsername'
password = 'EnterPassword'

# First and foremost, define a simple function for interacting with the API. 

def espa_api(endpoint, verb='get', body=None, uauth=None):
    """ Suggested simple way to interact with the ESPA JSON REST API """
    auth_tup = uauth if uauth else (username, password)
    response = getattr(requests, verb)(host + endpoint, auth=auth_tup, json=body)
    print('{} {}'.format(response.status_code, response.reason))
    data = response.json()
    if isinstance(data, dict):
        messages = data.pop("messages", None)  
        if messages:
            print((json.dumps(messages, indent=4)))
    try:
        response.raise_for_status()
    except Exception as e:
        print(e)
        return None
    else:
        return data

# # In-process orders
# print('GET /api/v1/list-orders')
# filters = {"status": ["ordered"]}
# orders = espa_api('list-orders', body=filters)
# # status = espa_api('order-status', body=filters)
# print orders
# # print status

## View orders
# filters = {"status": ["complete", "ordered", "orderid"]}  # Here, we ignore any purged orders
# resp = espa_api('list-orders', body=filters)
# print((json.dumps(resp, indent=4)))

# # General Interactions: Authentication
# print('GET /api/v1/user')
# resp = espa_api('user')
# print((json.dumps(resp, indent=4)))


#************************************************
# Check availability, this information may be needed to modify 
# the remove list around line 150.
# Make a dictionay item to put the list in 
avail_list = {}
avail_list['inputs'] = FilteredImages
# print avail_list
resp = espa_api('available-products', body=avail_list)
print(json.dumps(resp, indent=4))

# sys.exit()


################################## SHUT OFF

l8_ls = set(FilteredImages)
# print l8_ls


avail_list = {}
avail_list['inputs'] = FilteredImages


# Differing products
# l8_prods = ['sr_evi', 'sr_ndvi', 'sr_savi', 'sr_ndmi', 'pixel_qa']
l8_prods = ['sr_ndvi', 'pixel_qa']
# l8_prods = ['sr_ndvi', 'sr_savi', 'sr_ndmi', 'pixel_qa']

# epsg: 32614 * What we want maybe.... 5/10/20

# Standard Albers CONUS
projection = {'aea': {'standard_parallel_1': 29.5,
                      'standard_parallel_2': 45.5,
                      'central_meridian': -96.0,
                      'latitude_of_origin': 23.0,
                      'false_easting': 0,
                      'false_northing': 0,
                      'datum': 'nad83'}}

# Let available-products place the acquisitions under their respective sensors

order = espa_api('available-products', body=avail_list)
# print((json.dumps(order, indent=4)))

# sys.exit()


# **NOTE**: Here we will not need to know what the sensor names were 
# for the Product IDs, thanks to the response from this `available-products` resource. 

# Replace the available products that was returned with what we want
for sensor in order.keys():
    if isinstance(order[sensor], dict) and order[sensor].get('inputs'):
        # if set(l7_ls) & set(order[sensor]['inputs']):
        #     order[sensor]['products'] = l7_prods
        if set(l8_ls) & set(order[sensor]['inputs']):
            order[sensor]['products'] = l8_prods

# Add in the rest of the order information
order['projection'] = projection
order['format'] = 'gtiff'
order['resampling_method'] = 'cc'



# Place the order
# print('POST /api/v1/order')
resp = espa_api('order', verb='post', body=order)
# print((json.dumps(resp, indent=4)))

# If successful, we will get our order-id

orderid = resp['orderid']


############################# Additional code which could be used later:


# # Check the status of an order

# print(('GET /api/v1/order-status/{}'.format(orderid)))
# resp = espa_api('order-status/{}'.format(orderid))
# print((json.dumps(resp, indent=4)))

# # Now, we can check for any completed products, and get the download url's for completed scenes

# # In[24]:

# print(('GET /api/v1/item-status/{0}'.format(orderid)))
# resp = espa_api('item-status/{0}'.format(orderid), body={'status': 'complete'})
# print((json.dumps(resp[orderid], indent=4)))


# # In[25]:

# # Once the order is completed or partially completed, can get the download url's
# for item in resp[orderid]:
#     print(("URL: {0}".format(item.get('product_dload_url'))))



# # In-process orders
# print('GET /api/v1/list-orders')
# filters = {"status": ["ordered"]}
# orders = espa_api('list-orders', body=filters)
# status = espa_api('order-status', body=filters)
# print orders
# print status

# # Find previous orders 

# # List backlog orders for the authenticated user.

# print('GET /api/v1/list-orders')
# filters = {"status": ["complete", "ordered", "orderid"]}  # Here, we ignore any purged orders
# resp = espa_api('list-orders', body=filters)
# print((json.dumps(resp, indent=4)))
# orderid = 1


# Once the order is completed or partially completed, can get the download url's
# item = 'http://espa.cr.usgs.gov/ordering/order-status/espa-201306381@panthers.greenville.edu-05092020-174923-809'
# print(("URL: {0}".format(item.get('product_dload_url'))))
