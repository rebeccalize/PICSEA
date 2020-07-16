import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import numpy as np
from calendar import monthrange
from datetime import date
import operator
import math


	
#read in observation data
obs_data = np.genfromtxt("/home/users/emerton/MADAGASCAR_RR_daily_2006_2018_CSV.csv", delimiter=',',dtype=str)
	
station_names = obs_data[0,1:]
station_IDs = obs_data[1,1:]
station_lons = obs_data[2,1:]
station_lats = obs_data[3,1:]

obs_precip = obs_data[4:,:] #first value of each row [:,0] is the date! each column [x,:] is a different station

#set up a dictionary to hold the obs and trmm precip at each station location / grid point:
pcp_dict = {}
pcp_dict["obs"] = {}
pcp_dict["trmm"] = {}

for ID in station_IDs:
	pcp_dict["obs"][ID] = []
	pcp_dict["trmm"][ID] = []



no_stations = len(station_IDs)

#need to include something to exclude any negative rainfall obs (think missing value uses -99)


trmm_dir="/gws/nopw/j04/klingaman/datasets/TRMM_3B42/V7_NC_daily/"

for s,ID,name in zip(range(1,no_stations+1), station_IDs, station_names):



	for year in [2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018]: #,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018
	
		#read in the daily TRMM data for this year
	
		trmm_file = trmm_dir+"3B42_daily."+str(year)+".nc"
	
		fftrmm = Dataset(trmm_file,'r')
		#trmm_pcp = fftrmm.variables['precipitation'][:]
	
		if year < 2016:		
			trmm_pcp_array = fftrmm.variables['r'][:]
		
		elif year >= 2016:
			trmm_pcp_array = fftrmm.variables['precipitation'][:]
		
		trmm_lons = fftrmm.variables['longitude'][:]
		trmm_lats = fftrmm.variables['latitude'][:]
	
	
		#for s in range(1,no_stations+1):
	
			#print obs_precip[:,s]
			#print obs_precip[:,0]
		
		
		d0 = date(year,1,1)
		
	
	
		for month in [1,2,3,4,5,6,7,8,9,10,11,12]: #,2,3,4,5,6,7,8,9,10,11,12
		
			days=monthrange(year,month)[1]
		
			for day in range(1,days+1): #1,days+1
		
				print year,month,day
				#d = day+1
			
				#get the index of this date in the observation data array
				i = np.where(obs_precip[:,0] == str(year)+str(month).zfill(2)+str(day).zfill(2))
			
				#testing:	
				#print 'looking for ',str(year)+str(month).zfill(2)+str(day).zfill(2)	
				#print 'index: ', i
				#print 'date found: ',obs_precip[i,0]
			
			
				#get index of this date in the TRMM files (no date variable, file starts at 1/1/year, so find how many days this day is from 1st Jan)		
				d1 = date(year,month,day)
				delta = d1-d0
				j=delta.days
				#print j
			
			
				#print "station ID: ", ID
			
				station_pcp = obs_precip[i,s][0]
				
				#not sure why this loop is necessary to convert to a float, but it's the only thing I've found that works...				
				for item in station_pcp: 
				
					a = float(item)
					#print "station pcp: ", a
				
					if a >= 0.0:
						pcp_dict["obs"][ID].append(a)
					
					else:
						print "missing value?: ", a
						pcp_dict["obs"][ID].append(np.ma.masked) #if it doesn't append something, the dates won't match up
						
						
						
				slat = float(station_lats[s-1])
				slon = float(station_lons[s-1])
				
				#get nearest lat and lon values in TRMM to the station lat and lon, and the index of these in the trmm lat and lon arrays
				tlat = min(range(len(trmm_lats)), key=lambda y: abs(trmm_lats[y]-slat))
				tlon = min(range(len(trmm_lons)), key=lambda x: abs(trmm_lons[x]-slon))
				
				#print slat, trmm_lats[tlat]
				#print slon, trmm_lons[tlon]
				
				#get the TRMM precip at this grid point, to match up with the observed precip at the station at that location
				trmm_pcp = trmm_pcp_array[j,tlat,tlon] #[time,lat,lon] is the trmm structure
				
				#print "trmm_pcp: ",trmm_pcp
				pcp_dict["trmm"][ID].append(trmm_pcp)
				
				
				
				
				#########################################################################
				# up to this point, the code has all been about reading in the data, and getting the trmm data at the grid point closest
				# to each rainfall station, on the same dates
				#########################################################################
						

	def roundup(x):
    		return int(math.ceil(x / 10.0)) * 10	
			
					
	fig, ax = plt.subplots()
	fig.set_size_inches(6,4)

	#dictmax=max(pcp_dict.iteritems(), key=operator.itemgetter(1))[0]
	#key, value = max(pcp_dict.iteritems(), key=lambda x:x[1])
	omax = max(pcp_dict["obs"][ID][:])
	tmax = max(pcp_dict["trmm"][ID][:])

	dictmax = roundup(max(omax,tmax))

	print dictmax

	print omax,tmax
	#print value

	

	plt.scatter(pcp_dict["obs"][ID][:], pcp_dict["trmm"][ID][:],s=5,color='navy',label=str(ID)+"\n"+name)	
	
	

	plt.xlim(0,dictmax)
	plt.xlabel('Observed Daily Rainfall (mm)', fontsize=9)

	plt.ylabel('TRMM Daily Rainfall (mm)', fontsize=9)
	plt.ylim(0, dictmax)
	
	plt.legend()
	
	
	plt.savefig("/home/users/emerton/madagascar_obs_rainfall_figures/madagascar_obs_vs_trmm_rainfall_scatter_station_"+str(ID)+"_"+name+".png", dpi=400,bbox_inches='tight')	
	print "saved madagascar_obs_vs_trmm_rainfall_scatter_station_"+str(ID)+"_"+name+".png"	
					
					
#for ID in station_IDs:	
	#print ID	
	#print pcp_dict["obs"][ID]
	#print pcp_dict["trmm"][ID]


				
			
	
	



