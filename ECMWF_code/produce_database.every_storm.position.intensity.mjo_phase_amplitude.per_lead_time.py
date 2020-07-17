import track_statistics_library as ts
import picsea_library as pl
import sys
import os
import fnmatch
import numpy as np
from netCDF4 import Dataset, num2date, netcdftime
import datetime
import scipy.stats

#This script, for the UKMO deterministic forecast, loops over every storm in the specified seasons, and, for each lead time days 0 - 7, saves a file containing
#For every individual forecast, the observed location, wind speed and pressure of the storm (from ibtracs) on the forecats validation date, 
#the forecast initialisation date, the forecast location, wind speed and pressure, the location error, 
#and the MJO phase and amplitude on the date the forecast was initialised
#wind speeds saved in km/h

fcst_type = 'ctrl'


year1s=[2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019]
year2s=[2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020]

#res1
#year1s=[2010, 2011, 2012, 2013]
#year2s=[2011, 2012, 2013, 2014]

#res2
#year1s = [2014, 2015, 2016]
#year2s = [2015, 2016, 2017]

#res3
#year1s = [2017, 2018, 2019]
#year2s = [2018, 2019, 2020]


array_len=31 #number of forecast timesteps (7 days x 4 forecasts)


savedir = "/home/users/emerton/analysis_scripts/"


if not os.path.exists(savedir):
    os.makedirs(savedir)



#function for computing the mean and the confidence intervals of the data (omitting nan values from the calculations)
def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.nanmean(a), scipy.stats.sem(a, nan_policy='omit')
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h


#This MJO file just has the phase and amplitude on every day, not split up by phase, probably easiest to use for any analysis with individual phases
MJO_file = "/gws/nopw/j04/klingaman/datasets/MJO_INDICES/MJO_indices_rmm1_rmm2.jan1979-may2020_dmeans_ts.index_values.nc"

#MJO_file = "/gws/nopw/j04/klingaman/datasets/MJO_INDICES/MJO_rmm1_rmm2.jan-dec_dmeans_ts.1979-2019.nc"

print MJO_file

ffMJO = Dataset(MJO_file,'r')
MJOamp = ffMJO.variables['amplitude'][:]
MJOphase = ffMJO.variables['phase'][:]


start_time = datetime.datetime(1979,1,1)
MJO_tvalue = np.array([start_time + datetime.timedelta(days=i) for i in xrange(len(MJOamp))])
			
MJOdates = [j.strftime("%Y-%m-%d") for j in MJO_tvalue]


#Do the calculations comparing the forecasts to both the analysis track and the observed track
for lead_time in [0,1,2,3,4,5,6,7]:#, ["analysis","ibtracs"]:

	#'this category' actually refers to all storms across all years, for this obs_track, I just don't want to rename them...
	#empty arrays to hold the average errors at each lead time, across all of the storms and seasons, for this intensity category, vs ibtracs
	#this_category_trerr, this_category_mslp_bias, this_category_wind_bias, all_storms_trerr_sum, all_storms_trerr_wgt, all_storms_mslp_bias_sum, all_storms_mslp_bias_wgt, all_storms_wind_bias_sum, all_storms_wind_bias_wgt = (np.zeros(array_len) for i in range(9))


	every_TC_database = np.zeros(15) #Think I can just produce the number of columns and then use dstack?
	#every_TC_database[0] = "TC_season"
	#every_TC_database[1] = "track_number"
	#every_TC_database[2] = "obs_date"
	#every_TC_database[3] = "obs_lon"
	#every_TC_database[4] = "obs_lat"
	#every_TC_database[5] = "obs_wind_kmh"
	#every_TC_database[6] = "obs_mslp"
	#every_TC_database[7] = "fcst_date"
	#every_TC_database[8] = "fcst_lon"
	#every_TC_database[9]  = "fcst_lat"
	#every_TC_database[10] = "track_error"
	#every_TC_database[11] = "fcst_wind_kmh"
	#every_TC_database[12] = "fcst_mslp"
	#every_TC_database[13] = "MJO_phase"
	#every_TC_database[14] = "MJO_amplitude"
	

		
	for y1, y2 in zip(year1s,year2s):
	
		if fcst_type == 'hres':
		
			datadir="/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_storm/"+str(y1)+"_"+str(y2)+"/SIO_storms/"
			
		else:
			datadir="/gws/nopw/j04/klingaman/emerton/ukmo_TIGGE_analysis/reformatted_track_files_per_storm/Y"+str(y1)+str(y2)+"/SIO_storms/"
			
		#Find out the number of storms in this TC season and create array of the IDs so we can loop over them
		season_dirs=[]
		for root,dirs,files in os.walk(datadir):
			for dir in dirs:
					season_dirs.append(dir)
		NS = len(season_dirs) #total number of storms in this season in the SIO
	
		
		#Calculate the error of each foreacst at each lead time, for each storm in this TC season / directory, and save
		#Then average across all the forecasts of each storm, and save
		def storm_stats(ibtracs_file, analysis_file, nwp_files_list, array_len, track_number):
		
			global every_TC_database
			
			
			
			"""This function calculates the statistics for all the forecasts of one storm"""

			
			obs_file = ibtracs_file
			#Find out number of forecast tracks for this storm
			
			NT = len(nwp_files_list)


			#Get the date, lon, lat and vorticity data for the observed track
			obs_data=np.genfromtxt(obs_file, dtype=float, skip_header=1)
			obs_lon=obs_data[:,7]
			obs_lat=obs_data[:,8]
			obs_mslp=obs_data[:,9]
			obs_wind=obs_data[:,10]
			obs_datelist = pl.get_dates(obs_data)

			#compute the statistics for each forecast of this storm, at each lead time
			"""need to make sure we compare the forecast timestamp with the correct obs timestamp!!"""
			
				
			for ff,a in zip(nwp_files_list, range(NT)):
				print ff

				#get the forecast date, lat, lon and vorticity
				fcst_data=np.genfromtxt(ff,dtype=float,skip_header=1)
				fcst_lon=fcst_data[:,7]
				fcst_lat=fcst_data[:,8]
				fcst_mslp=fcst_data[:,9]
				fcst_wind=fcst_data[:,10]
				fcst_datelist = pl.get_dates(fcst_data)
			
				
				first_fcst_date = datetime.datetime.strptime(str(fcst_datelist[0]), "%Y%m%d%H")	
				
				#print fcst_datelist
				#print obs_datelist
				
				if first_fcst_date.strftime("%m-%d") == '02-29': #leap years missing in the MJO data?
			
					continue
				else:		
					z = MJOdates.index(first_fcst_date.strftime("%Y-%m-%d"))
			
			
				if len(fcst_datelist) <= lead_time*4:
					continue
					
				#print lead_time*4
				#print len(fcst_datelist)
				
				elif fcst_datelist[lead_time*4] in obs_datelist[:]:

						
					this_dates_array = np.zeros(15)
					
					this_dates_array[0] = int(y1)
					
					print track_number
					print track_number[2:6]
					
					this_dates_array[1] = int(track_number[2:6])
					
					
					this_dates_array[2] = int(fcst_datelist[lead_time*4]) #validation date / obs date - date we're forecasting for
					
					i = obs_datelist.index(fcst_datelist[lead_time*4])
					print i
					
					this_dates_array[3] = obs_lon[i]
					this_dates_array[4] = obs_lat[i]
					
					if obs_wind[i] > 10000:
						this_dates_array[5] = np.nan
					else:
						this_dates_array[5] = obs_wind[i] *3.6
						
					if obs_mslp[i] > 10000:
						this_dates_array[6] = np.nan
					else:
						this_dates_array[6] = obs_mslp[i]
					
			
					this_dates_array[7] = int(fcst_datelist[0]) #date forecast initialised
					
					this_dates_array[8] = fcst_lon[lead_time*4]
					this_dates_array[9] = fcst_lat[lead_time*4]
					this_dates_array[11] = fcst_wind[lead_time*4]*3.6
					this_dates_array[12] = fcst_mslp[lead_time*4]
					
					if MJOamp[z] >= 1.0:
						this_dates_array[13] = int(MJOphase[z])
						this_dates_array[14] = MJOamp[z]
					
					else:
						this_dates_array[13] = np.nan
						this_dates_array[14] = np.nan
	
					
					this_dates_array[10] = ts.trerr([obs_lon[i], obs_lat[i]], [fcst_lon[lead_time*4], fcst_lat[lead_time*4]])
					
					print this_dates_array
					
					every_TC_database = np.vstack([every_TC_database, this_dates_array])
					
					print every_TC_database
								
							
									
				else:
					continue			
	
			#############################################################################
			###################### CALL THE FUNCTION ####################################
			#############################################################################
	
		#loop over storm directories in this TC season, and in each directory find each file that matches the analysis, ibtracs and nwp filenames
		#then call the function for each storm using the identified files
		for dir, x in zip(season_dirs,range(len(season_dirs))):
			print datadir+dir
			nwp_files_list = []
			analysis_file=0
			ibtracs_file=0
			#make a list of all the files in this directory
			list_of_all_files = os.listdir(datadir+dir)
			#print list_of_all_files
			
			if fcst_type == 'hres':
			
				pattern="ukmo_nwp*.txt"
				
			elif fcst_type == 'mean':
				
				pattern="ukmo_tigge_mean*.txt"
				
			elif fcst_type == 'ctrl':
				
				pattern="ukmo_tigge_ctrl*.txt"
				
			pattern2="analysis_*.txt"
			pattern3="ibtracs_*.txt"
			#for each file in the directory, find the ones that match the forecast filename, analysis filename and ibtracs filename
			for entry in list_of_all_files:
				if fnmatch.fnmatch(entry,pattern):
					nwp_files_list.append(datadir+dir+"/"+entry)
				elif fnmatch.fnmatch(entry,pattern2):
					analysis_file = datadir+dir+"/"+entry
					#print analysis_file
				elif fnmatch.fnmatch(entry,pattern3):
					ibtracs_file = datadir+dir+"/"+entry
			#run the statistics for this storm
			storm_stats(ibtracs_file, analysis_file, nwp_files_list,array_len, dir)


	
	with open(savedir + str(year1s[0])+"_"+str(year2s[-1])+".ukmo_"+fcst_type+".database_of_each_individual_forecast.obs_and_fcst.location_intensity_trackerror_MJOphase.lead_time"+str(lead_time)+".using_ibtracs.txt", 'wb') as f:
	
		f.write('TC_season obs_validation_date obs_lon obs_lat obs_wind_kmh obs_mslp fcst_init_date fcst_lon fcst_lat track_error fcst_wind_kmh fcst_mslp MJO_phase MJO_amplitude\n')
	
		np.savetxt(f,  every_TC_database[1:,:], '%.4f')
	

	
	
	
	
