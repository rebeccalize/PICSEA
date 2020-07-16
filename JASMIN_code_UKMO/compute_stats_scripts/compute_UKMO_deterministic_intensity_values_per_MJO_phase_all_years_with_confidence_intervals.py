import track_statistics_library as ts
import picsea_library as pl
import sys
import os
import fnmatch
import numpy as np
from netCDF4 import Dataset, num2date, netcdftime
import datetime
import scipy.stats
import math

#This script, for the given set of TC seasons:
#Goes through each forecast of each storm in the SWIO, 
#Checks the amplitude of the given phase of the MJO on the date each forecast is initialised
#If the amplitude is >= 1.0, it includes this forecast in the sample
#It calculates the track error, wind speed error and MSLP error compared to ibtracs and/or ECMWF op analysis, 
#for each forecast and adds this to an array holding the full sample across all the seasons for this MJO phase
#Then calculates the mean across all the errors of the included forecasts for this MJO phase, and the 95% confidence intervals
#Saves these to a text file

year1s=[2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019] #
year2s=[2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020] # 


MJO1=int(sys.argv[1]) #MJO phase / phase pair [12,23,34,45,56,67,78,81]
MJO2=int(sys.argv[2])



#function for computing the mean and the confidence intervals of the data (omitting nan values from the calculations)
def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.nanmean(a), scipy.stats.sem(a, nan_policy='omit')
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h


array_len=31 #number of forecast timesteps (7 days x 4 forecasts)

#CURRENTLY RUNNING OTHER SH STORMS, NOT SIO - CHECK BEFORE RUNNING

#FOR DIFFERENT SEASONS, or across like 2006 - 2016, use this:

savedir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_and_intensity_errors/MJO_phase/"

#FOR LANDFALLING CYCLONES IN DIFFERENT COUNTRIES, USE THESE:
#datadir="/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_storm/mozambique_landfalling/" #SIO_storms/"
#savedir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_and_intensity_errors/mozambique_landfalling/"

if not os.path.exists(savedir):
    os.makedirs(savedir)



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

#print MJOamp[400:420]
#print MJOdates[400:420]



#Do the calculations comparing the forecasts to both the analysis track and the observed track
for obs_track in ["ibtracs"]:#, ["analysis","ibtracs"]:
	print obs_track

	#empty arrays to hold the average errors at each lead time, across all of the storms
	#this_mjo_trerr, this_mjo_mslp_bias, this_mjo_wind_bias, all_storms_trerr_sum, all_storms_trerr_wgt, all_storms_mslp_bias_sum, all_storms_mslp_bias_wgt, all_storms_wind_bias_sum, all_storms_wind_bias_wgt = (np.zeros(array_len) for i in range(9))
	this_mjo_obs_mslp, this_mjo_obs_wind, this_mjo_fcst_wind, this_mjo_fcst_mslp = (np.zeros(array_len) for i in range(4))
	
	for y1, y2 in zip(year1s, year2s):
	
		datadir="/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_storm/"+str(y1)+"_"+str(y2)+"/SIO_storms/" #SIO_storms/"
	 
	 	#Find out the number of storms in this TC season and create array of the IDs so we can loop over them
		season_dirs=[]
		for root,dirs,files in os.walk(datadir):
			for dir in dirs:
				season_dirs.append(dir)
		NS = len(season_dirs) #total number of storms in this season in the SIO
	
		#Calculate the error of each foreacst at each lead time, for each storm in this TC season / directory, and save
		#Then average across all the forecasts of each storm, and save
		def storm_stats(ibtracs_file, analysis_file, nwp_files_list, use_analysis_or_ibtracs, array_len):
			"""This function calculates the statistics for all the forecasts of one storm"""
			
			global this_mjo_obs_mslp
			global this_mjo_obs_wind
			global this_mjo_fcst_mslp
			global this_mjo_fcst_wind

			if use_analysis_or_ibtracs == "analysis":
				obs_file = analysis_file
			elif use_analysis_or_ibtracs == "ibtracs":
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
			#print NT
			for ff,a in zip(nwp_files_list, range(NT)):
				#print ff

				#get the forecast date, lat, lon and vorticity
				fcst_data=np.genfromtxt(ff,dtype=float,skip_header=1)
				fcst_lon=fcst_data[:,7]
				fcst_lat=fcst_data[:,8]
				fcst_mslp=fcst_data[:,9]
				fcst_wind=fcst_data[:,10]
				fcst_datelist = pl.get_dates(fcst_data)
			
				#check whether first date of this forecast file has an MJO amplitude > 1 (do we want to use it for this MJO phase?)
				#first get the forecast data date into datetime format to compare to the MJO dates
				first_fcst_date = datetime.datetime.strptime(str(fcst_datelist[0]), "%Y%m%d%H")	
			
			
				#print first_fcst_date.strftime("%m-%d")
			
				#then find the index of MJO dates at which the date is the first date of the forecast
				if first_fcst_date.strftime("%m-%d") == '02-29': #leap years missing in the MJO data?
			
					continue
				else:		
					z = MJOdates.index(first_fcst_date.strftime("%Y-%m-%d"))
			
					#then find out the MJO amplitude on this date (i.e. at the index of this date)
					amp = MJOamp[z]
					#print "MJO amp", amp
			
					#if the MJO amplitude is less than 1, don't include this foreacst
					#we only want to compute the errors for forecasts started during this MJO phase / phase pair, when the MJO amplitude is >1
					
				
					#if the MJO amplitude was >1 on this day, include this forecast in the stats calculations
					if (MJOphase[z] == MJO1 or MJOphase[z] == MJO2) and MJOamp[z] >= 1.0:			

						"""We need to get the indices of both the observed data and the forecast data, where the dates match"""
						#This is because the dates in the observed track file and forecast track files cover different ranges,
						#depending on the date the forecast was initialised and the period the forecast covers

						#find the indices of the forecast array, where the dates exist in the observed dates array
						indices_fcst = np.nonzero(np.in1d(fcst_datelist,obs_datelist))[0]
					
						# find the indices of the observed array, where the dates exist in the forecast dates array
						indices_obs = np.nonzero(np.in1d(obs_datelist, fcst_datelist))[0]
	
						#So the first few lead times of the forecast might have no observations to match against
						#BUT we don't want the first matching date to then be calculated as if it were lead time 1 of the forecast
						#So fill the arrays with empty values for the first few lead times where there are no obs to verify against
						#How many lead times do we need to skip? Up to the index of the first matched forecast date:

						#if using ibtracs, sometimes there are no observations for the entire length of the forecast track
						#so we tell it not to run the calculations if there were no matching dates:
						if not len(indices_fcst) == 0:

							lt_to_skip = indices_fcst[0]


							#empty arrays to hold the data at the matched timesteps
							matched_data_dates, matched_obs_lon_lat, matched_fcst_lon_lat = (np.ma.masked_all((array_len,2)) for i in range(3))
							matched_obs_mslp, matched_fcst_mslp, matched_obs_wind, matched_fcst_wind = (np.ma.masked_all((array_len, 1)) for i in range(4))

							#get the observed track data for the dates where this forecast track matches the observed track
							#at each lead time
							#print "indices_obs: ", indices_obs
							for i,z in zip(indices_obs,range(array_len-lt_to_skip)):
								matched_data_dates[z+lt_to_skip, 0] = obs_datelist[i]
								matched_obs_lon_lat[z+lt_to_skip,0] = obs_lon[i]
								matched_obs_lon_lat[z+lt_to_skip,1] = obs_lat[i]
								matched_obs_mslp[z+lt_to_skip] = obs_mslp[i]
								matched_obs_wind[z + lt_to_skip] = obs_wind[i]

							#get the forecast track data for the dates where this forecast track matches the observed track, at each lead time

							#print "indices_fcst: ", indices_fcst
							#print "indices_fcst[0:array_len-1]: ", indices_fcst[0:array_len-1]
							for i,z in zip(indices_fcst[0:array_len-1], range(array_len-lt_to_skip)):
								matched_data_dates[z+lt_to_skip,1] = fcst_datelist[i]
								matched_fcst_lon_lat[z+lt_to_skip, 0] = fcst_lon[i]
								matched_fcst_lon_lat[z+lt_to_skip, 1] = fcst_lat[i]
								matched_fcst_mslp[z+lt_to_skip] = fcst_mslp[i]
								matched_fcst_wind[z + lt_to_skip] = fcst_wind[i]
							#print matched_data_dates

							#calculate the track error (great circle distance, in km) for this forecast, at each lead time
							#calculate the mslp and wind biases for this forecast, at each lead time
							#trerr, mslp_bias, wind_bias = (np.zeros(array_len) for i in range(3))
							obs_wind_arr, obs_mslp_arr, fcst_wind_arr, fcst_mslp_arr = (np.zeros(array_len) for i in range(4))
						
							for lt in range(array_len):

								#trerr[lt]=ts.trerr(matched_obs_lon_lat[lt,:],matched_fcst_lon_lat[lt,:])
								#mslp_bias[lt]=ts.bias(matched_obs_mslp[lt],matched_fcst_mslp[lt])
								#wind_bias[lt]=ts.bias(matched_obs_wind[lt],matched_fcst_wind[lt])
								if matched_obs_wind[lt] > 10000:
				
									obs_wind_arr[lt] = np.nan
								else:
									obs_wind_arr[lt] = matched_obs_wind[lt]
								
								if matched_obs_mslp[lt] > 10000:
									obs_mslp_arr[lt] = np.nan
								else:
									obs_mslp_arr[lt] = matched_obs_mslp[lt]
							
								if matched_fcst_wind[lt] > 10000:
									print matched_fcst_wind[lt]
									fcst_wind_arr[lt] = np.nan
								else:
									fcst_wind_arr[lt] = matched_fcst_wind[lt]
								
								if matched_fcst_mslp[lt] > 10000:
									print matched_fcst_mslp[lt]
									fcst_mslp_arr[lt] = np.nan
								else:
									fcst_mslp_arr[lt] = matched_fcst_mslp[lt]
							

							#print trerr
							#print mslp_bias
							#print wind_bias
						
							print obs_wind_arr
							print fcst_wind_arr
							print obs_mslp_arr
							print fcst_mslp_arr
						
							#We want to know the sample size of forecats included, so for each lead time, as long as there is a value
							#(i.e. don't count one in the sample size if we couldn't verify the first 3 days of the forecast because there
							#wasn't yet a best track to verify against), add one to the sample size (weight)
						
							#for lt in range(array_len):
								#if not np.isnan(trerr[lt]):
									#all_storms_trerr_wgt[lt] += 1
								
							#If the array is jsut full of zeros, this is the first forecast we've run calcs for, so replace the empty array
							#with the errors for this forecast
							if np.all(this_mjo_obs_wind == 0):
								for lt in range(array_len):
							
									this_mjo_obs_wind[lt] = obs_wind_arr[lt]
									this_mjo_obs_mslp[lt] = obs_mslp_arr[lt]
									this_mjo_fcst_wind[lt] = fcst_wind_arr[lt]
									this_mjo_fcst_mslp[lt] = fcst_mslp_arr[lt]
								
									#this_mjo_trerr[lt] = trerr[lt]
									#this_mjo_mslp_bias[lt] = mslp_bias[lt]
									#this_mjo_wind_bias[lt] = wind_bias[lt]
								
							#Otherwise, if there are already values in that array, this isn't the first forecast we've run calcs for, 
							#so stack these errors on top of the existing array, so that we can calculate the mean later over all the forecasts
							#initialised in this MJO phase
				
							else:
								this_mjo_obs_wind = np.vstack([this_mjo_obs_wind, obs_wind_arr])
								this_mjo_obs_mslp = np.vstack([this_mjo_obs_mslp, obs_mslp_arr])
								this_mjo_fcst_wind = np.vstack([this_mjo_fcst_wind, fcst_wind_arr])
								this_mjo_fcst_mslp = np.vstack([this_mjo_fcst_mslp, fcst_mslp_arr])
							
								#this_mjo_trerr = np.vstack([this_mjo_trerr, trerr])
								#this_mjo_mslp_bias = np.vstack([this_mjo_mslp_bias, mslp_bias])
								#this_mjo_wind_bias = np.vstack([this_mjo_wind_bias, wind_bias])

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
			pattern="ukmo_nwp*.txt"
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
			storm_stats(ibtracs_file, analysis_file, nwp_files_list,obs_track,array_len)


	print "overall arrays: ###################################"
	#print this_mjo_trerr
	#print this_mjo_mslp_bias
	#print this_mjo_wind_bias
	#calculate the average errors at each lead time, across all the storms in this TC season / directory
	
	#avg_trerr_all_storms, avg_mslp_bias_all_storms, avg_wind_bias_all_storms = (np.zeros((array_len, 3)) for i in range(3))
	
	avg_obs_wind, avg_obs_mslp, avg_fcst_wind, avg_fcst_mslp = (np.zeros((array_len, 3)) for i in range(4))
	
	
	for lt in range(array_len):
	
		#This would output the mean, the min bound for 95% confidence interval, and the max bound for 95% confidence interval, in that order
			
		avg_obs_wind[lt,0], avg_obs_wind[lt,1], avg_obs_wind[lt,2] = mean_confidence_interval(this_mjo_obs_wind[:,lt])
		avg_obs_mslp[lt,0], avg_obs_mslp[lt,1], avg_obs_mslp[lt,2] = mean_confidence_interval(this_mjo_obs_mslp[:,lt])
		
		avg_fcst_wind[lt,0], avg_fcst_wind[lt,1], avg_fcst_wind[lt,2] = mean_confidence_interval(this_mjo_fcst_wind[:,lt])
		avg_fcst_mslp[lt,0], avg_fcst_mslp[lt,1], avg_fcst_mslp[lt,2] = mean_confidence_interval(this_mjo_fcst_mslp[:,lt])
		


	np.savetxt(savedir + str(year1s[0])+"_"+str(year2s[-1])+"_MJO_phase"+str(MJO1)+str(MJO2)+"_average_"+obs_track+"_wind_speed_values_per_lead_time.with_confidence_intervals.txt", avg_obs_wind, '%.4f')
	np.savetxt(savedir + str(year1s[0])+"_"+str(year2s[-1])+"_MJO_phase"+str(MJO1)+str(MJO2)+"_average_"+obs_track+"_mslp_values_per_lead_time.with_confidence_intervals.txt", avg_obs_mslp, '%.4f')
	
	np.savetxt(savedir + str(year1s[0])+"_"+str(year2s[-1])+"_MJO_phase"+str(MJO1)+str(MJO2)+"_average_UKMO_deterministic_wind_speed_values_per_lead_time_vs_"+obs_track+".with_confidence_intervals.txt", avg_fcst_wind, '%.4f')
	np.savetxt(savedir + str(year1s[0])+"_"+str(year2s[-1])+"_MJO_phase"+str(MJO1)+str(MJO2)+"_average_UKMO_deterministic_mslp_values_per_lead_time_vs_"+obs_track+".with_confidence_intervals.txt", avg_fcst_mslp, '%.4f')

	
	
	
	
	
