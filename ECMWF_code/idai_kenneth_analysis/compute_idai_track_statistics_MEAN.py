import track_statistics_library as ts
import picsea_library as pl
import sys
import os
import fnmatch
import numpy as np

#This script, for one TC season (e.g. 2015-2016) in the SIO:
#Uses functions in track_statistics_library.py, to:
#Compute the track location error for every consecutive deterministic forecast of each storm (great circle distance in km), at each lead time
#Compute the intensity bias (vorticity) for every consecutive deterministic forecast of each storm, at each lead time
#Compute the average track location error and intensity bias across all the forecast for each storm
#Compute the average track location error and intensity bias across all the forecasts and across all the storms in this TC season

array_len=29 #number of forecast timesteps (7 days x 4 forecasts)

#THIS SCRIPT PRODUCES TRACK ERRORS FOR IDAI ONLY!

datadir="/perm/mo/more/picsea/reformatted_idai_kenneth_tracks/"
savedir = "/perm/mo/more/picsea/idai_kenneth_track_statistics/"

if not os.path.exists(savedir):
    os.makedirs(savedir)
    



#Do the calculations comparing the forecasts to both the analysis track and the observed track
for obs_track in ["analysis", "ibtracs"]:#, ["analysis","ibtracs"]:
	print obs_track

	
	#Calculate the error of each foreacst at each lead time, for each storm in this TC season / directory, and save
	#Then average across all the forecasts of each storm, and save
	def storm_stats(ibtracs_file, analysis_file, nwp_files_list, use_analysis_or_ibtracs, array_len):
		"""This function calculates the statistics for all the forecasts of one storm"""

		if use_analysis_or_ibtracs == "analysis":
			obs_file = analysis_file
		elif use_analysis_or_ibtracs == "ibtracs":
			obs_file = ibtracs_file
		#Find out number of forecast tracks for this storm
		NT = len(nwp_files_list)

		#empty arrays to hold the track error and intensity bias statistics for this storm, for each forecast track
		storm_trerr, storm_mslp_bias, storm_wind_bias =(np.zeros((NT,array_len)) for i in range(3))

		#empty arrays to hold the error sums and counts, for calculating the average errors for this storm
		storm_trerr_sum, storm_trerr_wgt, storm_mslp_bias_sum, storm_mslp_bias_wgt, storm_wind_bias_sum, storm_wind_bias_wgt = (np.zeros(array_len) for i in range(6))

		#Get the date, lon, lat and vorticity data for the observed track
		obs_data=np.genfromtxt(obs_file, dtype=float, skip_header=1,usecols=np.arange(0,11))
		obs_lon=obs_data[:,7]
		obs_lat=obs_data[:,8]
		obs_mslp=obs_data[:,9]
		obs_wind=obs_data[:,10]
		obs_datelist = pl.get_dates(obs_data)
		
		print obs_lon
		print obs_lat

		#compute the statistics for each forecast of this storm, at each lead time
		"""need to make sure we compare the forecast timestamp with the correct obs timestamp!!"""
		print NT
		for ff,a in zip(nwp_files_list, range(NT)):
			#print ff

			#get the forecast date, lat, lon and vorticity
			fcst_data=np.genfromtxt(ff,dtype=float,skip_header=1,usecols=np.arange(0,11))
			fcst_lon=fcst_data[:,7]
			fcst_lat=fcst_data[:,8]
			fcst_mslp=fcst_data[:,9]
			fcst_wind=fcst_data[:,10]
			fcst_datelist = pl.get_dates(fcst_data)
			
			print fcst_lon
			print fcst_lat
			
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
				for i,z in zip(indices_obs,range(array_len-lt_to_skip)):
					matched_data_dates[z+lt_to_skip, 0] = obs_datelist[i]
					matched_obs_lon_lat[z+lt_to_skip,0] = obs_lon[i]
					matched_obs_lon_lat[z+lt_to_skip,1] = obs_lat[i]
					matched_obs_mslp[z+lt_to_skip] = obs_mslp[i]
					matched_obs_wind[z + lt_to_skip] = obs_wind[i]

				#get the forecast track data for the dates where this forecast track matches the observed track, at each lead time


				for i,z in zip(indices_fcst[0:array_len-1], range(array_len-lt_to_skip)):
					matched_data_dates[z+lt_to_skip,1] = fcst_datelist[i]
					matched_fcst_lon_lat[z+lt_to_skip, 0] = fcst_lon[i]
					matched_fcst_lon_lat[z+lt_to_skip, 1] = fcst_lat[i]
					matched_fcst_mslp[z+lt_to_skip] = fcst_mslp[i]
					matched_fcst_wind[z + lt_to_skip] = fcst_wind[i]
				print matched_data_dates

				#calculate the track error (great circle distance, in km) for this forecast, at each lead time
				trerr, mslp_bias, wind_bias = (np.zeros(array_len) for i in range(3))
				for lt in range(array_len):

					print matched_obs_lon_lat[lt,:], matched_fcst_lon_lat[lt,:]
					trerr[lt]=ts.trerr(matched_obs_lon_lat[lt,:],matched_fcst_lon_lat[lt,:])
					#print trerr
					mslp_bias[lt]=ts.bias(matched_obs_mslp[lt],matched_fcst_mslp[lt])
					wind_bias[lt]=ts.bias(matched_obs_wind[lt],matched_fcst_wind[lt])

				#add the errors for this forecast track, to the arrays holding all the errors for this storm
				#and add one to the "weight" for each error, which counts the number of forecasts contributing to the error calculation
				#this is because some forecasts are shorter than others, and we want to divide by the correct sample size

				#print wind_bias
				for lt in range(array_len):
					if not np.isnan(trerr[lt]):
						storm_trerr[a,lt] = trerr[lt]
						storm_trerr_sum[lt] += trerr[lt]
						storm_trerr_wgt[lt] += 1
				
					if not np.isnan(mslp_bias[lt]):
						storm_mslp_bias[a,lt] = mslp_bias[lt]
						storm_mslp_bias_sum[lt] += mslp_bias[lt]
						storm_mslp_bias_wgt[lt] += 1
				
					if not np.isnan(wind_bias[lt]):
						storm_wind_bias[a,lt] = wind_bias[lt]
						storm_wind_bias_sum[lt] += wind_bias[lt]
						storm_wind_bias_wgt[lt] += 1
					

		#calculate the average error at each lead time, across all the forecasts of this storm
		storm_trerr_mean, storm_mslp_bias_mean, storm_wind_bias_mean = (np.zeros(array_len) for i in range(3))
		for lt in range(array_len):
			storm_trerr_mean[lt] = storm_trerr_sum[lt] / storm_trerr_wgt[lt]
			storm_mslp_bias_mean[lt] = storm_mslp_bias_sum[lt] / storm_mslp_bias_wgt[lt]
			storm_wind_bias_mean[lt] = storm_wind_bias_sum[lt] / storm_wind_bias_wgt[lt]

		#print storm_wind_bias_sum
		#print storm_wind_bias_wgt
		#print storm_wind_bias_mean

		np.savetxt(savedir + "idai_MEAN_each_forecast_location_error_per_lead_time_vs_"+use_analysis_or_ibtracs+".txt", storm_trerr[:,:], '%.4f')
		np.savetxt(savedir + "idai_MEAN_each_forecast_mslp_bias_per_lead_time_vs_"+use_analysis_or_ibtracs+".txt", storm_mslp_bias[:,:], '%.4f')
		np.savetxt(savedir + "idai_MEAN_each_forecast_wind_bias_per_lead_time_vs_" + use_analysis_or_ibtracs + ".txt", storm_wind_bias[:, :], '%.4f')

		np.savetxt(savedir + "idai_MEAN_average_location_error_per_lead_time_vs_"+use_analysis_or_ibtracs+".txt", storm_trerr_mean[:], '%.4f')
		np.savetxt(savedir + "idai_MEAN_average_mslp_bias_per_lead_time_vs_"+use_analysis_or_ibtracs+".txt", storm_mslp_bias_mean[:], '%.4f')
		np.savetxt(savedir + "idai_MEAN_average_wind_bias_per_lead_time_vs_" + use_analysis_or_ibtracs + ".txt", storm_wind_bias_mean[:], '%.4f')


	### CALL THE FUNCTION ###
	#loop over storm directories in this TC season, and in each directory find each file that matches the analysis, ibtracs and nwp filenames
	#then call the function for each storm using the identified files
	print datadir
	nwp_files_list = []
	analysis_file=0
	ibtracs_file=0
	#make a list of all the files in this directory
	list_of_all_files = os.listdir(datadir)
	pattern="*_idai_ECMWF_mean_3way_matched_reformatted.txt"
	pattern2="idai_analysis_*.txt"
	pattern3="idai_ibtracs_*.txt"
	#for each file in the directory, find the ones that match the forecast filename, analysis filename and ibtracs filename
	for entry in list_of_all_files:
		if fnmatch.fnmatch(entry,pattern):
			nwp_files_list.append(datadir+entry)
		elif fnmatch.fnmatch(entry,pattern2):
			analysis_file = datadir+"idai_analysis_track_reformatted.txt"
			#print analysis_file
		elif fnmatch.fnmatch(entry,pattern3):
			ibtracs_file = datadir+"idai_ibtracs_reformatted.txt"
	#run the statistics for this storm
	print nwp_files_list
	storm_stats(ibtracs_file, analysis_file, nwp_files_list,obs_track,array_len)
	
	
	
	
	


