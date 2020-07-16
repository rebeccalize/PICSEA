import track_statistics_library as ts
import picsea_library as pl
import sys
import os
import fnmatch
import numpy as np
from netCDF4 import Dataset, num2date, netcdftime
import datetime
import scipy.stats



#year1s=[2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019]
#year2s=[2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020]

#res1
#year1s=[2010, 2011, 2012, 2013]
#year2s=[2011, 2012, 2013, 2014]

#res2
#year1s = [2014, 2015, 2016]
#year2s = [2015, 2016, 2017]

#res3
year1s = [2017, 2018, 2019]
year2s = [2018, 2019, 2020]


array_len=31

savedir = "/gws/nopw/j04/klingaman/emerton/ukmo_TIGGE_analysis/translation_speed_errors/"


if not os.path.exists(savedir):
    os.makedirs(savedir)
    
    
#function for computing the mean and the confidence intervals of the data (omitting nan values from the calculations)
def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.nanmean(a), scipy.stats.sem(a, nan_policy='omit')
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h


#Do the calculations comparing the forecasts to both the analysis track and the observed track
for obs_track in ["analysis", "ibtracs"]:#, ["analysis","ibtracs"]:

	for fcst_type in ["mean", "ctrl"]:

		all_storms_speed_err, all_storms_speed_bias, all_storms_wgt = (np.zeros(array_len) for i in range(3))
	
		for y1, y2 in zip(year1s, year2s):
			datadir="/gws/nopw/j04/klingaman/emerton/ukmo_TIGGE_analysis/reformatted_track_files_per_storm/Y"+str(y1)+str(y2)+"/SIO_storms/" #SIO_storms/"
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
			
				global all_storms_speed_err
				global all_storms_speed_bias

				if use_analysis_or_ibtracs == "analysis":
					obs_file = analysis_file
				elif use_analysis_or_ibtracs == "ibtracs":
					obs_file = ibtracs_file
				#Find out number of forecast tracks for this storm
				NT = len(nwp_files_list)
	
				#empty arrays to hold the track error and intensity bias statistics for this storm, for each forecast track
				storm_err, storm_bias =(np.zeros((NT,array_len)) for i in range(2))

				#empty arrays to hold the error sums and counts, for calculating the average errors for this storm
				storm_err_sum, storm_err_wgt, storm_bias_sum, storm_bias_wgt = (np.zeros(array_len) for i in range(4))

				#Get the date, lon, lat and vorticity data for the observed track
				obs_data=np.genfromtxt(obs_file, dtype=float, skip_header=1)
				obs_lon=obs_data[:,7]
				obs_lat=obs_data[:,8]
				obs_datelist = pl.get_dates(obs_data)
		
				#compute the statistics for each forecast of this storm, at each lead time
				"""need to make sure we compare the forecast timestamp with the correct obs timestamp!!"""
				print NT
				for ff,a in zip(nwp_files_list, range(NT)):
					#print ff

					#get the forecast date, lat, lon and vorticity
					fcst_data=np.genfromtxt(ff,dtype=float,skip_header=1)
					fcst_lon=fcst_data[:,7]
					fcst_lat=fcst_data[:,8]
					fcst_datelist = pl.get_dates(fcst_data)

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

						#get the observed track data for the dates where this forecast track matches the observed track
						#at each lead time
						print "indices_obs: ", indices_obs
						for i,z in zip(indices_obs,range(array_len-lt_to_skip)):
							matched_data_dates[z+lt_to_skip, 0] = obs_datelist[i]
							matched_obs_lon_lat[z+lt_to_skip,0] = obs_lon[i]
							matched_obs_lon_lat[z+lt_to_skip,1] = obs_lat[i]
	
						#get the forecast track data for the dates where this forecast track matches the observed track, at each lead time

						print "indices_fcst: ", indices_fcst
						print "indices_fcst[0:array_len-1]: ", indices_fcst[0:array_len-1]
						for i,z in zip(indices_fcst[0:array_len-1], range(array_len-lt_to_skip)):
							matched_data_dates[z+lt_to_skip,1] = fcst_datelist[i]
							matched_fcst_lon_lat[z+lt_to_skip, 0] = fcst_lon[i]
							matched_fcst_lon_lat[z+lt_to_skip, 1] = fcst_lat[i]
					
					
				
						#calculate the translation (propagation) speed along the whole forecast track, and along the whole observed track
						obs_speed = ts.prop_speed_vals(matched_obs_lon_lat[:,0], matched_obs_lon_lat[:,1])
						fcst_speed = ts.prop_speed_vals(matched_fcst_lon_lat[:,0], matched_fcst_lon_lat[:,1])
				
						print len(obs_speed)
						print len(fcst_speed)
				
						err, bias = (np.zeros(array_len) for i in range(2))
						for lt in range(array_len-1):
					
							#trerr[lt]=ts.trerr(matched_obs_lon_lat[lt,:],matched_fcst_lon_lat[lt,:])
							err[lt] = ts.prop_speed_abs_err(obs_speed[lt], fcst_speed[lt])
							bias[lt] = ts.prop_speed_bias(obs_speed[lt], fcst_speed[lt])
					
					
						#add the errors for this forecast track, to the arrays holding all the errors for this storm
						#and add one to the "weight" for each error, which counts the number of forecasts contributing to the error calculation
						#this is because some forecasts are shorter than others, and we want to divide by the correct sample size

				
						for lt in range(array_len):
							if not np.isnan(err[lt]):
								#storm_err[a,lt] = err[lt]
								#storm_err_sum[lt] += err[lt]
								#storm_err_wgt[lt] += 1
								#all_storms_err_sum[lt] += err[lt]
								all_storms_wgt[lt] += 1
							
							#if not np.isnan(bias[lt]):
								#storm_bias[a,lt] = bias[lt]
								#storm_bias_sum[lt] += bias[lt]
								#storm_bias_wgt[lt] += 1
								#all_storms_bias_sum[lt] += bias[lt]
								#all_storms_bias_wgt[lt] += 1
							
							
						if np.all(all_storms_speed_err == 0):
							for lt in range(array_len):
								all_storms_speed_err[lt] = err[lt]
							
								all_storms_speed_bias[lt] = bias[lt]
							
						else:
							all_storms_speed_err = np.vstack([all_storms_speed_err, err])
							all_storms_speed_bias = np.vstack([all_storms_speed_bias, bias])
				
		
				#calculate the average error at each lead time, across all the forecasts of *this* storm
				#storm_err_mean, storm_bias_mean = (np.zeros(array_len) for i in range(2))
				#for lt in range(array_len):
					#storm_err_mean[lt] = storm_err_sum[lt] / storm_err_wgt[lt]
					#storm_bias_mean[lt] = storm_bias_sum[lt] / storm_bias_wgt[lt]
			
			
				#print storm_bias_sum
				#print storm_bias_wgt
				#print storm_bias_mean

				#np.savetxt(savedir + dir + "_each_forecast_translation_speed_error_per_lead_time_vs_"+use_analysis_or_ibtracs+".txt", storm_err[:,:], '%.4f')
				#np.savetxt(savedir + dir + "_each_forecast_translation_speed_bias_per_lead_time_vs_"+use_analysis_or_ibtracs+".txt", storm_bias[:,:], '%.4f')

				#np.savetxt(savedir + dir + "_average_translation_speed_error_per_lead_time_vs_"+use_analysis_or_ibtracs+".txt", storm_err_mean[:], '%.4f')
				#np.savetxt(savedir + dir + "_average_translation_speed_bias_per_lead_time_vs_"+use_analysis_or_ibtracs+".txt", storm_bias_mean[:], '%.4f')
				
			
			### CALL THE FUNCTION ###
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
				pattern="ukmo_tigge_"+fcst_type+"*.txt"
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


		#calculate the average errors at each lead time, across all the storms in this TC season / directory
		avg_err_all_storms, avg_bias_all_storms = (np.zeros((array_len,3)) for i in range(2))
	
		for lt in range(array_len):
			avg_err_all_storms[lt,0], avg_err_all_storms[lt,1], avg_err_all_storms[lt,2] = mean_confidence_interval(all_storms_speed_err[:,lt])
			avg_bias_all_storms[lt,0], avg_bias_all_storms[lt,1], avg_bias_all_storms[lt,2] = mean_confidence_interval(all_storms_speed_bias[:,lt])
				
		
		np.savetxt(savedir + str(year1s[0])+"-"+str(year2s[-1])+".average_translation_speed_error_per_lead_time.vs_" + obs_track + "."+fcst_type+".with_confidence_intervals.txt",  avg_err_all_storms, '%.4f')
		np.savetxt(savedir + str(year1s[0])+"-"+str(year2s[-1])+".average_translation_speed_bias_per_lead_time.vs_" + obs_track + "."+fcst_type+".with_confidence_intervals.txt", avg_bias_all_storms[:], '%.4f')
	
		np.savetxt(savedir + str(year1s[0])+"-"+str(year2s[-1])+".number_of_forecasts_included.vs_" + obs_track + "."+fcst_type+".txt",  all_storms_wgt[:], '%.4f')	
		
		
		
		
		
		
		
		
		
		
		
		
		
	
