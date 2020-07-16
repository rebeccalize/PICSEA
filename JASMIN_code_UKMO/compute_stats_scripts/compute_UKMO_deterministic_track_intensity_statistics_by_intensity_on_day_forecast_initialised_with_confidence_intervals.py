import track_statistics_library as ts
import picsea_library as pl
import sys
import os
import fnmatch
import numpy as np
from netCDF4 import Dataset, num2date, netcdftime
import datetime
import scipy.stats

#This script calculates the track location errors, wind speed and mslp bias, by intensity of the storm on the day the forecast is initialised
#It does this by taking the IBTrACS wind speed on the first day of the forecast, and checking if it's in the category we're currently running
#If not, it moves on to the next forecast, but if it is in the right category, it calculates the error at each lead time and adds the errors to the array for this category
#It loops over all the seasons listed below, and at the end, computes the mean error at each lead time, and the standard deviation (confidence intervals, currently using 95% confidence)
#It saves files of the mean and the standard deviation [mean, lower, upper] in that order, averaged over all the forecasts that fall in each category

year1s=[2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019]
year2s=[2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020]

array_len=31 #number of forecast timesteps (7 days x 4 forecasts)


savedir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_and_intensity_errors/intensity_category/"


if not os.path.exists(savedir):
    os.makedirs(savedir)



#function for computing the mean and the confidence intervals of the data (omitting nan values from the calculations)
def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.nanmean(a), scipy.stats.sem(a, nan_policy='omit')
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h


	
for category in ['pre-TC','all_tropical_storms', 'intense_and_very_intense_tropical_cyclones','moderate_tropical_storm', 'strong_tropical_storm', 'tropical_cyclone', 'intense_tropical_cyclone', 'very_intense_tropical_cyclone']: #'moderate_tropical_storm', 'strong_tropical_storm', 'tropical_cyclone', 'intense_tropical_cyclone', 'very_intense_tropical_cyclone'

	if category == 'moderate_tropical_storm':
		cat_wind_min = 63
		cat_wind_max = 89
	elif category == 'strong_tropical_storm':
		cat_wind_min = 89
		cat_wind_max = 118
	elif category == 'tropical_cyclone':
		cat_wind_min = 118
		cat_wind_max = 166
	elif category == 'intense_tropical_cyclone':
		cat_wind_min = 166
		cat_wind_max = 213
	elif category == 'very_intense_tropical_cyclone':
		cat_wind_min = 213
		cat_wind_max = 1000
		
	elif category == 'pre-TC':
		cat_wind_min = 0
		cat_wind_max = 63
	elif category == 'all_tropical_storms':
		cat_wind_min = 63
		cat_wind_max=118
	elif category == 'intense_and_very_intense_tropical_cyclones':
		cat_wind_min = 166
		cat_wind_max = 1000
		
		
	print category


	#Do the calculations comparing the forecasts to both the analysis track and the observed track
	for obs_track in ["ibtracs"]:#, ["analysis","ibtracs"]:
		print obs_track

		#empty arrays to hold the average errors at each lead time, across all of the storms and seasons, for this intensity category, vs ibtracs
		this_category_trerr, this_category_mslp_bias, this_category_wind_bias, all_storms_trerr_sum, all_storms_trerr_wgt, all_storms_mslp_bias_sum, all_storms_mslp_bias_wgt, all_storms_wind_bias_sum, all_storms_wind_bias_wgt = (np.zeros(array_len) for i in range(9))

		
		for y1, y2 in zip(year1s,year2s):
		
			datadir="/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_storm/"+str(y1)+"_"+str(y2)+"/SIO_storms/"
			
			#Find out the number of storms in this TC season and create array of the IDs so we can loop over them
			season_dirs=[]
			for root,dirs,files in os.walk(datadir):
				for dir in dirs:
					season_dirs.append(dir)
			NS = len(season_dirs) #total number of storms in this season in the SIO
	
		
			#Calculate the error of each foreacst at each lead time, for each storm in this TC season / directory, and save
			#Then average across all the forecasts of each storm, and save
			def storm_stats(ibtracs_file, analysis_file, nwp_files_list, use_analysis_or_ibtracs, array_len):
		
				global this_category_trerr
				global this_category_wind_bias
				global this_category_mslp_bias
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
				obs_data=np.genfromtxt(obs_file, dtype=float, skip_header=1)
				obs_lon=obs_data[:,7]
				obs_lat=obs_data[:,8]
				obs_mslp=obs_data[:,9]
				obs_wind=obs_data[:,10]
				obs_datelist = pl.get_dates(obs_data)

				#compute the statistics for each forecast of this storm, at each lead time
				"""need to make sure we compare the forecast timestamp with the correct obs timestamp!!"""
			
				
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
				
					#if the first forecast date isn't in the obs date list, we can't use it as we don't know the intensity (also too weak to be classed on this scale using ib)
					if fcst_datelist[0] in obs_datelist:
			
					
						z = obs_datelist.index(fcst_datelist[0])

						#if the observed wind speed on the date the forecast was initialised, is in the range for the category we're currently looking at in the first loop, 
						#then include this track in the calculations for this intensity category
					
						if cat_wind_min <= ((obs_wind[z])*3.6) < cat_wind_max:
					
							print obs_wind[z]*3.6


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
								trerr, mslp_bias, wind_bias = (np.zeros(array_len) for i in range(3))
								for lt in range(array_len):

									trerr[lt]=ts.trerr(matched_obs_lon_lat[lt,:],matched_fcst_lon_lat[lt,:])
									mslp_bias[lt]=ts.bias(matched_obs_mslp[lt],matched_fcst_mslp[lt])
									wind_bias[lt]=ts.bias(matched_obs_wind[lt],matched_fcst_wind[lt])
								

								#print wind_bias
								for lt in range(array_len):
									if not np.isnan(trerr[lt]):
										
										#all_storms_trerr_sum[lt] += trerr[lt]
										all_storms_trerr_wgt[lt] += 1
									
									if not np.isnan(mslp_bias[lt]):
										
										#all_storms_mslp_bias_sum[lt] += mslp_bias[lt]
										all_storms_mslp_bias_wgt[lt] += 1
										
									if not np.isnan(wind_bias[lt]):
									
										#all_storms_wind_bias_sum[lt] += wind_bias[lt]
										all_storms_wind_bias_wgt[lt] += 1
							
								if np.all(this_category_trerr == 0.):
									for lt in range(array_len):
										this_category_trerr[lt] = trerr[lt]
										this_category_mslp_bias[lt] = mslp_bias[lt]
										this_category_wind_bias[lt] = wind_bias[lt]
							
								
							
								else:
									
									
									this_category_trerr = np.vstack([this_category_trerr, trerr])
				
									this_category_mslp_bias = np.vstack([this_category_mslp_bias, mslp_bias])
				
									this_category_wind_bias = np.vstack([this_category_wind_bias, wind_bias])
									
							
	
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


		#calculate the average errors at each lead time, across all the storms in this TC season / directory
		avg_trerr_all_storms, avg_mslp_bias_all_storms, avg_wind_bias_all_storms = (np.zeros((array_len,3)) for i in range(3))

		
		print this_category_trerr
		print np.shape(this_category_trerr)
		
		#calculate the mean across allthe forecasts of all the years we've included, for this category, alongside the confidence intervals	 
		for lt in range(array_len):
			
			#This would output the mean, the min bound for 95% confidence interval, and the max bound for 95% confidence interval, in that order
			
			avg_trerr_all_storms[lt,0], avg_trerr_all_storms[lt,1], avg_trerr_all_storms[lt,2] = mean_confidence_interval(this_category_trerr[:,lt])
			avg_mslp_bias_all_storms[lt,0], avg_mslp_bias_all_storms[lt,1], avg_mslp_bias_all_storms[lt,2] = mean_confidence_interval(this_category_mslp_bias[:,lt])
			avg_wind_bias_all_storms[lt,0], avg_wind_bias_all_storms[lt,1], avg_wind_bias_all_storms[lt,2] = mean_confidence_interval(this_category_wind_bias[:,lt])
			
		print 'mean, lower bound, upper bound, for location error using 95% confidence interval: '
		print avg_trerr_all_storms
		
		print avg_mslp_bias_all_storms
		print avg_wind_bias_all_storms
			
		print np.shape(avg_trerr_all_storms)
		
		
		

		np.savetxt(savedir + str(year1s[0])+"_"+str(year2s[-1])+"."+str(category)+".average_location_error_per_lead_time_vs_"+obs_track+".txt", avg_trerr_all_storms, '%.4f')
		np.savetxt(savedir + str(year1s[0])+"_"+str(year2s[-1])+"."+str(category)+".average_mslp_bias_per_lead_time_vs_"+obs_track+".txt", avg_mslp_bias_all_storms, '%.4f')
		np.savetxt(savedir + str(year1s[0])+"_"+str(year2s[-1])+"."+str(category)+".average_wind_bias_per_lead_time_vs_" + obs_track + ".txt", avg_wind_bias_all_storms, '%.4f')


	
		np.savetxt(savedir + str(year1s[0])+"_"+str(year2s[-1])+"."+str(category)+"_number_of_forecasts_included_vs_" + obs_track + ".txt",  all_storms_trerr_wgt[:], '%.4f')
	

	
	
	
	
