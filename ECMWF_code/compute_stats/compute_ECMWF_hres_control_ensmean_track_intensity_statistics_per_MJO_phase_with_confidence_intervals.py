import track_statistics_library as ts
import picsea_library as pl
import sys
import os
import fnmatch
import numpy as np
from netCDF4 import Dataset, num2date
import datetime
import scipy.stats

#This script, for one TC season (e.g. 2015-2016) in the SIO:
#Uses functions in track_statistics_library.py, to:
#Compute the track location error for every consecutive UKMO control and ensemble mean forecast of each storm (great circle distance in km), at each lead time
#Compute the intensity bias (MSLP and 10m wind speed) for every consecutive UKMO control and ensemble mean forecast of each storm, at each lead time
#Compute the average track location error and intensity bias across all the forecasts for each storm
#Compute the average track location error and intensity bias across all the forecasts and across all the storms in this TC season

year1s=[2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019] #2017, 2019
year2s=[2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020] #2018, 2020


MJO=int(sys.argv[1]) #MJO phase 



#function for computing the mean and the confidence intervals of the data (omitting nan values from the calculations)
def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.nanmean(a), scipy.stats.sem(a, nan_policy='omit')
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h

array_len=30 #number of forecast timesteps (7 days x 4 forecasts)

#CURRENTLY RUNNING OTHER SH STORMS, NOT SIO - CHECK BEFORE RUNNING

#FOR DIFFERENT SEASONS, or across like 2006 - 2016, use this:

savedir = "/perm/mo/more/picsea/ecmwf_track_and_intensity_errors/MJO_phase/"


if not os.path.exists(savedir): 
    os.makedirs(savedir)


#With the new MJO files that include up to May 2020, don't have individual phase files
#So may need to add an extra bit to the if statement where it checks the amplitude, to also check the phase. Easy enough, just don't forget!!

MJO_file = "/home/mo/more/PICSEA/MJO_INDICES/MJO_indices_rmm1_rmm2.jan1979-may2020_dmeans_ts.index_values.nc"
print MJO_file

print MJO_file

ffMJO = Dataset(MJO_file,'r')
MJOamp = ffMJO.variables['amplitude'][:]
MJOphase = ffMJO.variables['phase'][:]


start_time = datetime.datetime(1979,1,1)
MJO_tvalue = np.array([start_time + datetime.timedelta(days=i) for i in xrange(len(MJOamp))])
			
MJOdates = [j.strftime("%Y-%m-%d") for j in MJO_tvalue]

ffMJO.close()



#Do the calculations comparing the forecasts to both the analysis track and the observed track
for obs_track in ["ibtracs"]:#, ["analysis","ibtracs"]:
	#print obs_track
	
	for fcst_type in ["det"]:
	
		#empty arrays to hold the average errors at each lead time, across all of the storms
		this_mjo_trerr, all_storms_trerr_wgt, this_mjo_mslp_bias, this_mjo_wind_bias = (np.zeros(array_len) for i in range(4))
	
	
		for y1, y2 in zip(year1s, year2s):
		
			datadir="/perm/mo/more/TIGGE/Y"+str(y1)+str(y2)+"/reformatted_track_files_per_storm_correct/SIO_storms/"
			#print datadir
			#Find out the number of storms in this TC season and create array of the IDs so we can loop over them
			season_dirs=[]
			for root,dirs,files in os.walk(datadir):
				for dir in dirs:
					season_dirs.append(dir)
			NS = len(season_dirs) #total number of storms in this season in the SIO
			
			#print fcst_type
			#print NS

			

			#Calculate the error of each foreacst at each lead time, for each storm in this TC season / directory, and save
			#Then average across all the forecasts of each storm, and save
			def storm_stats(ibtracs_file, analysis_file, nwp_files_list, use_analysis_or_ibtracs, array_len):
				"""This function calculates the statistics for all the forecasts of one storm"""
				
				global this_mjo_trerr
				global this_mjo_wind_bias
				global this_mjo_mslp_bias

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
				
					first_fcst_date = datetime.datetime.strptime(str(fcst_datelist[0]), "%Y%m%d%H")	
				
					#print first_fcst_date.strftime("%m-%d")
					
					if first_fcst_date.strftime("%m-%d") == '02-29':
				
						continue #leap years missing in MJO data
					
					else: 
					
						z = MJOdates.index(first_fcst_date.strftime("%Y-%m-%d"))
				
						amp = MJOamp[z]
				
					
						if MJOphase[z] == MJO and MJOamp[z] >= 1.0: 
				
				
	
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

								#add the errors for this forecast track, to the arrays holding all the errors for this storm
								#and add one to the "weight" for each error, which counts the number of forecasts contributing to the error calculation
								#this is because some forecasts are shorter than others, and we want to divide by the correct sample size
	
								#print wind_bias
								for lt in range(array_len):
									if not np.isnan(trerr[lt]):
							
										all_storms_trerr_wgt[lt] += 1
								
								if np.all(this_mjo_trerr == 0):
									for lt in range(array_len):
										this_mjo_trerr[lt] = trerr[lt]
										this_mjo_mslp_bias[lt] = mslp_bias[lt]
										this_mjo_wind_bias[lt] = wind_bias[lt]
									
								else:
									this_mjo_trerr = np.vstack([this_mjo_trerr, trerr])
									this_mjo_mslp_bias = np.vstack([this_mjo_mslp_bias, mslp_bias])
									this_mjo_wind_bias = np.vstack([this_mjo_wind_bias, wind_bias])
		

						else:
						
							continue
							
							
			### CALL THE FUNCTION ###
			#loop over storm directories in this TC season, and in each directory find each file that matches the analysis, ibtracs and nwp filenames
			#then call the function for each storm using the identified files
			for dir, x in zip(season_dirs,range(len(season_dirs))):
				#print datadir+dir
				nwp_files_list = []
				analysis_file=0
				ibtracs_file=0
				#make a list of all the files in this directory
				list_of_all_files = os.listdir(datadir+dir)
				#print list_of_all_files
				pattern="ecmwf_"+fcst_type+"*.txt"
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
		avg_trerr_all_storms, avg_mslp_bias_all_storms, avg_wind_bias_all_storms = (np.zeros((array_len, 3)) for i in range(3))
		
		for lt in range(array_len):
			avg_trerr_all_storms[lt,0], avg_trerr_all_storms[lt,1], avg_trerr_all_storms[lt,2] = mean_confidence_interval(this_mjo_trerr[:,lt])
			avg_mslp_bias_all_storms[lt,0], avg_mslp_bias_all_storms[lt,1], avg_mslp_bias_all_storms[lt,2] = mean_confidence_interval(this_mjo_mslp_bias[:,lt])
			avg_wind_bias_all_storms[lt,0], avg_wind_bias_all_storms[lt,1], avg_wind_bias_all_storms[lt,2] = mean_confidence_interval(this_mjo_wind_bias[:,lt])



		np.savetxt(savedir + str(year1s[0])+"_"+str(year2s[-1])+"_MJO_phase"+str(MJO)+"_average_location_error_per_lead_time_vs_" + obs_track + "_"+fcst_type+".with_confidence_intervals_95.txt",  avg_trerr_all_storms, '%.4f')
		np.savetxt(savedir + str(year1s[0])+"_"+str(year2s[-1])+"_MJO_phase"+str(MJO)+"_average_mslp_bias_per_lead_time_vs_" + obs_track + "_"+fcst_type+".with_confidence_intervals_95.txt", avg_mslp_bias_all_storms, '%.4f')
		np.savetxt(savedir + str(year1s[0])+"_"+str(year2s[-1])+"_MJO_phase"+str(MJO)+"_average_wind_bias_per_lead_time_vs_" + obs_track + "_"+fcst_type+".with_confidence_intervals_95.txt", avg_wind_bias_all_storms, '%.4f')
	
		np.savetxt(savedir + str(year1s[0])+"_"+str(year2s[-1])+"_MJO_phase"+str(MJO)+"_number_of_forecasts_included_vs_" + obs_track + "_"+fcst_type+".txt",  all_storms_trerr_wgt, '%.4f')

	
	
	
	
