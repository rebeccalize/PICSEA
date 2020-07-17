import track_statistics_library as ts
import picsea_library as pl
import sys
import os
import fnmatch
import numpy as np

#This script, for one TC season (e.g. 2015-2016) in the SIO:
#Uses functions in track_statistics_library.py, to:
#Compute the track location error and the wind and pressure bias of ECMWF forecasts that were initiated in a certain phase of the MJO
#For the start date of each forecast, it checks the MJO amplitude of the given phase, and if greater than 1, uses this forecast in the analysis
#If <1, it moves on and doesn't include the forecast

y1=sys.argv[1]
y2=sys.argv[2]
MJO=sys.argv[3]

array_len=42 #number of forecast timesteps (7 days x 4 forecasts + analysis)

#CURRENTLY RUNNING OTHER SH STORMS, NOT SIO - CHECK BEFORE RUNNING

#FOR DIFFERENT SEASONS, or across like 2006 - 2016, use this:
datadir="/perm/mo/more/TIGGE/Y"+str(y1)+str(y2)+"/reformatted_track_files_per_storm_WITH_INTENSITY/SIO_storms/" #SIO_storms/"
savedir = "/perm/mo/more/picsea/MJO_phase_errors/"

#FOR LANDFALLING CYCLONES IN DIFFERENT COUNTRIES, USE THESE:
#datadir="/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_storm/mozambique_landfalling/" #SIO_storms/"
#savedir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_and_intensity_errors/mozambique_landfalling/"





if not os.path.exists(savedir): 
    os.makedirs(savedir)

#Find out the number of storms in this TC season and create array of the IDs so we can loop over them
season_dirs=[]
for root,dirs,files in os.walk(datadir):
	for dir in dirs:
		season_dirs.append(dir)
NS = len(season_dirs) #total number of storms in this season in the SIO




MJO_file = "/MJO_INDICES/MJO_phase"+str(MJO)+".jan-dec_dmeans_ts.1979-2019.nc"
print MJO_file

ffMJO = Dataset(MJO_file,'r')
MJOamp = ffMJO.variables['phase_ts'][:]
MJOdatesnc = ffMJO.variables['time'][:]

t_unit = ffMJO.variables['time'].units
t_cal = ffMJO.variables['time'].calendar
tvalue = num2date(MJOdatesnc,units=t_unit, calendar=t_cal)

MJOdates = [i.strftime("%Y-%m-%d") for i in tvalue]





#Do the calculations comparing the forecasts to both the analysis track and the observed track
for obs_track in ["analysis", "ibtracs"]:#, ["analysis","ibtracs"]:
	print obs_track
	
	
	for fcst_type in ["det", "mean", "ctrl"]:
		print fcst_type

		#empty arrays to hold the average errors at each lead time, across all of the storms
		all_storms_trerr_sum, all_storms_trerr_wgt, all_storms_mslp_bias_sum, all_storms_mslp_bias_wgt, all_storms_wind_bias_sum, all_storms_wind_bias_wgt = (np.zeros(array_len) for i in range(6))
	
		#Calculate the error of each foreacst at each lead time, for each storm in this TC season / directory, and save
		#Then average across all the forecasts of each storm, and save
		def storm_stats(ibtracs_file, analysis_file, nwp_files_list, use_analysis_or_ibtracs, array_len):
		
			print ibtracs_file
			print analysis_file
		
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
			print NT
			for ff,a in zip(nwp_files_list, range(NT)):
				#print ff

				#get the forecast date, lat, lon and vorticity
				fcst_data=np.genfromtxt(ff,dtype=float,skip_header=1)
				fcst_lon=fcst_data[:,7]
				fcst_lat=fcst_data[:,8]
				fcst_mslp=fcst_data[:,9]
				fcst_wind=fcst_data[:,10]
				
				fcst_datelist = pl.get_dates(fcst_data)
				
				#what's the initialisation date of this forecast?
				first_fcst_date = datetime.datetime.strptime(str(fcst_datelist[0]), "%Y%m%d%H")	
				
				print first_fcst_date.strftime("%m-%d")
				
				if first_fcst_date.strftime("%m-%d") == '02-29':
				
					continue #leap years missing in MJO data
					
				#Find the index of this date in the MJO data file
				else: z = MJOdates.index(first_fcst_date.strftime("%Y-%m-%d"))
				
				#Find the MJO amplitude on this date
				amp = MJOamp[z]
				
				#If the MJO amplitude is <0, we don't want to use this forecast, so continue to the next forecast
				if amp < 1.0:
					continue
				
				#If the amplitude is greater than 1, we want to include this forecast in the analysis for this MJO phase	
				elif amp >= 1.0:
				
				

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
							matched_obs_wind[z+lt_to_skip] = obs_wind[i]


						#get the forecast track data for the dates where this forecast track matches the observed track, at each lead time


						for i,z in zip(indices_fcst[0:array_len-1], range(array_len-lt_to_skip)):
							matched_data_dates[z+lt_to_skip,1] = fcst_datelist[i]
							matched_fcst_lon_lat[z+lt_to_skip, 0] = fcst_lon[i]
							matched_fcst_lon_lat[z+lt_to_skip, 1] = fcst_lat[i]
							matched_fcst_mslp[z+lt_to_skip] = fcst_mslp[i]
							matched_fcst_wind[z+lt_to_skip] = fcst_wind[i]

						#print matched_data_dates
					
						#calculate the translation (propagation) speed along the whole forecast track, and along the whole observed track
						obs_speed = ts.prop_speed_vals(matched_obs_lon_lat[:,0], matched_obs_lon_lat[:,1])
						fcst_speed = ts.prop_speed_vals(matched_fcst_lon_lat[:,0], matched_fcst_lon_lat[:,1])

						#calculate the track error (great circle distance, in km) for this forecast, at each lead time
						trerr, mslp_bias, wind_bias = (np.zeros(array_len-1) for i in range(3))
						for lt in range(array_len-1):

							trerr[lt] = ts.trerr(matched_obs_lon_lat[lt,:], matched_fcst_lon_lat[lt,:])
							mslp_bias[lt]=ts.bias(matched_obs_mslp[lt],matched_fcst_mslp[lt])
							wind_bias[lt]=ts.bias(matched_obs_wind[lt],matched_fcst_wind[lt])

						#add the errors for this forecast track, to the arrays holding all the errors for this storm
						#and add one to the "weight" for each error, which counts the number of forecasts contributing to the error calculation
						#this is because some forecasts are shorter than others, and we want to divide by the correct sample size

						for lt in range(array_len-1):
							if not np.isnan(trerr[lt]):
								#storm_trerr[a,lt] = trerr[lt]
								#storm_trerr_sum[lt] += trerr[lt]
								#storm_trerr_wgt[lt] += 1
								all_storms_trerr_sum[lt] += trerr[lt]
								all_storms_trerr_wgt[lt] += 1
							if not np.isnan(mslp_bias[lt]):
								#storm_mslp_bias[a,lt] = mslp_bias[lt]
								#storm_mslp_bias_sum[lt] += mslp_bias[lt]
								#storm_mslp_bias_wgt[lt] += 1
								all_storms_mslp_bias_sum[lt] += mslp_bias[lt]
								all_storms_mslp_bias_wgt[lt] += 1
							if not np.isnan(wind_bias[lt]):
								#storm_wind_bias[a,lt] = wind_bias[lt]
								#storm_wind_bias_sum[lt] += wind_bias[lt]
								#storm_wind_bias_wgt[lt] += 1
								all_storms_wind_bias_sum[lt] += wind_bias[lt]
								all_storms_wind_bias_wgt[lt] += 1

			#If you wanted to save the stats for each individual forecast of each storm, you'd do that here - but this is a repetition of data we've already saved
			#calculate the average error at each lead time, across all the forecasts of *this* storm
			storm_err_mean, storm_bias_mean = (np.zeros(array_len-1) for i in range(2))
			#for lt in range(array_len-1):
				#storm_err_mean[lt] = storm_err_sum[lt] / storm_err_wgt[lt]
				#storm_bias_mean[lt] = storm_bias_sum[lt] / storm_bias_wgt[lt]

			#print storm_bias_sum
			#print storm_bias_wgt
			#print storm_bias_mean

			
			#np.savetxt(savedir + dir + "_each_forecast_translation_speed_error_per_lead_time_vs_"+use_analysis_or_ibtracs+"_"+fcst_type+".txt", storm_err[:,:], '%.4f')
			#np.savetxt(savedir + dir + "_each_forecast_translation_speed_bias_per_lead_time_vs_" + use_analysis_or_ibtracs + "_"+fcst_type+".txt", storm_bias[:, :], '%.4f')


			#np.savetxt(savedir + dir + "_average_translation_speed_error_per_lead_time_vs_"+use_analysis_or_ibtracs+"_"+fcst_type+".txt", storm_err_mean[:], '%.4f')
			#np.savetxt(savedir + dir + "_average_translation_speed_bias_per_lead_time_vs_" + use_analysis_or_ibtracs + "_"+fcst_type+".txt", storm_bias_mean[:], '%.4f')


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
				else:
					continue
			#run the statistics for this storm
			
			if ibtracs_file==0 or analysis_file==0:
				continue
			else:
			
				storm_stats(ibtracs_file, analysis_file, nwp_files_list,obs_track,array_len)


		#calculate the average errors at each lead time, across all the storms in this TC season / directory
		avg_trerr_all_storms, avg_mslp_bias_all_storms, avg_wind_bias_all_storms = (np.zeros(array_len) for i in range(3))
		
		for lt in range(array_len-1):
			avg_trerr_all_storms[lt] = all_storms_trerr_sum[lt] / all_storms_trerr_wgt[lt]
			avg_mslp_bias_all_storms[lt] = all_storms_mslp_bias_sum[lt] / all_storms_mslp_bias_wgt[lt]
			avg_wind_bias_all_storms[lt] = all_storms_wind_bias_sum[lt] / all_storms_wind_bias_wgt[lt]

		print all_storms_wind_bias_sum
		print all_storms_wind_bias_wgt
		print avg_wind_bias_all_storms

		
		np.savetxt(savedir + str(y1)+"_"+str(y2)+"_MJO_phase"+str(MJO)+"_average_location_error_per_lead_time_vs_" + obs_track + "_"+fcst_type+".txt",  avg_trerr_all_storms[:], '%.4f')
		np.savetxt(savedir + str(y1)+"_"+str(y2)+"_MJO_phase"+str(MJO)+"_average_mslp_bias_per_lead_time_vs_" + obs_track + "_"+fcst_type+".txt", avg_mslp_bias_all_storms[:], '%.4f')
		np.savetxt(savedir + str(y1)+"_"+str(y2)+"_MJO_phase"+str(MJO)+"_average_wind_bias_per_lead_time_vs_" + obs_track + "_"+fcst_type+".txt", avg_wind_bias_all_storms[:], '%.4f')
	
		np.savetxt(savedir + str(y1)+"_"+str(y2)+"_MJO_phase"+str(MJO)+"_number_of_forecasts_included_vs_" + obs_track + "_"+fcst_type+".txt",  all_storms_trerr_wgt[:], '%.4f')
		
	
	
	
